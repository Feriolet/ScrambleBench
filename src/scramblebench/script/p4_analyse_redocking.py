import yaml
import argparse
import logging
import sys
import subprocess
import os
import tempfile
import rdkit
import pandas as pd

from rdkit import Chem
from typing import Any, Callable
from pathlib import Path
from copy import deepcopy

from scramblebench.script.config_preparation import config_constant, config_input, config_redocking, config_parameter
from scramblebench.script.utils.error_handler import DirNotFoundError
from scramblebench.script.utils.process_data import read_input, fetch_model_file_from_model_dir
from scramblebench.script.utils.docking_utils.prepare_protein import VinaProtein, GlideProtein
from scramblebench.script.utils.docking_utils.prepare_docking import calculate_easydock_cpu, fetch_glide_sp_cpu, fetch_ligprep_cpu, prepare_easydock_grid,\
                                                                     prepare_easydock_grid, prepare_easydock_vina_parameter, prepare_ligprep_inp_file,\
                                                                     prepare_glide_inp_file, prepare_easydock_ligand_input, process_easydock_docking_output


logger = logging.getLogger(__name__)


def run_ligprep_protonation(schrodinger_dir, output_dir, valid_molecule_file_dict):

    cmd = [f'{schrodinger_dir}/ligprep']

    current_dir = Path.cwd() 
    with tempfile.TemporaryDirectory() as temp_protonation_dir:      
        os.chdir(temp_protonation_dir)
        for model, fname in valid_molecule_file_dict.items():
            model_cmd = deepcopy(cmd)
            
            ligprep_inp = str(Path(temp_protonation_dir) / f'{Path(fname).stem}_ligprep.inp')
            ligprep_output_dir = Path(Path(output_dir) / model)
            ligprep_output_dir.mkdir(parents=True, exist_ok=True)
            ligprep_output_fname = ligprep_output_dir / f'{Path(fname).stem}_protonated.sdf'

            if Path(ligprep_output_fname).is_file():
                logging.info(f'Ligprep protonation with {fname} has been done. Skipping...')
                continue

            subprocess.run(['cp', fname, temp_protonation_dir], text=True)
            prepare_ligprep_inp_file(ligand_fname=fname,
                                     output_fname=ligprep_output_fname,
                                     write_fname=ligprep_inp)
            
            protonation_cpu = str(fetch_ligprep_cpu(schrodinger_dir=schrodinger_dir))
            model_cmd += ['-inp', ligprep_inp, 
                          '-NJOBS', protonation_cpu, '-JOBNAME', 
                          f'ligprep_{Path(temp_protonation_dir).name}', 
                          '-HOST', f'localhost:{protonation_cpu}', '-WAIT']

            logging.info(f'Ligprep is running with input: {fname}, output: {ligprep_output_fname}, cmd: {model_cmd}')
            subprocess.run(model_cmd, text=True)
            subprocess.run(['cp', Path(ligprep_output_fname).name, ligprep_output_fname], text=True)

    os.chdir(current_dir)
    
    return ligprep_output_fname

def run_easydock_protonation(protonation_data, valid_molecule_file_dict):

    cmd = ['conda', 'run', '-n', protonation_data.environment_value,
            'easydock', '-c', str(calculate_easydock_cpu())]

    for model, fname in valid_molecule_file_dict.items():

        protonation_output_dir = Path(Path(protonation_data.output_value) / model)
        protonation_output_dir.mkdir(parents=True, exist_ok=True)
        output_easydock_sdf = str(f'{protonation_output_dir / Path(fname).stem}_protonated.sdf')
        if Path(output_easydock_sdf).is_file():
            logging.info(f'Protonation with {output_easydock_sdf} has been done. Skipping...')
            continue
        
        with tempfile.TemporaryDirectory() as temp_protonation_dir:      
            output_easydock_db = str(Path(temp_protonation_dir) / f'{Path(fname).stem}_protonated.db')
            model_cmd = deepcopy(cmd)                               
            model_cmd += ['-i', fname, '-o', output_easydock_db, '--protonation', protonation_data.method_value]

            logging.info(f'Running easydock protonation with cmd: {model_cmd=}')
            subprocess.run(model_cmd, text=True)

            Path(Path(protonation_data.output_value) / model).mkdir(parents=True, exist_ok=True)
            FETCH_SDF_SCRIPT = Path(__file__).parent / 'utils/docking_utils' / 'easydock_fetch_sdf_from_db.py'
            protonated_cmd = ['conda', 'run', '-n', protonation_data.environment_value,
                            'python', FETCH_SDF_SCRIPT, '-i', output_easydock_db, '-o', output_easydock_sdf]

            subprocess.run(protonated_cmd, text=True)


def run_protonation(config_data):
    redocking_data = config_redocking.RedockingConfig(config_data[config_constant.ANALYSIS_KEY])
    parameter_data = config_parameter.ParameterConfig(config_data=config_data)
    redocking_data.validate_config()
    protonation_data = redocking_data.protonation_value
    valid_molecule_file_dict = fetch_model_file_from_model_dir(dir_path=protonation_data.input_value,
                                                               model_list=parameter_data.model_list_value)

    if 'schrodinger' in str(Path(protonation_data.environment_value).name).lower():
        run_ligprep_protonation(schrodinger_dir=protonation_data.environment_value,
                                output_dir=protonation_data.output_value,
                                valid_molecule_file_dict=valid_molecule_file_dict)
    else:
        run_easydock_protonation(protonation_data=protonation_data,
                                 valid_molecule_file_dict=valid_molecule_file_dict)


def prepare_easydock_files(docking_data: config_redocking.EasyDockConfig, input_data: config_input.InputStructure, output_dir):
    config_fname = docking_data.config_value
    protein_pdb = input_data.pdb_value
    grid_fname = str(Path(protein_pdb).parent / f'{input_data.protein_value}_grid.txt')
    
    prepare_easydock_grid(pocket_coordinates=input_data.pocket_coord_value,
                          write_fname=grid_fname)

    with open(config_fname, 'r') as config_f:
        config_data = yaml.load(config_f, Loader=yaml.SafeLoader)
    
    config_data['protein'] = VinaProtein(pdb_filepath=protein_pdb, 
                                         prepare_receptor_bin_path=docking_data.protein_pdbqt_executable_value,
                                         protonation_method=docking_data.protein_preparation_value,
                                         preparation_method=docking_data.protein_pdbqt_preparation_value).pdbqt_filepath 
    
    config_data['protein_setup'] = grid_fname

    if str(docking_data.docking_value).lower().strip() == 'vina':
        config_data = prepare_easydock_vina_parameter(config_data=config_data)
    else:
        logging.warning(f'You are using a different docking program other than vina. Ensure that the config at {new_config_fname} has \
                        the correct format. Please check easydock github for reference!')

    new_config_fname = Path(output_dir) / f'easydock_{Path(input_data.protein_value).stem}.yml'
    with open(new_config_fname, 'w') as config_f:
        yaml.dump(config_data, config_f)

    return new_config_fname


def run_easydock_plif( environment, db_fname, protein_fname, output_fname, 
                      reference_fname, ncpu=1, reference_plif=None, similarity=0.8):
    
    if not Path(output_fname).parent.is_dir():
        Path(output_fname).parent.mkdir(parents=True, exist_ok=True)

    if not reference_plif:
    
        reference_output = Path(db_fname).parent / 'ref_similarity.csv'
        subprocess.run(['conda', 'run', '-n', environment,
                        'easydock_plif', '-i', reference_fname, '-p', protein_fname,
                        '-o', reference_output, '-c', str(4)])

        reference_plif = [col for col in pd.read_csv(reference_output, delimiter='\t').columns.tolist() if col != 'Name']
    
    plif_result = Path(db_fname).parent / 'similarity.csv'
    cmd = ['conda', 'run', '-n', environment,
           'easydock_plif', '-i', db_fname, '-p', protein_fname, '-o',
            plif_result, '--ref_plif']
    
    cmd += reference_plif
    cmd += ['-c', str(ncpu)]

    logging.info(f'Running easydock PLIF with cmd: {cmd}')
    subprocess.run(cmd)
    
    plif_df = pd.read_csv(plif_result, delimiter='\t').set_index('id')

    filtered_plif_df = plif_df[plif_df['plif_sim'] >= 0.2]

    total_mol = []
    for pose, df in filtered_plif_df.groupby('pose'):
        output_pose_sdf = Path(db_fname).parent / f'filtered_sdf_{pose}.sdf'
        mol_name_list = df.index.tolist() # assuming only 1 stereoisomer

        cmd = ['conda', 'run', '-n', environment,
               'get_sdf_from_easydock', '-i', db_fname]
        cmd += ['-d'] + mol_name_list + ['--poses', str(pose)]
        cmd += ['-o', output_pose_sdf, '--fields', 'docking_score']
        
        subprocess.run(cmd)

        total_mol += [mol for mol in Chem.SDMolSupplier(output_pose_sdf, removeHs=False)]


    if len(total_mol) > 0:
        temp_output = Path(db_fname).parent / f'filtered_sdf_{pose}.sdf'
        with Chem.SDWriter(temp_output) as output_f:
            for mol in total_mol:
                similarity_score = filtered_plif_df.loc[mol.GetProp('_Name').rsplit(':', 1)[0], 'plif_sim']
                mol.SetProp('plif_similarity', str(similarity_score))
                output_f.write(mol)
        
        processed_temp_output = process_easydock_docking_output(temp_output)

        subprocess.run(['cp', processed_temp_output, output_fname], text=True)
    else:
        logging.warning('No molecule matching with PLIF')


def run_easydock_docking(docking_data: config_redocking.EasyDockConfig, parameter_data, input_data):
    
    cmd = ['conda', 'run', '-n', docking_data.environment_value,
            'easydock',  '--sdf']

    Path(docking_data.output_value).mkdir(parents=True, exist_ok=True)

    config_fname = prepare_easydock_files(docking_data, input_data, output_dir=docking_data.output_value)
    cmd += ['--config', str(config_fname)]

    valid_molecule_file_dict = fetch_model_file_from_model_dir(dir_path=docking_data.input_value,
                                                               model_list=parameter_data.model_list_value)
    
    if docking_data.docking_value:
        cmd += ['--program', docking_data.docking_value]
    else:
        logging.warning('You did not specify which docking program to use for easydock')

    if docking_data.cpu_value:
        if hasattr(os, 'sched_getaffinity'):
            num_cpus = len(os.sched_getaffinity(0))
        else:
            num_cpus = os.cpu_count() or 1
        cmd += ['-c', str(min(int(docking_data.cpu_value), num_cpus - 5))]

    else:
        cmd += ['-c', str(calculate_easydock_cpu())]

    for model, fname in valid_molecule_file_dict.items():
        if model == 'DiffSBDD':
            continue
        output_easydock_sdf =  str(Path(docking_data.output_value) / model / f'{Path(fname).stem}_easydock_redocked.sdf')
        Path(output_easydock_sdf).parent.mkdir(parents=True, exist_ok=True)
        
        if Path(output_easydock_sdf).is_file():
            if docking_data.plif_value:
                output_easydock_plif =  str(Path(docking_data.output_value) / model / 'plif' / f'{Path(fname).stem}_easydock_plif_processed.sdf')
                if Path(output_easydock_plif).is_file():
                    logging.info(f'Easydock Docking with {output_easydock_sdf} has been done. Skipping...')
                    continue

            else:
                logging.info(f'Easydock Docking with {output_easydock_sdf} has been done. Skipping...')
                continue

        with tempfile.TemporaryDirectory() as temp_docking_dir:      
            temp_output_easydock_db = str(Path(temp_docking_dir) / f'{Path(fname).stem}_redocked.db')
            temp_output_easydock_sdf = str(Path(temp_docking_dir) / f'{Path(fname).stem}_redocked.sdf')
            
            model_cmd = deepcopy(cmd)  
            
            preprocessed_fname = deepcopy(fname) 
            if docking_data.protonation_value:
                if 'schrodinger' not in docking_data.protonation_value:
                    model_cmd += ['--protonation', docking_data.protonation_value]            
                else:
                    preprocessed_fname = run_ligprep_protonation(schrodinger_dir=docking_data.protonation_value, 
                                            output_dir=temp_docking_dir, 
                                            valid_molecule_file_dict={model: fname})   

            easydock_input_fname = prepare_easydock_ligand_input(input_fname=preprocessed_fname,
                                                                 output_dir=temp_docking_dir)

            model_cmd += ['-i', str(easydock_input_fname), '-o', temp_output_easydock_db]

            logging.info(f'Easydock docking with cmd: {model_cmd=}')
            subprocess.run(model_cmd, text=True)

            post_processed_docking_fname = process_easydock_docking_output(temp_output_easydock_sdf)
            subprocess.run(['cp', post_processed_docking_fname, output_easydock_sdf], text=True)

            if docking_data.plif_value:
                run_easydock_plif(environment=docking_data.environment_value,
                                  db_fname=temp_output_easydock_db,
                                  protein_fname=input_data.pdb_value,
                                  output_fname=output_easydock_plif,
                                  reference_fname=input_data.sdf_value,
                                  similarity=docking_data.plif_similarity_value
                                  )

def prepare_hypothesis(schrodinger_env, rec_file, lig_file, center_coord):

    current_tmpdir = Path.cwd()
    with tempfile.TemporaryDirectory() as tempdir:
        os.chdir(tempdir)
        rec_tempfile = str(Path(tempdir) / 'rec_file.pdb')
        lig_tempfile = str(Path(tempdir) / 'lig_file.sdf')

        subprocess.run(['cp', rec_file, rec_tempfile])
        subprocess.run(['cp', lig_file, lig_tempfile])

        jobname = f'{Path(rec_file).stem}_phypo'

        cmd = [f'{schrodinger_env}/utilities/epharmacophores', '-rec_file', 
            rec_tempfile, '-lig_file', lig_tempfile, 
            f'-site_center={center_coord}', '-in_place', '-j', jobname, 
            '-f', '7', '-site_dist', '2.0', '-pair_dist', '4.0', '-xvol', '-scale', 
            '0.5', '-buff', '2.0', '-limit', '5.0', '-HOST', 'localhost', '-WAIT']
    
        subprocess.run(cmd)
        
        phypo_file = f'{jobname}.phypo'

        output_fname = Path(rec_file).parent / phypo_file
        subprocess.run(['cp', phypo_file, output_fname])

    os.chdir(current_tmpdir)

    return str(output_fname)



def run_schrodinger_phase_plif(working_dir, environment, protein_fname, library_fname, output_fname,
                      reference_fname, center_coord, hypo_input=None):

    if not Path(output_fname).parent.is_dir():
        Path(output_fname).parent.mkdir(parents=True, exist_ok=True)

    if not hypo_input or isinstance(hypo_input, bool):
        hypo_input = prepare_hypothesis(schrodinger_env=environment,
                        rec_file=protein_fname,
                        lig_file=reference_fname,
                        center_coord=center_coord)
    elif not isinstance(hypo_input, str):
        raise ValueError(f'hypothesis for Phase screening expecting bool (generate own file), user phypo zip input, not {type(hypo_input)}')
    
    jobname = f'phase_screen_1'
    temp_output = f'{Path(jobname).stem}-hits.sdfgz'


    subprocess.run(['cp', hypo_input, 'hypo.phypo'])

    cmd = [f'{environment}/phase_screen', library_fname, 'hypo.phypo', 
           jobname, '-inplace', '-keep', '999999999', '-report', '1', 
           '-ex', '-HOST', 'localhost:10', '-TMPLAUNCHDIR', '-osd', '-WAIT']

    subprocess.run(cmd)

    subprocess.run(['cp', temp_output, f'{output_fname}.gz'])
    subprocess.run(['gunzip', f'{output_fname}.gz'], text=True)

def run_glide_docking(docking_data: config_redocking.GlideConfig, parameter_data, input_data):

    grid_filepath = GlideProtein(pdb_filepath=input_data.pdb_value,
                                native_ligand=list(Chem.SDMolSupplier(input_data.sdf_value, removeHs=False))[0],
                                grid_output_dirpath=str(Path(input_data.pdb_value).parent),
                                schrodinger_dirpath=docking_data.dir_value,
                                protein_preparation=docking_data.protein_preparation_value).grid_filepath

    valid_molecule_file_dict = fetch_model_file_from_model_dir(dir_path=docking_data.input_value,
                                                               model_list=parameter_data.model_list_value)

    current_dir = Path.cwd()
    for model, fname in valid_molecule_file_dict.items():

        if model == 'DiffSBDD':
            continue
        if docking_data.protonation_value:
            output_glide_fname = Path(docking_data.output_value) / model / f'{Path(fname).stem}_protonated_glide.sdf.gz'
        else:
            output_glide_fname = Path(docking_data.output_value) / model / f'{Path(fname).stem}_glide.sdf.gz'

        Path(output_glide_fname).parent.mkdir(parents=True, exist_ok=True)

        if Path(Path(output_glide_fname).parent / output_glide_fname.stem).is_file():

            if docking_data.plif_value:
                output_phase_plif = str(Path(docking_data.output_value) / model / 'plif' / f'{Path(fname).stem}_phase_plif.sdf')
                if Path(output_phase_plif).is_file():
                    logging.info(f'Phase Pharmacophore  with {output_phase_plif} has been done. Skipping...')
                    continue
            
            else:
                logging.info(f'Glide docking with {Path(fname).name} has been done. Skipping...')
                continue

        with tempfile.TemporaryDirectory() as tempdir:  
            os.chdir(tempdir)

            if docking_data.protonation_value:  
                run_ligprep_protonation(schrodinger_dir=docking_data.dir_value, 
                                        output_dir=tempdir, 
                                        valid_molecule_file_dict={model: fname})
                protonated_fname = Path(tempdir) / model / f'{Path(fname).stem}_protonated.sdf'
            else:
                protonated_fname = fname

            os.chdir(Path(protonated_fname).parent)
            jobname = f'glide_{Path(fname).stem}'
            inp_fname = Path(tempdir) / f'{jobname}.inp'
            prepare_glide_inp_file(grid_fname=grid_filepath,
                                    ligand_fname=protonated_fname,
                                    intra_hbonds=docking_data.reward_intra_hbonds_value,
                                    write_fname=inp_fname)

            cmd = [f"{docking_data.dir_value}/glide", inp_fname, '-OVERWRITE', 
                '-adjust', '-HOST', f'localhost:{fetch_glide_sp_cpu(docking_data.dir_value)}', '-TMPLAUNCHDIR', '-WAIT']

            temp_glide_sdf_output = f'{jobname}_lib.sdfgz' 
            logging.info(f'Glide docking with cmd: {cmd=}')
            subprocess.run(cmd, text=True)
            subprocess.run(['cp', temp_glide_sdf_output, output_glide_fname], text=True)
            subprocess.run(['gunzip', output_glide_fname], text=True)

            if docking_data.plif_value:
                run_schrodinger_phase_plif(tempdir, docking_data.dir_value, protein_fname=input_data.pdb_value, output_fname=output_phase_plif, library_fname=temp_glide_sdf_output,
                                           reference_fname=input_data.sdf_value, center_coord=input_data.pocket_coord_value, hypo_input=docking_data.plif_value)
    os.chdir(current_dir)           


def run_docking(config_data):
    redocking_data = config_redocking.RedockingConfig(config_data[config_constant.ANALYSIS_KEY])
    parameter_data = config_parameter.ParameterConfig(config_data=config_data)
    input_data = config_input.InputStructure(config_data[config_constant.INPUT_KEY])
    docking_data_list = redocking_data.docking_value.valid_key_list

    for docking_data in docking_data_list:
        if isinstance(docking_data, config_redocking.EasyDockConfig):
            run_easydock_docking(docking_data, parameter_data=parameter_data, input_data=input_data)
        elif isinstance(docking_data, config_redocking.GlideConfig):
            run_glide_docking(docking_data, parameter_data=parameter_data, input_data=input_data)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Redocking analysis after molecule preparation")

    parser.add_argument("-i", "--input", help="config yaml input file or txt file containing yaml filepath", required=True, type=str)
    args = parser.parse_args()

    logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='%(asctime)s - %(module)s: - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
    
    logging.info('Running p4_analyse_redocking.py')
    logging.info('Reading the config filename :\)')


    yaml_file_list = read_input(args.input)

    for yaml_file in yaml_file_list:
        checkpoint_json = Path(yaml_file).parent
        data_input = yaml.safe_load(open(yaml_file, 'r'))

        if config_constant.ANALYSIS_REDOCKING_PROTONATION_KEY in data_input[config_constant.ANALYSIS_KEY][config_constant.ANALYSIS_REDOCKING_KEY]:
            run_protonation(config_data=data_input)
        if config_constant.ANALYSIS_REDOCKING_DOCKING_KEY in data_input[config_constant.ANALYSIS_KEY][config_constant.ANALYSIS_REDOCKING_KEY]:
            run_docking(config_data=data_input)