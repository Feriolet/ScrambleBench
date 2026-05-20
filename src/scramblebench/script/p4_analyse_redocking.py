from typing import Any, Callable
from pathlib import Path

import yaml
import argparse
import logging
import sys
import subprocess
import os

from scramblebench.script.config_preparation import config_constant, config_input, config_redocking, config_parameter
from scramblebench.script.utils.error_handler import DirNotFoundError
from scramblebench.script.utils.process_data import read_input, fetch_model_file_from_model_dir
from scramblebench.script.utils.docking_utils.prepare_protein import VinaProtein, GlideProtein
from scramblebench.script.utils.docking_utils.prepare_docking import calculate_easydock_cpu, fetch_glide_sp_cpu, fetch_ligprep_cpu, prepare_easydock_grid,\
                                                                     prepare_easydock_grid, prepare_easydock_vina_parameter, prepare_ligprep_inp_file,\
                                                                     prepare_glide_inp_file

from copy import deepcopy
import tempfile

import rdkit
from rdkit import Chem

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
                                         prepare_receptor_bin_path=docking_data.protein_preparation_executable_value,
                                         protonation_method=docking_data.protein_protonation_value,
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


def run_easydock_docking(docking_data: config_redocking.EasyDockConfig, parameter_data, input_data):
    
    cmd = ['conda', 'run', '-n', docking_data.environment_value,
            'easydock',  '--sdf']

    Path(docking_data.output_value).mkdir(parents=True, exist_ok=True)

    config_fname = prepare_easydock_files(docking_data, input_data, output_dir=docking_data.output_value)
    cmd += ['--config', config_fname]

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

        output_easydock_sdf =  str(Path(docking_data.output_value) / model / f'{Path(fname).stem}_easydock_redocked.sdf')
        Path(output_easydock_sdf).parent.mkdir(parents=True, exist_ok=True)
        if Path(output_easydock_sdf).is_file():
            logging.info(f'Easydock Docking with {output_easydock_sdf} has been done. Skipping...')
            continue

        with tempfile.TemporaryDirectory() as temp_docking_dir:      
            temp_output_easydock_db = str(Path(temp_docking_dir) / f'{Path(fname).stem}_redocked.db')
            temp_output_easydock_sdf = str(Path(temp_docking_dir) / f'{Path(fname).stem}_redocked.sdf')
            
            model_cmd = deepcopy(cmd)  
            
            if docking_data.protonation_value and docking_data.protonation_value != 'ligprep':
                model_cmd += ['--protonation', docking_data.protonation_value]               

            model_cmd += ['-i', fname, '-o', temp_output_easydock_db]

            logging.info(f'Easydock docking with cmd: {model_cmd=}')
            subprocess.run(model_cmd, text=True)
            subprocess.run(['cp', temp_output_easydock_sdf, output_easydock_sdf], text=True)


def run_glide_docking(docking_data: config_redocking.GlideConfig, parameter_data, input_data):

    grid_filepath = GlideProtein(pdb_filepath=input_data.pdb_value,
                                native_ligand=list(Chem.SDMolSupplier(input_data.sdf_value))[0],
                                grid_output_dirpath=str(Path(input_data.pdb_value).parent),
                                schrodinger_dirpath=docking_data.dir_value,
                                protein_preparation=docking_data.protein_preparation_value).grid_filepath

    valid_molecule_file_dict = fetch_model_file_from_model_dir(dir_path=docking_data.input_value,
                                                               model_list=parameter_data.model_list_value)

    current_dir = Path.cwd()
    for model, fname in valid_molecule_file_dict.items():
        if docking_data.protonation_value:
            output_glide_fname = Path(docking_data.output_value) / model / f'{Path(fname).stem}_protonated_glide.sdf.gz'
        else:
            output_glide_fname = Path(docking_data.output_value) / model / f'{Path(fname).stem}_glide.sdf.gz'

        Path(output_glide_fname).parent.mkdir(parents=True, exist_ok=True)
        if Path(Path(output_glide_fname).parent / output_glide_fname.stem).is_file():
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