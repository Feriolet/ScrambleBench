from typing import Any, Callable
from pathlib import Path

import yaml
import argparse
import logging
import sys
import subprocess
import os

from scramblebench.script.config_preparation import config_constant, config_genbench3d, config_input, config_redocking, config_parameter
from scramblebench.script.utils.error_handler import DirNotFoundError
from scramblebench.script.utils.process_data import read_input, fetch_model_folder_name, find_file_name_through_regex
from scramblebench.script.utils.prepare_protein import VinaProtein, GlideProtein


from copy import deepcopy
import tempfile

import rdkit
from rdkit import Chem

logger = logging.getLogger(__name__)


def fetch_valid_prepared_molecule_file(dir_path, parameter_class: config_parameter.ParameterConfig) -> dict[str, str]:
    model_for_generation_list = parameter_class.model_list_value

    generation_folder_dirpath = Path(dir_path)
    if not generation_folder_dirpath.is_dir():
        logging.exception('You have prepared your molecule yet!')
        raise DirNotFoundError(f'{generation_folder_dirpath} is not found! Please make sure you run p3_prepare_molecule.py')
    
    valid_molecule_file_dict = {}
    for model in model_for_generation_list:

        matched_fname_list = find_file_name_through_regex(character=model, file_format='.sdf', dirname=Path(generation_folder_dirpath) / model)
        if len(matched_fname_list) > 1:
            logging.exception(f'We found more than 1 matching file for {model} model: {matched_fname_list}. Please ensure only 1 is detected')
            raise ValueError(f'We found more than 1 matching file for {model} model: {matched_fname_list}. Please ensure only 1 is detected')
        elif len(matched_fname_list) == 0:
            logging.warning(f'There are no matched file for {model} model in {generation_folder_dirpath}. Make sure this is intended')
        else:
            valid_molecule_file_dict[model] = str(matched_fname_list[0])
    
    return valid_molecule_file_dict


def run_ligprep_protonation(schrodinger_dir, output_dir, valid_molecule_file_dict):

    cmd = [f'{schrodinger_dir}/ligprep']
    if not Path(output_dir).is_dir():
        Path(output_dir).mkdir(parents=True, exist_ok=True)

    current_dir = Path.cwd()
    
    with tempfile.TemporaryDirectory() as tempfile_dir:      
        os.chdir(tempfile_dir)
        for model, fname in valid_molecule_file_dict.items():
            model_cmd = deepcopy(cmd)
            
            ligprep_inp = str(Path(tempfile_dir) / f'{Path(fname).stem}_ligprep.inp')
            ligprep_output_dir = Path(Path(output_dir) / model)
            ligprep_output_dir.mkdir(parents=True, exist_ok=True)
            ligprep_output_fname = ligprep_output_dir / f'{Path(fname).stem}_protonated.sdf'

            if Path(ligprep_output_fname).is_file():
                logging.info(f'Ligprep protonation with {fname} has been done. Skipping...')
                continue

            subprocess.run(['cp', fname, tempfile_dir], text=True)
            ligprep_inp_data = f'''
            INPUT_FILE_NAME   {Path(fname).name}
            MAX_ATOMS   500
            FORCE_FIELD   16
            PH   7.4
            PH_THRESHOLD   1.0
            EPIK   yes
            EPIKX   no
            EPIK_METAL_BINDING   no
            INCLUDE_ORIGINAL_STATE   no
            DETERMINE_CHIRALITIES   no
            IGNORE_CHIRALITIES   no
            NUM_STEREOISOMERS   32
            OUT_SD  {Path(ligprep_output_fname).name}

            '''

            with open(ligprep_inp, 'w') as ligprep_fname:
                ligprep_fname.write(ligprep_inp_data)
            
            protonation_cpu = str(fetch_ligprep_cpu(schrodinger_dir=schrodinger_dir))
            model_cmd += ['-inp', ligprep_inp, '-NJOBS', protonation_cpu, '-JOBNAME', f'ligprep_{Path(tempfile_dir).name}', '-HOST', f'localhost:{protonation_cpu}', '-WAIT']

            logging.info(f'Ligprep is running with input: {fname}, output: {ligprep_output_fname}, cmd: {model_cmd}')
            subprocess.run(model_cmd, text=True)
            subprocess.run(['cp', Path(ligprep_output_fname).name, ligprep_output_fname], text=True)

    os.chdir(current_dir)
    

def run_easydock_protonation(protonation_data, valid_molecule_file_dict):

    cmd = ['conda', 'run', '-n', protonation_data.environment_value,
            'easydock', '-c', str(calculate_easydock_cpu())]

    if not Path(protonation_data.output_value).is_dir():
        Path(protonation_data.output_value).mkdir(parents=True, exist_ok=True)

    for model, fname in valid_molecule_file_dict.items():

        Path(Path(protonation_data.output_value) / model).mkdir(parents=True, exist_ok=True)
        output_easydock_sdf = str(f'{Path(protonation_data.output_value) / model / Path(fname).stem}_protonated.sdf')
        if Path(output_easydock_sdf).is_file():
            logging.info(f'Protonation with {output_easydock_sdf} has been done. Skipping...')
            continue
        with tempfile.TemporaryDirectory() as tempfile_dir:      
            output_easydock = str(Path(tempfile_dir) / f'{Path(fname).stem}_protonated.db')
            model_cmd = deepcopy(cmd)                               
            model_cmd += ['-i', fname, '-o', output_easydock, '--protonation', protonation_data.method_value]

            logging.info(f'Running easydock protonation with cmd: {model_cmd=}')
            subprocess.run(model_cmd, text=True)

            Path(Path(protonation_data.output_value) / model).mkdir(parents=True, exist_ok=True)
            FETCH_SDF_SCRIPT = Path(__file__).parent / 'utils' / 'easydock_fetch_sdf_from_db.py'
            protonated_cmd = ['conda', 'run', '-n', protonation_data.environment_value,
                            'python', FETCH_SDF_SCRIPT, '-i', output_easydock, '-o', output_easydock_sdf]

            subprocess.run(protonated_cmd, text=True)


def run_protonation(config_data):
    redocking_data = config_redocking.RedockingConfig(config_data[config_constant.ANALYSIS_KEY])
    parameter_data = config_parameter.ParameterConfig(config_data=config_data)
    redocking_data.validate_config()
    protonation_data = redocking_data.protonation_value
    valid_molecule_file_dict = fetch_valid_prepared_molecule_file(dir_path=protonation_data.input_value,
                                                                  parameter_class=parameter_data)

    if 'schrodinger' in str(protonation_data.environment_value).lower():
        run_ligprep_protonation(schrodinger_dir=protonation_data.environment_value,
                                output_dir=protonation_data.output_value,
                                valid_molecule_file_dict=valid_molecule_file_dict)
    else:
        run_easydock_protonation(protonation_data=protonation_data,
                                 valid_molecule_file_dict=valid_molecule_file_dict)


def prepare_easydock_config(docking_data: config_redocking.EasyDockConfig, input_data: config_input.InputStructure, output_dir):
    config_fname = docking_data.config_value
    protein_pdb = input_data.pdb_value
    pocket_coordinate = input_data.pocket_coord_value.split(',')

    with open(config_fname, 'r') as config_f:
        config_data = yaml.load(config_f, Loader=yaml.SafeLoader)
    
    new_config_fname = Path(output_dir) / f'easydock_{Path(input_data.protein_value).stem}.yml'

    grid_txt = f'''
    center_x = {pocket_coordinate[0]}
    center_y = {pocket_coordinate[1]}
    center_z = {pocket_coordinate[2]}
    size_x = 25.0
    size_y = 25.0
    size_z = 25.0
    '''
    grid_fname = str(Path(protein_pdb).parent / f'{input_data.protein_value}_grid.txt')
    with open(grid_fname, 'w') as grid_f:
        grid_f.write(grid_txt)
    

    config_data['protein'] = VinaProtein(pdb_filepath=protein_pdb, 
                                         prepare_receptor_bin_path=docking_data.protein_preparation_executable_value,
                                         preparation_method=docking_data.protein_preparation_value).pdbqt_filepath 
    config_data['protein_setup'] = grid_fname

    if str(docking_data.docking_value).lower() == 'vina':
        if config_data.get('exhaustiveness') is None:
            config_data['exhaustiveness'] = 32
        if config_data.get('n_poses') is None:
            config_data['n_poses'] = 5
        if config_data.get('ncpu') is None:
            config_data['ncpu'] = max(1, min(len(os.sched_getaffinity(0)) - 5, 5))
        if config_data.get('seed') is None:
            config_data['seed'] = 42
    else:
        logging.warning(f'You are using a different docking program other than vina. Ensure that the config at {new_config_fname} has \
                        the correct format. Please check easydock github for reference!')
    
    with open(new_config_fname, 'w') as config_f:
        yaml.dump(config_data, config_f)

    return new_config_fname


def calculate_easydock_cpu():
    cpu_available = len(os.sched_getaffinity(0))
    CPU_FOR_VINA = 5
    CPU_BUFFER = 5
    return max(1, int((cpu_available - CPU_BUFFER)/ CPU_FOR_VINA))


def run_easydock_docking(docking_data: config_redocking.EasyDockConfig, parameter_data, input_data):
    
    cmd = ['conda', 'run', '-n', docking_data.environment_value,
            'easydock',  '--sdf']

    if not Path(docking_data.output_value).is_dir():
        Path(docking_data.output_value).mkdir(parents=True, exist_ok=True)

    config_fname = prepare_easydock_config(docking_data, input_data, output_dir=docking_data.output_value)
    cmd += ['--config', config_fname]

    valid_molecule_file_dict = fetch_valid_prepared_molecule_file(dir_path=docking_data.input_value,
                                                                  parameter_class=parameter_data)
    if docking_data.protonation_value:
        cmd += ['--protonation', docking_data.protonation_value]
    
    if docking_data.docking_value:
        cmd += ['--program', docking_data.docking_value]
    else:
        logging.warning('You did not specify which docking program to use for easydock')

    cmd += ['-c', str(calculate_easydock_cpu())]
    for model, fname in valid_molecule_file_dict.items():

        output_easydock_sdf =  str(Path(docking_data.output_value) / model / f'{Path(fname).stem}_easydock_redocked.sdf')
        Path(output_easydock_sdf).parent.mkdir(parents=True, exist_ok=True)
        if Path(output_easydock_sdf).is_file():
            logging.info(f'Easydock Docking with {output_easydock_sdf} has been done. Skipping...')
            continue

        with tempfile.TemporaryDirectory() as tempfile_dir:      
            temp_output_easydock_db = str(Path(tempfile_dir) / f'{Path(fname).stem}_redocked.db')
            temp_output_easydock_sdf = str(Path(tempfile_dir) / f'{Path(fname).stem}_redocked.sdf')
            
            model_cmd = deepcopy(cmd)                               
            model_cmd += ['-i', fname, '-o', temp_output_easydock_db, ]


            logging.info(f'Easydock docking with cmd: {model_cmd=}')
            subprocess.run(model_cmd, text=True)
            subprocess.run(['cp', temp_output_easydock_sdf, output_easydock_sdf], text=True)


def fetch_glide_sp_cpu(schrodinger_dir):

    cmd = [f'{schrodinger_dir}/run', 'lictool', 'status']

    schrodinger_license_strings = subprocess.run(cmd, capture_output=True, text=True).stdout
    glide_sp_license_available = int(schrodinger_license_strings.split('GLIDE_SP_DOCKING')[-1].split('\n')[0].split('licenses available')[0].split('Total of')[-1])

    MAX_CPU_USED = 50
    CPU_BUFFER = 5
    GLIDE_SP_LICENSE_PER_CPU = 4

    AVAILABLE_CPU = len(os.sched_getaffinity(0)) - CPU_BUFFER

    return max(1, min(int(glide_sp_license_available / GLIDE_SP_LICENSE_PER_CPU) - CPU_BUFFER, MAX_CPU_USED, AVAILABLE_CPU))


def fetch_ligprep_cpu(schrodinger_dir):

    cmd = [f'{schrodinger_dir}/run', 'lictool', 'status']

    schrodinger_license_strings = subprocess.run(cmd, capture_output=True, text=True).stdout
    glide_sp_license_available = int(schrodinger_license_strings.split('LIGPREP_MAIN')[-1].split('\n')[0].split('licenses available')[0].split('Total of')[-1])

    MAX_CPU_USED = 50
    GLIDE_SP_LICENSE_PER_CPU = 1
    CPU_BUFFER = 5

    AVAILABLE_CPU = len(os.sched_getaffinity(0)) - CPU_BUFFER

    return max(1, min(int(glide_sp_license_available / GLIDE_SP_LICENSE_PER_CPU), MAX_CPU_USED, AVAILABLE_CPU))


def run_glide_docking(docking_data: config_redocking.GlideConfig, parameter_data, input_data):

    grid_filepath = GlideProtein(pdb_filepath=input_data.pdb_value,
                                native_ligand=list(Chem.SDMolSupplier(input_data.sdf_value))[0],
                                grid_output_dirpath=str(Path(input_data.pdb_value).parent),
                                schrodinger_dirpath=docking_data.dir_value,
                                protein_preparation=docking_data.protein_preparation_value).grid_filepath

    valid_molecule_file_dict = fetch_valid_prepared_molecule_file(dir_path=docking_data.input_value,
                                                                    parameter_class=parameter_data)

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

        with tempfile.TemporaryDirectory() as tempfile_dir:  
            os.chdir(tempfile_dir)
            if docking_data.protonation_value:  
                run_ligprep_protonation(schrodinger_dir=docking_data.dir_value, 
                                        output_dir=tempfile_dir, 
                                        valid_molecule_file_dict={model: fname})
                protonated_fname = Path(tempfile_dir) / model / f'{Path(fname).stem}_protonated.sdf'

            with tempfile.TemporaryDirectory() as tempfile_dir2: 
                os.chdir(tempfile_dir2)
                subprocess.run(['cp', protonated_fname, tempfile_dir2], text=True) 
                glide_inp = f'''
                    FORCEFIELD   OPLS_2005
                    FORCEPLANAR   True
                    GRIDFILE   {grid_filepath}
                    LIGANDFILE   {Path(protonated_fname).name}
                    POSE_OUTTYPE   ligandlib_sd
                    POSES_PER_LIG   5
                    POSTDOCK_NPOSE   20
                    POSTDOCKSTRAIN   True
                    PRECISION   SP
                    REWARD_INTRA_HBONDS   {True if docking_data.reward_intra_hbonds_value else False}
                    WRITE_VSDB   False
                '''

                tempfile_dir2_path = Path(tempfile_dir2)
                jobname = f'glide_{Path(fname).stem}'
                with open(tempfile_dir2_path / f'{jobname}.inp', 'w') as glide_inp_f:
                    glide_inp_f.write(glide_inp)

                cmd = [f"{docking_data.dir_value}/glide", tempfile_dir2_path / f'{jobname}.inp', '-OVERWRITE', 
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