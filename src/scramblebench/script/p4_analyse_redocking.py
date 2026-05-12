from typing import Any, Callable
from pathlib import Path

import yaml
import argparse
import logging
import sys
import subprocess
import os

from scramblebench.script.config_preparation import config_constant, config_genbench3d, config_input, config_redocking
from scramblebench.script.utils.error_handler import DirNotFoundError
from scramblebench.script.utils.process_data import read_input, fetch_model_folder_name, find_file_name_through_regex

from copy import deepcopy
import json
from collections import defaultdict
from enum import Enum
from multiprocessing import Pool
import tempfile

def fetch_valid_prepared_molecule_file(config_data) -> dict[str, str]:
    model_for_generation_list = fetch_model_folder_name(config_data=config_data)

    generation_folder_dirpath = Path(config_data[config_constant.ANALYSIS_KEY][config_constant.ANALYSIS_REDOCKING_KEY][config_constant.ANALYSIS_REDOCKING_PROTONATION_KEY]['input'])
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

        valid_molecule_file_dict[model] = str(matched_fname_list[0])
    
    return valid_molecule_file_dict


def run_ligprep_protonation(protonation_data, valid_molecule_file_dict):

    cmd = [f'{protonation_data.environment_value}/ligprep']
    if not Path(protonation_data.output_value).is_dir():
        Path(protonation_data.output_value).mkdir(parents=True, exist_ok=True)

    for fname in valid_molecule_file_dict.values():

        with tempfile.TemporaryDirectory() as tempfile_dir:      

            model_cmd = deepcopy(cmd)

            ligprep_inp = str(Path(tempfile_dir) / f'{Path(fname).stem}_ligprep.inp')
            ligprep_output_fname = str(f'{Path(protonation_data.output_value) / Path(fname).stem}_protonated.sdf')
            ligprep_inp_data = f'''
            INPUT_FILE_NAME   {fname}
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
            OUT_SD  {ligprep_output_fname}

            '''

            with open(ligprep_inp, 'w') as ligprep_fname:
                ligprep_fname.write(ligprep_inp_data)
            
            model_cmd += ['-inp', ligprep_inp, '-NJOBS', '1', '-JOBNAME', f'ligprep_{Path(tempfile_dir).name}', '-HOST', 'localhost:1', '-WAIT']
            subprocess.run(model_cmd, text=True)

def run_easydock_protonation(protonation_data, valid_molecule_file_dict):
    print(protonation_data.environment_value)
    cmd = ['conda', 'run', '-n', protonation_data.environment_value,
            'easydock']

    if not Path(protonation_data.output_value).is_dir():
        Path(protonation_data.output_value).mkdir(parents=True, exist_ok=True)

    for fname in valid_molecule_file_dict.values():

        with tempfile.TemporaryDirectory() as tempfile_dir:      
            output_easydock = str(Path(tempfile_dir) / f'{Path(fname).stem}_protonated.db')
            model_cmd = deepcopy(cmd)                               
            model_cmd += ['-i', fname, '-o', output_easydock, '--protonation', protonation_data.method_value]

            print(f'{model_cmd=}')
            subprocess.run(model_cmd, text=True)

            FETCH_SDF_SCRIPT = Path(__file__).parent / 'utils' / 'easydock_fetch_sdf_from_db.py'
            protonated_cmd = ['conda', 'run', '-n', protonation_data.environment_value,
                            'python', FETCH_SDF_SCRIPT, '-i', output_easydock, '-o', str(f'{Path(protonation_data.output_value) / Path(fname).stem}_protonated.sdf')]

            subprocess.run(protonated_cmd, text=True)

def run_protonation(config_data):
    redocking_data = config_redocking.RedockingConfig(config_data[config_constant.ANALYSIS_KEY])
    redocking_data.validate_config()
    protonation_data = redocking_data.protonation_value
    valid_molecule_file_dict = fetch_valid_prepared_molecule_file(config_data=config_data)

    if 'schrodinger' in str(protonation_data.environment_value).lower():
        run_ligprep_protonation(protonation_data=protonation_data,
                                valid_molecule_file_dict=valid_molecule_file_dict)
    else:
        run_easydock_protonation(protonation_data=protonation_data,
                                 valid_molecule_file_dict=valid_molecule_file_dict)

            
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Prepare generated molecules for downstream analysis")

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
        run_protonation(config_data=data_input)