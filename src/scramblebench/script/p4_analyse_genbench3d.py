from typing import Any, Callable
from pathlib import Path

import yaml
import argparse
import logging
import sys
import subprocess
import os

from scramblebench.script.config_preparation import config_constant, config_genbench3d, config_input
from scramblebench.script.utils.error_handler import DirNotFoundError
from scramblebench.script.utils.process_data import read_input, fetch_model_folder_name, find_file_name_through_regex

from copy import deepcopy
import json
from collections import defaultdict
from enum import Enum
from multiprocessing import Pool

logger = logging.getLogger(__name__)


class CheckpointStatus(Enum):
    COMPLETED = "COMPLETED"
    PENDING = "PENDING"
    FAILED = "FAILED"


class CheckpointManager():

    def __init__(self):
        self.checkpoint_fname = 'checkpoint.json'
        self.state = None

    def from_dict(self, state_dict):

        assert isinstance(state_dict, dict)
        self.state = state_dict
    
        return self
    
    def from_json(self, json_fname):

        assert Path(json_fname).is_file() and Path(json_fname).suffix == '.json'
        self.checkpoint_fname = json_fname
        self.state = self.load_state()

        return self
    
    def load_state(self):
        with open(self.checkpoint_fname, 'r') as checkpoint_f:
            return json.load(checkpoint_f)
        
    def save_state(self, checkpoint_fname=None):
        if checkpoint_fname:
            assert Path(checkpoint_fname).suffix == '.json'
            checkpoint_destination = checkpoint_fname
        else:
            checkpoint_destination = self.checkpoint_fname

        with open(checkpoint_destination, 'w') as checkpoint_f:
            json.dump(self.state, checkpoint_f, indent=4)
    
        return self


def fetch_valid_prepared_molecule_file(config_data) -> dict[str, str]:
    model_for_generation_list = fetch_model_folder_name(config_data=config_data)

    generation_folder_dirpath = Path(config_data[config_constant.ANALYSIS_KEY][config_constant.ANALYSIS_GENBENCH3D_KEY]['input'])
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


def check_schrodinger_status(schrodinger_dir):
    if 'schrodinger' not in Path(schrodinger_dir).name:
        raise ValueError('schrodinger dirname should end in schrodinger string (e.g., /path/to/schrodinger-year-quarter)')
    
    schrodinger_version_year = str(schrodinger_dir.split('schrodinger')[-1].split('-')[0])
    schrodinger_version_quarter = str(schrodinger_dir.split('schrodinger')[-1].split('-')[-1])

    cmd = [f'{str(schrodinger_dir)}/jsc', 'local-server-status']

    # schrodinger jsc is implemented from 2025-4 onwards
    if int(schrodinger_version_year) < 2025:
        return None
    elif schrodinger_version_year == '2025' and schrodinger_version_quarter != '4':
        return None

    schrodinger_status = subprocess.run(cmd, capture_output=True, text=True).stdout
    if 'STOPPED' in schrodinger_status:
        subprocess.run([f'{str(schrodinger_dir)}/jsc', 'local-server-start'], capture_output=True, text=True)


def prepare_genbench3d_input(input_data, genbench3d_data):
    GENBENCH3D_INPUT_DIR = 'structural_input'
    genbench3d_input_dirpath = Path(genbench3d_data.output_value) / GENBENCH3D_INPUT_DIR
    genbench3d_input_dirpath.mkdir(exist_ok=True, parents=True)

    input_pdb = str(genbench3d_input_dirpath / Path(input_data.pdb_value).name)
    input_sdf = str(genbench3d_input_dirpath / Path(input_data.sdf_value).name)
    subprocess.run(['cp', input_data.pdb_value, genbench3d_input_dirpath], text=True)
    subprocess.run(['cp', input_data.sdf_value, genbench3d_input_dirpath], text=True)

    return input_pdb, input_sdf


def prepare_genbench3d_cmd(config_data):
    genbench_data = config_genbench3d.GenBench3DConfig(config_data[config_constant.ANALYSIS_KEY])

    do_sbdd_analysis = False


    valid_molecule_file_dict = fetch_valid_prepared_molecule_file(config_data=config_data)
    genbench_output_dirpath = Path(genbench_data.output_value)
    genbench_data.validate_config()
    if genbench_data.schrodinger_dir_value:
        check_schrodinger_status(genbench_data.schrodinger_dir_value)


    complex_minimisation = ['unminimised']
    if config_constant.INPUT_KEY in config_data:
        do_sbdd_analysis = True
        input_data = config_input.InputStructure(config_data[config_constant.INPUT_KEY])
        input_pdb, input_sdf = prepare_genbench3d_input(input_data=input_data, genbench3d_data=genbench_data)
    
        cmd  = ["conda", "run", "-n", genbench_data.conda_env_value, 
                'python', str(Path(genbench_data.genbench_dir_value) / 'sb_benchmark_mols.py'),
                '-c', genbench_data.genbench_config_value,
                '-p', input_pdb,
                '-n', input_sdf,
                '--do_conf_analysis']

        if genbench_data.do_cancel_protonation_by_genbench_value:
            cmd += ['--cancel_protonation']

        if genbench_data.do_complex_forcefield_minimisation_value:
            complex_minimisation.append('minimised')

    else:
        cmd = ["conda", "run", "-n", genbench_data.conda_env_value, 
               'python', str(Path(genbench_data.genbench_dir_value) / 'benchmark_mols.py'),
               '-c', genbench_data.genbench_config_value]


    json_output_dirpath = Path(genbench_output_dirpath) / 'json_output'
    json_output_dirpath.mkdir(parents=True, exist_ok=True)


    checkpoint_filepath = Path(genbench_output_dirpath) / 'genbench3d_checkpoint.json'
    checkpoint_manager = CheckpointManager()
    if not checkpoint_filepath.is_file():
        checkpoint = defaultdict(dict)
        for model in valid_molecule_file_dict:
            for minimisation in complex_minimisation:
                checkpoint[model][minimisation] = CheckpointStatus.PENDING.value

        checkpoint_manager = checkpoint_manager.from_dict(checkpoint)
        checkpoint_manager.checkpoint_fname = checkpoint_filepath
        checkpoint_manager.save_state()
    else:
        checkpoint_manager = checkpoint_manager.from_json(checkpoint_filepath)

    
    cmd_list = []
    for model, fname in valid_molecule_file_dict.items():
        for minimisation in complex_minimisation:
            cmd_data = {}
            cmd_data['minimisation'] = minimisation
            cmd_data['model'] = model
            if checkpoint_manager.state[model][minimisation] == CheckpointStatus.COMPLETED.value:
                logging.info(f'{fname} {model} has been completed with {minimisation} config. Skipping...')
                continue


            model_cmd = deepcopy(cmd)
            fname_stem = Path(fname).stem

            model_cmd += ['-i', fname]
            model_cmd += ['-o', f'{json_output_dirpath}/{fname_stem}_{minimisation}.json']
            model_cmd += ['--log_output', f'{json_output_dirpath}/{fname_stem}_{minimisation}_sb_benchmark.log']

            if do_sbdd_analysis:
                if minimisation == 'minimised':
                    model_cmd += ['-m']
                if genbench_data.schrodinger_dir_value:
                    glide_output =  Path(genbench_output_dirpath) / 'Glide' / model / minimisation
                    glide_output.mkdir(parents=True, exist_ok=True)
                    model_cmd += ['--glide', '--output_glide_dir', str(glide_output)]

                vina_output = Path(genbench_output_dirpath) / 'Vina' / model / minimisation
                vina_output.mkdir(parents=True, exist_ok=True)
                model_cmd += ['--output_vina_dir', str(vina_output)]

            cmd_data['cmd'] = model_cmd
            cmd_list.append(cmd_data)

    return cmd_list, checkpoint_filepath


def run_single_genbench3d(cmd_data):
    try:
        model = cmd_data['model']
        minimisation = cmd_data['minimisation']
        cmd = cmd_data['cmd']
        logging.info(f"Running Genbench. cmd: {' '.join(cmd)}")
        subprocess.run(cmd, capture_output=True, text=True)

        return CheckpointStatus.COMPLETED.value, model, minimisation

    except (subprocess.CalledProcessError, PermissionError, KeyboardInterrupt) as e:
        return CheckpointStatus.FAILED.value, model, minimisation

def run_genbench3d(config_data):
    cmd_list, checkpoint_fname = prepare_genbench3d_cmd(config_data=config_data)
    checkpoint_manager = CheckpointManager().from_json(checkpoint_fname)

    MAX_CPU_USED = 5
    CPU_BUFFER = 5
    num_cpu_available = len(os.sched_getaffinity(0)) - CPU_BUFFER

    used_cpu = max(1, min(MAX_CPU_USED, num_cpu_available))
    with Pool(used_cpu) as pool:
        for genbench_status, model, minimisation in pool.imap_unordered(run_single_genbench3d, cmd_list):
            checkpoint_manager.state[model][minimisation] = genbench_status
            checkpoint_manager.save_state()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Perform genbench3d analysis after molecule preparation")

    parser.add_argument("-i", "--input", help="config yaml input file or txt file containing yaml filepath", required=True, type=str)
    args = parser.parse_args()

    logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='%(asctime)s - %(module)s: - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
    
    logging.info('Running p4_analyse_genbench3d.py')
    logging.info('Reading the config filename :\)')


    yaml_file_list = read_input(args.input)

    for yaml_file in yaml_file_list:
        checkpoint_json = Path(yaml_file).parent
        data_input = yaml.safe_load(open(yaml_file, 'r'))
        run_genbench3d(config_data=data_input)