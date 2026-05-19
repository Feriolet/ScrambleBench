from typing import Any, Callable
from pathlib import Path

import yaml
import argparse
import logging
import sys
import subprocess
import json
from scramblebench.script.config_preparation import config_constant, config_parameter, config_diversity
from scramblebench.script.utils.error_handler import DirNotFoundError
from scramblebench.script.utils.process_data import read_input, fetch_model_file_from_model_dir

import os

import tempfile
from copy import deepcopy
from collections import defaultdict


logger = logging.getLogger(__name__)


def run_diversity(config_data):
    diversity_data = config_diversity.DiversityConfig(config_data[config_constant.ANALYSIS_KEY])
    parameter_data = config_parameter.ParameterConfig(config_data=config_data)

    valid_molecule_file_dict = fetch_model_file_from_model_dir(dir_path=diversity_data.input_value,
                                                               model_list=parameter_data.model_list_value)
    
    diversity_output_dirpath = diversity_data.output_value

    cmd = ['conda', 'run', '-n', diversity_data.environment_value,
            'python', str(Path(__file__).parent / 'utils/diversity_utils/analyse_tsp.py')]
    
    for model, fname in valid_molecule_file_dict.items():
        diversity_output_dict = defaultdict(list)
        diversity_model_dirpath = Path(Path(diversity_output_dirpath) / model)
        diversity_model_dirpath.mkdir(parents=True, exist_ok=True)

        diversity_model_output_fname = diversity_model_dirpath / 'diversity_output.json'

        if diversity_model_output_fname.is_file():
            logging.info(f'Diversity analysis is done using {fname}. Skipping...')
            continue

        for method in diversity_data.method_value:
            diversity = method[diversity_data.method_diversity_name]
            distance = method[diversity_data.method_distance_name]
        
            with tempfile.TemporaryDirectory() as tempfile_dir:

                model_cmd = deepcopy(cmd)
                temp_output =  str(Path(tempfile_dir) / 'json_result.json')
                model_cmd += ['-i', fname, '-o', temp_output,
                            '--diversity', diversity]
                
                if distance:
                    model_cmd += ['--distance', distance]

                print(f'{" ".join(model_cmd)=}')
                subprocess.run(model_cmd, text=True)

                with open(temp_output, 'r') as diversity_f:
                    data = json.load(diversity_f)
                
                diversity_output_dict[model].append(data)
                
        with open(diversity_model_output_fname, 'w') as diversity_model_f:
            json.dump(diversity_output_dict, diversity_model_f, indent=4)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate Diversity of generated compound")

    parser.add_argument("-i", "--input", help="config yaml input file or txt file containing yaml filepath", required=True, type=str)
    args = parser.parse_args()

    logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='%(asctime)s - %(module)s: - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
    
    logging.info('Running p4_analyse_diversity.py')
    logging.info('Reading the config filename :)')


    yaml_file_list = read_input(args.input)

    for yaml_file in yaml_file_list:
        checkpoint_json = Path(yaml_file).parent
        data_input = yaml.safe_load(open(yaml_file, 'r'))
        run_diversity(data_input)
        