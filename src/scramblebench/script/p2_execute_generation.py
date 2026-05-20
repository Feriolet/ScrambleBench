from typing import Any, Callable
from pathlib import Path

import yaml
import argparse
import logging
import sys
import subprocess
from scramblebench.script.config_preparation import config_constant, config_generation, config_model
from scramblebench.script.utils.process_data import read_input


logger = logging.getLogger(__name__)
    

def validate_generation_script(config_data: dict) -> None:
    generation_data = config_generation.GenerationConfig(config_data=config_data)
    model_data = config_model.ModelConfig(config_data=config_data)

    logging.debug('Validating the template generation (or custom) script file')

    script_file = generation_data.script_value
    script_pathfile = Path(script_file)
    assert script_pathfile.is_file() and script_pathfile.suffix == '.sh'
    script_text = script_pathfile.read_text()

    for model, modelstruct in model_data.modelstructure_dict.items():
        modelstruct.conda_env_value
        if modelstruct.conda_env_value is None or modelstruct.conda_env_value == 'non_applicable':
            logging.warning(f'{model} is detected not to be inferred.')
            continue

        if f'${model_data.name}_{model}_{modelstruct.conda_env_name}' not in script_text:
            raise ValueError(f"We did not detect {model} for the inference in the {script_file}. \
                             We detect each model by checking the string '${model_data.name}_{model}_{modelstruct.conda_env_name}'")
    
    return True


def run_inference(yaml_file, config_data, log_status=False) -> None:
    generation_data = config_generation.GenerationConfig(config_data=config_data)
    cmd = ['/bin/bash', generation_data.script_value, yaml_file]
    
    if log_status:
        log_msg = subprocess.run(cmd, check=True, text=True, capture_output=True)
        logging.info(log_msg)
    else:
        subprocess.run(cmd, check=True, text=True)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run de novo molecule generation after p1_generate_config.py")

    parser.add_argument("-i", "--input", help="config yaml input file or txt file containing yaml filepath", required=True, type=str)
    parser.add_argument('--log', help="config yaml input file or txt file containing yaml filepath", type=str)
    args = parser.parse_args()

    log_inference_status = False
    if args.log:
        logging.basicConfig(filename=args.log,
                        level=logging.INFO,
                        format='%(asctime)s - %(module)s: - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
        log_inference_status = True
    else:
        logging.basicConfig(stream=sys.stdout,
                        level=logging.INFO,
                        format='%(asctime)s - %(module)s: - %(levelname)s - %(message)s',
                        datefmt='%Y-%m-%d %H:%M:%S')
    
    logging.info('Running p2_execute_generation.py')
    logging.info('Reading the config filename :)')

    yaml_file_list = read_input(args.input)
    for yaml_file in yaml_file_list:
        data_input = yaml.safe_load(open(yaml_file, 'r'))

        if validate_generation_script(config_data=data_input):
            logging.info('Running inference!')
            run_inference(yaml_file, data_input, args.log)