from typing import Any, Callable
from pathlib import Path

import yaml
import argparse
import logging
import sys
import subprocess

from scramblebench.script.config_preparation import config_constant

logger = logging.getLogger(__name__)

def read_input(input_fname: str) -> list[str]:
    input_filepath = Path(input_fname)
    
    if not input_filepath.is_file():
        raise FileNotFoundError(f'The file {input_fname} is not found. Please check your directory')
    
    if input_filepath.suffix.lower() in ['.yaml', '.yml']:
        return [input_filepath.resolve()]
    elif input_filepath.suffix.lower() in ['.txt', '.text']:
        logging.debug(f'Input is in suffix {input_filepath.suffix}. Reading content to fetch yaml files.')
        yaml_file_list = []
        input_content = input_filepath.read_text().splitlines()

        for yaml_file in input_content:
            yaml_filepath = Path(yaml_file)
            if not yaml_filepath.is_file():
                raise FileNotFoundError(f'The file {yaml_file} is not found. Please check your directory')
    
            if yaml_filepath.suffix.lower() not in ['.yaml', '.yml']:
                raise ValueError(f'Incorrect file {yaml_file}. We only support yaml file (i.e., .yaml or .yml). Please use the p1_generate_config.py to prepare your config file.')

            yaml_file_list.append(yaml_file)

        return yaml_file_list
    

def fetch_model(config_data):
    return list(config_data[config_constant.MODEL_KEY].keys())


def validate_generation_script(config_data: dict) -> None:

    logging.debug('Validating the template generation (or custom) script file')
    script_file = config_data[config_constant.GENERATION_KEY][config_constant.GENERATION_TEMPLATE_SCRIPT_FILE_KEY]
    script_pathfile = Path(script_file)
    assert script_pathfile.is_file() and script_pathfile.suffix == '.sh'

    script_text = script_pathfile.read_text()

    model_for_generation_list = fetch_model(config_data=config_data)

    for model in model_for_generation_list:
        # here, we check whether the model is executed through each model conda environment
        if config_data[model][config_constant.MODEL_CONDA_ENV_KEY] is None or config_data[model][config_constant.MODEL_CONDA_ENV_KEY] == 'non_applicable':
            logging.warning(f'{model} is detected not to be inferred.')
            continue

        if f'${config_constant.MODEL_KEY}_{model}_{config_constant.MODEL_CONDA_ENV_KEY}' not in script_text:
            raise ValueError(f"We did not detect {model} for the inference in the {script_file}. We detect each model by checking the string '${config_constant.MODEL_KEY}_{model}_{config_constant.MODEL_CONDA_ENV_KEY}'")
    
    return True


def run_inference(yaml_file, config_data) -> None:
    cmd = [config_data[config_constant.GENERATION_KEY][config_constant.MODEL_CONDA_ENV_KEY],
           yaml_file]
    
    subprocess.run(cmd, shell=True, capture_output=True)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Run de novo molecule generation after p1_generate_config.py")

    parser.add_argument("-i", "--input", help="config yaml input file or txt file containing yaml filepath", required=True, type=str)
    args = parser.parse_args()

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
            run_inference(yaml_file, data_input)