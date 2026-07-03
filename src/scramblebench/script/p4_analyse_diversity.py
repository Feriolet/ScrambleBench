"""Calculate diversity of valid de novo ligand."""

import argparse
import logging
import sys
import subprocess
import tempfile
import json

from collections import defaultdict
from typing import Any
from pathlib import Path

import yaml

from scramblebench.script.config_preparation import config_constant, config_parameter, config_diversity
from scramblebench.script.utils.process_data import read_input, fetch_model_file_from_model_dir


logger = logging.getLogger(__name__)


def run_model_diversity(diversity_data: config_diversity.DiversityConfig,
                        cmd: list[str], input_fname: str) -> list[dict[str, Any]]:
    """run diversity for a single model

    Args:
        diversity_data (config_diversity.DiversityConfig): DiversityConfig
        cmd (list[str]): command line for analyse_tsp.py
        input_fname (str): sdf filename for valid de novo ligand

    Returns:
        list[dict[str, Any]]: list of diversity result detailing score, time, and method
    """

    output_list = []
    for method in diversity_data.method_value:
        diversity = method[diversity_data.method_diversity_name]
        distance = method[diversity_data.method_distance_name]

        with tempfile.TemporaryDirectory() as tempfile_dir:

            temp_output =  str(Path(tempfile_dir) / 'json_result.json')
            model_cmd = cmd + ['-i', input_fname, '-o', temp_output,
                        '--diversity', diversity]

            if distance:
                model_cmd += ['--distance', distance]

            print(f'{" ".join(model_cmd)=}')
            subprocess.run(model_cmd, text=True, check=True)

            with open(temp_output, 'r', encoding='utf-8') as diversity_f:
                data = json.load(diversity_f)

            output_list.append(data)

    return output_list


def run_diversity(config_data: dict[str, dict]) -> None:
    """calculate the score, time of each valid diversity metric of each model. This function will send a
    a command line to the actual python file responsible to calculate diversity named analyse_tsp.py,
    which will write a JSON output in the end. I did this way, because the analyse_tsp.py or python_tsp.py,
    require a separate dependencies outside of ScrambleBench.

    Args:
        config_data (dict[str, dict]): user's prepared YAML config as dictionary
    """

    diversity_data = config_diversity.DiversityConfig(config_data[config_constant.ANALYSIS_KEY])
    parameter_data = config_parameter.ParameterConfig(config_data=config_data)

    cmd = ['conda', 'run', '-n', diversity_data.environment_value,
            'python', str(Path(__file__).parent / 'utils/diversity_utils/analyse_tsp.py')]

    for model, input_fname in fetch_model_file_from_model_dir(dir_path=diversity_data.input_value,
                                                        model_list=parameter_data.model_list_value).items():


        model_output_dirpath = Path(Path(diversity_data.output_value) / model)
        model_output_dirpath.mkdir(parents=True, exist_ok=True)

        output_fname = model_output_dirpath / 'diversity_output.json'

        if output_fname.is_file():
            logging.info(f'Diversity analysis is done using {input_fname}. Skipping...')
            continue

        output_dict = defaultdict(list)
        output_dict[model] += run_model_diversity(diversity_data=diversity_data,
                                                  cmd=cmd,
                                                  input_fname=input_fname)

        with open(output_fname, 'w', encoding='utf-8') as diversity_model_f:
            json.dump(output_dict, diversity_model_f, indent=4)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Calculate Diversity of generated compound")

    parser.add_argument("-i", "--input", help="config yaml input file or txt file containing yaml filepath",
                        required=True, type=str)
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
        with open(yaml_file, 'r', encoding='utf-8') as yaml_f:
            data_input = yaml.safe_load(yaml_f)
        run_diversity(data_input)
