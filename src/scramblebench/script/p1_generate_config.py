
from typing import Any, Callable
from pathlib import Path
from dataclasses import dataclass
import yaml
import argparse
import logging
import sys
from copy import deepcopy

from scramblebench.script.config_preparation.config_model import ModelConfig
from scramblebench.script.config_preparation.config_generation import GenerationConfig
from scramblebench.script.config_preparation.config_post_generation import PostGenerationConfig
from scramblebench.script.config_preparation.config_genbench3d import GenBench3DConfig
from scramblebench.script.config_preparation.config_analysis import AnalysisConfig
from scramblebench.script.config_preparation.config_parameter import ParameterConfig
from scramblebench.script.config_preparation import config_constant
import itertools 
import subprocess

logger = logging.getLogger(__name__)


def load_config(config_fname: str) -> dict[str, Any]:
    with open(config_fname, 'r') as config_fn:
        return yaml.safe_load(config_fn)


def check_correct_input_output_folder(prestep: Callable, poststep: Callable) -> None:

    if not str(prestep.output_value) == str(poststep.input_value):
        raise TypeError(f'The output folder of {prestep.name} and input folder of {poststep.name} is not the same')


def validate_config(config_data: dict[str, Any]) -> None:

    compulsory_config_key_list = [
        config_constant.INPUT_KEY,
        config_constant.MODEL_KEY,
        config_constant.GENERATION_KEY,
        config_constant.POST_GENERATION_KEY
    ]

    for config_key in compulsory_config_key_list:
        if config_key not in config_data.keys():
            if config_key == config_constant.INPUT_KEY and config_constant.INPUT_DIR_KEY not in config_data.keys():
                raise ValueError(f'{config_key} not in your config file. Please add it with case sensitive letters.')
        
    InputConfig(config_data).validate_config()
    ModelConfig(config_data).validate_config()
    GenerationConfig(config_data).validate_config()
    PostGenerationConfig(config_data).validate_config()

    check_correct_input_output_folder(prestep=GenerationConfig(config_data),
                                      poststep=PostGenerationConfig(config_data))
    
    if config_constant.ANALYSIS_KEY in config_data.keys():
        AnalysisConfig(config_data).validate_config()
    
    return True

#https://github.com/yaml/pyyaml/issues/127
class MyDumper(yaml.SafeDumper):
    # HACK: insert blank lines between top-level objects
    # inspired by https://stackoverflow.com/a/44284819/3786245
    def write_line_break(self, data=None):
        super().write_line_break(data)

        if len(self.indents) == 1:
            super().write_line_break()


def deep_get(dictionary: dict, nested_key: list):
    copied_dict = deepcopy(dictionary)
    for key in nested_key:
        copied_dict = copied_dict.get(key)

        if copied_dict is None:
            logging.exception(f'The dictionary {dictionary} has no value for key {key}')
            raise KeyError(f'Key {key} not found for dictionary {dictionary}')
        if not isinstance(copied_dict, dict) and key != nested_key[-1]:
            logging.warning(f'Your dictionary {dictionary} overshoot the nested key. The value for {key} is not a dictionary. Ignoring subsequent keys')
            return copied_dict
    
    return copied_dict

#https://stackoverflow.com/questions/13687924/setting-a-value-in-a-nested-python-dictionary-given-a-list-of-indices-and-value
def deep_assign(dictionary: dict, nested_key: list, value: Any, create_missing=False):
    d = dictionary
    for key in nested_key[:-1]:
        if key in d:
            d = d[key]
        elif create_missing:
            d = d.setdefault(key, {})
        else:
            return dictionary
    if nested_key[-1] in d or create_missing:
        d[nested_key[-1]] = value
    return dictionary


def forcetype(value: Any, dtype='int'):
    if dtype == 'int':
        return int(value)
    elif dtype == 'float':
        return float(value)
    elif dtype == 'str':
        return str(value)
    elif dtype == 'dict':
        if isinstance(value, dict):
            return value
        else:
            logging.exception(f'{value} type is not dict, but instead {type(value)}')
    else:
        logging.exception(f'{value} has unsupported type of {dtype} requested')
        raise ValueError(f'unsupported dtype {dtype}. Please enter int, float, or str only')


def reassign_input_config(config_data):
    'to prevent subsequent script to identify the protein name of the input key'

    input_data = config_data[config_constant.INPUT_KEY]
    assert len(list(input_data.keys())) == 1

    config_data[config_constant.INPUT_KEY] = {}
    for values in input_data.values():
        config_data[config_constant.INPUT_KEY].update(values)

    config_data[config_constant.INPUT_KEY]['name'] = list(input_data.keys())[0]

    return config_data


def write_config(config_data: dict[str, Any], output_fname: str) -> None:

    repeat_parameter = [{'key':[config_constant.INPUT_KEY],
                         'type': 'dict'},
                        {'key':[config_constant.GENERATION_KEY, config_constant.GENERATION_PARAMETER_KEY, 'num_sample'],
                        'type': 'int'},
                        {'key':[config_constant.GENERATION_KEY, config_constant.GENERATION_PARAMETER_KEY, 'box_size'],
                        'type': 'float'}]
    
    config_data |= InputConfig(config_data).write()
    config_data |= GenerationConfig(config_data).write()
    
    for repeat_dict in repeat_parameter:
        nested_value = deep_get(config_data, repeat_dict['key'])
        if isinstance(nested_value, dict):
            repeat_dict['value'] = [{key: value} for key, value in nested_value.items()]
        else:
            repeat_dict['value'] = nested_value

    repeat_parameter = [repeat_dict for repeat_dict in repeat_parameter if isinstance(repeat_dict['value'], list)]   

    repeat_parameter_dict = {}
    for parameter in repeat_parameter:
        parameter = parameter['key']
        if parameter[-1] == config_constant.INPUT_KEY:
            continue
        # add repeat parameter key in config_generation.py
        repeat_parameter_dict[parameter[-1]] = ','.join(parameter) 
    
    generation_data = GenerationConfig(config_data)
    generation_dirpath = Path(generation_data.input_value).resolve()
    
    config_data |= generation_data.update('repeat_parameter', repeat_parameter_dict).write()
    config_data |= PostGenerationConfig(config_data).write()
    
    nested_value_lists = [repeat_dict['value'] for repeat_dict in repeat_parameter]
    nested_key_lists = [repeat_dict['key'] for repeat_dict in repeat_parameter]
    nested_type_lists = [repeat_dict['type'] for repeat_dict in repeat_parameter]
    combinatorial_value_end_list = list(itertools.product(*nested_value_lists))

    yaml_list = []
    for each_combination_list in combinatorial_value_end_list:
        
        assigned_config_data = deepcopy(config_data)
        for key, value, dtype in zip(nested_key_lists, each_combination_list, nested_type_lists):
            assigned_config_data = deep_assign(assigned_config_data, key, value=forcetype(value, dtype))

        analysis_dirpath = generation_dirpath
        for key, val in zip(nested_key_lists, each_combination_list):
            if isinstance(val, dict):
                val = list(val.keys())[0]
            analysis_dirpath = Path(analysis_dirpath) / f'{key[-1]}_{val}'

        config_output = {}
        logging.debug('Writing Config for Input key')
        config_output = config_output | InputConfig(assigned_config_data).write(cutoff=10)
        logging.debug('Writing Config for Model key')
        config_output = config_output | ModelConfig(assigned_config_data).write()
        logging.debug('Writing Config for Generation key')
        config_output = config_output | GenerationConfig(assigned_config_data).write(prefix_dir=analysis_dirpath)
        logging.debug('Writing Config for Parameter key')
        config_output = config_output | ParameterConfig().create(config_data=config_output).write()
        logging.debug('Writing Config for PostGeneration key')
        config_output = config_output | PostGenerationConfig(assigned_config_data).write(prefix_dir=analysis_dirpath)
        logging.debug('Writing Config for Analysis key')
        config_output = config_output | AnalysisConfig(assigned_config_data).write(prefix_dir=analysis_dirpath)


        config_output = reassign_input_config(config_output)
        Path(analysis_dirpath).mkdir(parents=True, exist_ok=True)
        yaml_fname = Path(analysis_dirpath) / Path(output_fname).name
        with open(yaml_fname, 'w') as yaml_f:
            yaml.dump(config_output, yaml_f, Dumper=MyDumper, sort_keys=False)
            yaml_list.append(str(yaml_fname))
            logging.info(f"Config saved in {yaml_fname}")

    logging.info(f"Preparing Script for Batch Generation")

    with open(Path(generation_dirpath) / 'yaml_list.txt', 'w') as generation_fname:
        for yaml_f in yaml_list:
            generation_fname.write(f'{yaml_f} \n')

    logging.info(f"Yaml file list saved in {str(Path(generation_dirpath) / 'yaml_list.txt')} for job manager and p2_execute_generation.py")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Prepare config file for ScrambleBench")

    parser.add_argument("-i", "--input", help="config yaml input file", type=str, default='/mnt/sod/Veincent/manuscript/ScrambleBench/tests/config_test_run.yml')
    parser.add_argument("-o", "--output", help="config yaml output file prefix (i.e., not path directory)", type=str)
    parser.add_argument("--dirpath_input", action='store_true', help='write input key as a single directory', default=False)
    args = parser.parse_args()

    logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='%(asctime)s - %(module)s: - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
    
    logging.info('Running p1_generate_config.py')
    logging.info('Reading the config filename :)')
    data_input = yaml.safe_load(open(args.input, 'r'))

    if args.dirpath_input:
        logging.info('Dirpath_input mode: expecting to read directories given a directory')
        from scramblebench.script.config_preparation.config_input import InputDirConfig as InputConfig
    else:
        from scramblebench.script.config_preparation.config_input import InputConfig
    if not args.output:
        args.output = Path(args.input).parent / f'{Path(args.input).stem}_clean_config.yml'
    if Path(args.output).suffix not in ['.yaml', '.yml']:
        logging.exception('Error failed in config output filename format')
        raise ValueError(f'{args.output} ends with {Path(args.output).suffix}. Only .yaml and .yml extension is allowed')
    
    args.output = str(Path(args.output))
    if validate_config(data_input):
        write_config(data_input, args.output)