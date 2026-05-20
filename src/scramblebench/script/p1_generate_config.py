
from typing import Any, Callable
from pathlib import Path
import yaml
import argparse
import logging
import sys
from copy import deepcopy

from scramblebench.script.config_preparation.config_model import ModelConfig
from scramblebench.script.config_preparation.config_generation import GenerationConfig
from scramblebench.script.config_preparation.config_post_generation import PostGenerationConfig
from scramblebench.script.config_preparation.config_analysis import AnalysisConfig
from scramblebench.script.config_preparation.config_parameter import ParameterConfig
from scramblebench.script.config_preparation import config_constant
import itertools 


logger = logging.getLogger(__name__)


def load_config(config_fname: str) -> dict[str, Any]:
    with open(config_fname, 'r') as config_fn:
        return yaml.safe_load(config_fn)


def check_correct_input_output_folder(prestep: Callable, poststep: Callable) -> None:

    if not str(prestep.output_value) == str(poststep.input_value):
        raise TypeError(f'The output folder of {prestep.name} and input folder of {poststep.name} is not the same')


def validate_config(config_data: dict[str, Any]) -> None:

    compulsory_config_list = [InputConfig, ModelConfig, GenerationConfig, PostGenerationConfig]

    for config_class in compulsory_config_list:
        config_class(config_data).validate_config()

    check_correct_input_output_folder(prestep=GenerationConfig(config_data),
                                      poststep=PostGenerationConfig(config_data))
    
    if config_constant.ANALYSIS_KEY in config_data:
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


def prepare_config(config_data: dict[str, Any]) -> dict[str, Any]:

    config_data |= InputConfig(config_data).write()
    config_data |= GenerationConfig(config_data).write()

    return config_data


def fetch_batch_parameters(config_data):

    batch_parameters : list[dict[str]] = deepcopy(config_constant.CONFIG_BATCH_PARAMETERS)
    for parameter_dict in batch_parameters:
        parameter_values = deep_get(config_data, parameter_dict['key'])
        if isinstance(parameter_values, dict):
            parameter_dict['value'] = [{key: value} for key, value in parameter_values.items()]
        elif isinstance(parameter_values, list):
            parameter_dict['value'] = parameter_values

    batch_parameters = [parameter_dict for parameter_dict in batch_parameters if parameter_dict.get('value')]   

    return batch_parameters


def write_new_config(config_data, prefix_dir, batch_parameter):

    config_output = {}
    logging.debug('Writing Config for Input key')
    config_output = config_output | InputConfig(config_data).write_inputstruct(cutoff=10)
    logging.debug('Writing Config for Model key')
    config_output = config_output | ModelConfig(config_data).write()
    logging.debug('Writing Config for Generation key')
    config_output = config_output | GenerationConfig(config_data).write(prefix_dir=prefix_dir)
    logging.debug('Writing Config for Parameter key')
    config_output = config_output | ParameterConfig().create(config_data=config_output, batch_parameter=batch_parameter).write()
    logging.debug('Writing Config for PostGeneration key')
    config_output = config_output | PostGenerationConfig(config_data).write(prefix_dir=prefix_dir)
    
    if config_constant.ANALYSIS_KEY in config_data:
        logging.debug('Writing Config for Analysis key')
        config_output = config_output | AnalysisConfig(config_data).write(prefix_dir=prefix_dir)

    return config_output


def write_config(config_data: dict[str, Any], output_fname: str) -> None:

    config_data = prepare_config(config_data=config_data)

    batch_parameters: list[dict[str, list]] = fetch_batch_parameters(config_data=config_data) 

    parameter_value_lists = [parameter_dict['value'] for parameter_dict in batch_parameters]
    parameter_key_lists = [parameter_dict['key'] for parameter_dict in batch_parameters]
    
    yaml_list = []
    generation_dirpath = Path(GenerationConfig(config_data).input_value).resolve()

    for assigned_parameter_values in list(itertools.product(*parameter_value_lists)): # [ [val, val, val], [val, val, val] ]
        
        assigned_config_data = deepcopy(config_data)
        param_config_batch_parameter_dict = {}
        analysis_dirpath = generation_dirpath

        for key, value, dtype in zip(parameter_key_lists, assigned_parameter_values, [parameter_dict['type'] for parameter_dict in batch_parameters]):
            assigned_config_data = deep_assign(assigned_config_data, key, value=forcetype(value, dtype))

            if key[0] != config_constant.INPUT_KEY:
                param_config_batch_parameter_dict[key[-1]] = value

            if isinstance(value, dict):
                value = list(value.keys())[0]
            analysis_dirpath = Path(analysis_dirpath) / f'{key[-1]}_{value}'

        config_output = write_new_config(config_data=assigned_config_data, 
                                         prefix_dir=analysis_dirpath, 
                                         batch_parameter=param_config_batch_parameter_dict)

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
    if True:
        write_config(data_input, args.output)