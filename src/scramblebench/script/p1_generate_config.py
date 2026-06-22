"""Preparing ScrambleBench Config before running generation and analysis"""
from typing import Any, Callable, Optional
from pathlib import Path
from copy import deepcopy

import argparse
import logging
import sys
import itertools

import yaml

from scramblebench.script.config_preparation.config_model import ModelConfig
from scramblebench.script.config_preparation.config_generation import GenerationConfig
from scramblebench.script.config_preparation.config_post_generation import PostGenerationConfig
from scramblebench.script.config_preparation.config_analysis import AnalysisConfig
from scramblebench.script.config_preparation.config_parameter import ParameterConfig
from scramblebench.script.config_preparation.config_report import ReportConfig
from scramblebench.script.config_preparation.config_input import InputConfig as InputUserConfig
from scramblebench.script.config_preparation import config_constant


logger = logging.getLogger(__name__)


def load_config(config_fname: str) -> dict[str, Any]:
    """generic function to load YAML config

    Args:
        config_fname (str): config filename by user

    Returns:
        dict[str, Any]: yaml dictionary
    """
    with open(config_fname, 'r', encoding='utf-8') as config_fn:
        return yaml.safe_load(config_fn)


def check_correct_input_output_folder(prestep: Callable, poststep: Callable) -> None:
    """Check if the name of output folder in prestep matches the name of the input folder in poststep

    Args:
        prestep (Callable): ScrambleBench Config class responsible for previous step
        poststep (Callable): ScrambleBench Config class responsible for subsequent step

    Raises:
        TypeError: raised when the absolute directory does not matched
    """

    if str(prestep.output_value) != str(poststep.input_value):
        raise TypeError(f'The output folder of {prestep.name} and input folder of {poststep.name} is not the same')


def validate_config(config_data: dict[str, Any]) -> Optional[bool]:

    """validate user YAML config before running the ScrambleBench

    Returns:
        _type_: returns True if a user config is validated
    """
    compulsory_config_list = [InputUserConfig, ModelConfig, GenerationConfig, PostGenerationConfig]

    for config_class in compulsory_config_list:
        config_class(config_data).validate_config()

    check_correct_input_output_folder(prestep=GenerationConfig(config_data),
                                      poststep=PostGenerationConfig(config_data))

    if config_constant.ANALYSIS_KEY in config_data:
        AnalysisConfig(config_data).validate_config()

    if config_constant.REPORT_KEY in config_data:
        ReportConfig(config_data).validate_config()

    return True


#https://github.com/yaml/pyyaml/issues/127
class MyDumper(yaml.SafeDumper):
    """YAML Dumper to provide spacing between key/field
    """
    # HACK: insert blank lines between top-level objects
    # inspired by https://stackoverflow.com/a/44284819/3786245
    def write_line_break(self, data=None):
        super().write_line_break(data)

        if len(self.indents) == 1:
            super().write_line_break()


def deep_get(dictionary: dict, nested_key: list) -> Any:

    """fetch the value of nested dictionaries given a key/subkeys in a list

    Raises:
        KeyError: if the key/subkey does not exist in the dictionary

    Returns:
        Any: value of subkey/key
    """
    copied_dict = deepcopy(dictionary)
    for key in nested_key:
        copied_dict = copied_dict.get(key)

        if copied_dict is None:
            logging.exception(f'The dictionary {dictionary} has no value for key {key}')
            raise KeyError(f'Key {key} not found for dictionary {dictionary}')
        if not isinstance(copied_dict, dict) and key != nested_key[-1]:
            logging.warning(f'Your dictionary {dictionary} overshoot the nested key. \
                            The value for {key} is not a dictionary. Ignoring subsequent keys')
            return copied_dict

    return copied_dict

#https://stackoverflow.com/questions/13687924/setting-a-value-in-a-nested-python-dictionary-given-a-list-of-indices-and-value
def deep_assign(dictionary: dict, nested_key: list, value: Any, create_missing=False) -> dict:
    """assign a value to nested dictionaries given key/subkeys as a list

    Args:
        dictionary (dict): the dictionary that needs to be updated
        nested_key (list): the key or subkey that exists in the dictionary
        value (Any): the value that will update the dictionary
        create_missing (bool, optional): whether to create the key/subkey into the dictionary. Defaults to False.

    Returns:
        dict: the updated dictionary
    """
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


def forcetype(value: Any, dtype='int') -> Any:
    """Force the type of the value given the dtype

    Args:
        value (Any): the value that will be updated to new dtype
        dtype (str, optional): dtype. Defaults to 'int'.

    Raises:
        ValueError: if dtype does not match [int, float, str, dict]

    Returns:
        Any: updated value with the assigned dtype
    """
    if dtype == 'int':
        return int(value)
    if dtype == 'float':
        return float(value)
    if dtype == 'str':
        return str(value)
    if dtype == 'dict':
        if isinstance(value, dict):
            return value

        logging.exception(f'{value} type is not dict, but instead {type(value)}')
        raise ValueError(f'value {value} is not a dictionary, but {type(value)}')

    logging.exception(f'{value} has unsupported type of {dtype} requested')
    raise ValueError(f'unsupported dtype {dtype}. Please enter int, float, or str only')


def prepare_config(config_data: dict[str, Any]) -> dict[str, Any]:
    """prepare ScrambleBench config to match the "write_config" functions.

    Args:
        config_data (dict[str, Any]): the user config

    Returns:
        dict[str, Any]: the prepared config with the correct structure of input and generation config
    """
    config_data |= InputUserConfig(config_data).write()
    config_data |= GenerationConfig(config_data).write()

    return config_data


def fetch_batch_parameters(config_data: dict[str, Any]) -> list[dict[str, Any]]:
    """fetch the parameters that have multiple values (named batch parameters).

    Args:
        config_data (dict[str, Any]): user YAML dictionary

    Returns:
        list[dict[str, Any]]: list of batch parameters with their multiple values
    """

    batch_parameters : list[dict[str, Any]] = deepcopy(config_constant.CONFIG_BATCH_PARAMETERS)
    for parameter_dict in batch_parameters:
        parameter_values = deep_get(config_data, parameter_dict['key'])
        if isinstance(parameter_values, dict):
            parameter_dict['value'] = [{key: value} for key, value in parameter_values.items()]
        elif isinstance(parameter_values, list):
            parameter_dict['value'] = parameter_values

    batch_parameters = [parameter_dict for parameter_dict in batch_parameters if parameter_dict.get('value')]

    return batch_parameters


def write_new_config(config_data: dict[str, Any],
                     prefix_dir: str,
                     batch_parameter: list[dict[str, Any]]) -> dict[str, Any]:
    """prepare ScrambleBench config based on the user YAML config data

    Args:
        config_data (dict[str, Any]): prepared config data with assigned batch parameters (should only have 1 value)
        prefix_dir (str): root directory for the ScrambleBench output 
        batch_parameter (list[dict[str, Any]]): list of batch parameters with single value. 
                                                This will be saved for ParameterConfig.

    Returns:
        dict[str, Any]: prepared ScrambleBench config as dictionary
    """
    config_output = {}
    logging.debug('Writing Config for Input key')
    config_output = config_output | InputUserConfig(config_data).write_inputstruct(cutoff=10)
    logging.debug('Writing Config for Model key')
    config_output = config_output | ModelConfig(config_data).write()
    logging.debug('Writing Config for Generation key')
    config_output = config_output | GenerationConfig(config_data).write(prefix_dir=prefix_dir)
    logging.debug('Writing Config for Parameter key')
    config_output = config_output | ParameterConfig().create(config_data=config_output,
                                                             batch_parameter=batch_parameter).write()
    logging.debug('Writing Config for PostGeneration key')
    config_output = config_output | PostGenerationConfig(config_data).write(prefix_dir=prefix_dir)

    if config_constant.ANALYSIS_KEY in config_data:
        logging.debug('Writing Config for Analysis key')
        config_output = config_output | AnalysisConfig(config_data).write(prefix_dir=prefix_dir)

    if config_constant.REPORT_KEY in config_data:
        logging.debug('Writing Config for Report key')
        config_output = config_output | ReportConfig(config_data).write()


    return config_output


def process_single_batch_combination(root_dir: str,
                                     config_data: dict[str, Any],
                                     batch_parameters: list[dict[str, Any]],
                                     assigned_parameter_values: list[list[Any]]) \
                                    -> tuple[dict[str, Any], str, dict[str,Any]]:
    """assigned config data according to the each combination of batch parameters

    Args:
        root_dir (str): output root dir
        config_data (dict[str, Any]): user YAML config data that has been prepared into the correct structure
        batch_parameters (list[dict[str, Any]]): list of dictionary of batch parameter (1 value)
        assigned_parameter_values (list[list[Any]]): assigned values of the batch parameter

    Returns:
        tuple[dict[str, Any], str, dict[str,Any]]: (assigned/final config data, 
                                                    analysis root directory, 
                                                    parameter config as dictionary)
    """
    param_config_batch_parameter_dict = {}

    parameter_key_lists = [parameter_dict['key'] for parameter_dict in batch_parameters]
    parameter_dtype_lists = [parameter_dict['type'] for parameter_dict in batch_parameters]

    for key, value, dtype in zip(parameter_key_lists, assigned_parameter_values, parameter_dtype_lists):
        config_data = deep_assign(config_data, key, value=forcetype(value, dtype))

        if key[0] != config_constant.INPUT_KEY:
            param_config_batch_parameter_dict[key[-1]] = value

        if isinstance(value, dict):
            value = list(value.keys())[0]
        root_dir = Path(root_dir) / f'{key[-1]}_{value}'

    return config_data, root_dir, param_config_batch_parameter_dict


def write_config(config_data: dict[str, Any], output_fname: str) -> None:
    """Main function for the script. This function will prepare and write the config provided by user 
    to run ScrambleBench.

    This function will perform the following:
    1) Re-assign the structure to follow the same template 
       (either from user supplying multiple protein targets or single)
    2) Find parameters that have multiple values (e.g., num_sample: 100, 200) and split them 
       into individual files with combinatorial combination
    3) Save individuals config as YAML
    4) Save location of such configs into txt file

    Args:
        config_data (dict[str, Any]): user's YAML config data as dictionary
        output_fname (str): output filename (not absolute directory, will be appended by prefix)
    """

    config_data = prepare_config(config_data=config_data)

    batch_parameters: list[dict[str, list]] = fetch_batch_parameters(config_data=config_data)

    parameter_value_lists = [parameter_dict['value'] for parameter_dict in batch_parameters]

    yaml_list = []
    generation_dirpath = Path(GenerationConfig(config_data).input_value).resolve()

    # [ [val, val, val], [val, val, val] ]
    for assigned_parameter_values in list(itertools.product(*parameter_value_lists)):

        result = process_single_batch_combination(root_dir=generation_dirpath,
                                                assigned_parameter_values=assigned_parameter_values,
                                                batch_parameters=batch_parameters,
                                                config_data=deepcopy(config_data))

        assigned_config_data, analysis_dirpath, param_config_batch_parameter_dict = result

        config_output = write_new_config(config_data=assigned_config_data,
                                         prefix_dir=analysis_dirpath,
                                         batch_parameter=param_config_batch_parameter_dict)

        Path(analysis_dirpath).mkdir(parents=True, exist_ok=True)
        yaml_fname = Path(analysis_dirpath) / Path(output_fname).name
        with open(yaml_fname, 'w', encoding='utf-8') as yaml_f:
            yaml.dump(config_output, yaml_f, Dumper=MyDumper, sort_keys=False)
            yaml_list.append(str(yaml_fname))
            logging.info(f"Config saved in {yaml_fname}")

    logging.info("Preparing Script for Batch Generation")

    with open(Path(generation_dirpath) / 'yaml_list.txt', 'w', encoding='utf-8') as generation_fname:
        for yaml_f in yaml_list:
            generation_fname.write(f'{yaml_f} \n')

    logging.info(f"Yaml file list saved in {str(Path(generation_dirpath) / 'yaml_list.txt')} \
                 for job manager and p2_execute_generation.py")


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Prepare config file for ScrambleBench")

    parser.add_argument("-i", "--input", help="config yaml input file", type=str)
    parser.add_argument("-o", "--output", help="config yaml output file prefix (i.e., not path directory)", type=str)
    parser.add_argument("--dirpath_input", action='store_true', help='write input key as a single directory', default=False)
    args = parser.parse_args()

    logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='%(asctime)s - %(module)s: - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

    logging.info('Running p1_generate_config.py')
    logging.info('Reading the config filename :)')

    with open(args.input, 'r', encoding='utf-8') as input_f:
        data_input = yaml.safe_load(input_f)

    if args.dirpath_input:
        logging.info('Dirpath_input mode: expecting to read directories given a directory')
        from scramblebench.script.config_preparation.config_input import InputDirpathConfig as InputUserConfig
    else:
        from scramblebench.script.config_preparation.config_input import InputNonDirpathConfig as InputUserConfig

    if not args.output:
        args.output = Path(args.input).parent / f'{Path(args.input).stem}_clean_config.yml'
    if Path(args.output).suffix not in ['.yaml', '.yml']:
        logging.exception('Error failed in config output filename format')
        raise ValueError(f'{args.output} ends with {Path(args.output).suffix}. \
                         Only .yaml and .yml extension is allowed')

    args.output = str(Path(args.output))

    if validate_config(config_data=data_input):
        write_config(data_input, args.output)
