"""This file contains functions that are shared by similar Config classes and
functions that are used in p1_generate_config.py
"""

import subprocess
import logging


from copy import deepcopy
from typing import Any
from pathlib import Path

import yaml

logger = logging.getLogger(__name__)


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


def check_file_start_with_number(fname: str):
    """check if filename starts with a number. Having numbers at start may cause an issue in processing

    Args:
        fname (str): filename

    Raises:
        TypeError: if filename starts with a number
    """
    try:
        float(Path(fname).name[0])
        raise TypeError(f'please do not use number as starting filename. change {fname}')
    except ValueError:
        pass


def check_conda_env(conda_env: str) -> None:
    """check if conda environment is active and present

    Args:
        conda_env (str): name of conda environment

    Raises:
        ValueError: if conda environment is not in the 'conda env list'
    """

    # Run 'conda env list' command and capture the output
    output = subprocess.run(
        "conda env list | awk '{ print $1}'",
        shell=True,
        capture_output=True,
        text=True,
        check=True
    ).stdout.split('\n')

    if not conda_env in output:
        if conda_env == 'not_applicable' or conda_env is None:
            logging.warning('Please note that your conda is explicitly left out')
        else:
            raise ValueError(f'Conda environment {conda_env} does not exist')


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
