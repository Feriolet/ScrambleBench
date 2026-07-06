"""This script perform files manipulation or search"""
import logging
import re

from pathlib import Path
from scramblebench.script.utils.error_handler import DirNotFoundError


def read_input(input_fname: str) -> list[str]:
    """Fetch .yaml, .yml, .txt, or .text files that are related to the ScrambleBench YAML config

    Args:
        input_fname (str): input filenames

    Raises:
        FileNotFoundError: if input_fname is not a file that exists
        FileNotFoundError: if suffixes are not supported
        ValueError: if suffix is not supported

    Returns:
        list[str]: list of yaml files
    """
    input_filepath = Path(input_fname)

    if not input_filepath.is_file():
        raise FileNotFoundError(f'The file {input_fname} is not found. Please check your directory')

    if input_filepath.suffix.lower() in ['.yaml', '.yml']:
        return [input_filepath.resolve()]

    if input_filepath.suffix.lower() not in ['.txt', '.text']:
        return []

    logging.debug(f'Input is in suffix {input_filepath.suffix}. Reading content to fetch yaml files.')
    yaml_file_list = []
    input_content = input_filepath.read_text(encoding='utf-8').splitlines()

    for yaml_file in input_content:
        yaml_file = yaml_file.strip()
        yaml_filepath = Path(yaml_file).resolve()

        if not yaml_filepath.is_file():
            raise FileNotFoundError(f'The file {yaml_file} is not found. Please check your directory')

        if yaml_filepath.suffix.lower() not in ['.yaml', '.yml']:
            raise ValueError(f'Incorrect file {yaml_file}.'
                                'We only support yaml file (i.e., .yaml or .yml).'
                                'Please use the p1_generate_config.py to prepare your config file.')

        yaml_file_list.append(yaml_file)

    return yaml_file_list


def fetch_model_file_from_dir(dir_path: str, model_list: list[str]) -> dict[str, str]:
    """Fetch SDF files from a folder containing the model name (case insensitive).
    For example,

    folderA
    |---- hello_world_model1.sdf (fetched if model1 is in the list)
    |____ hello_world_model2.sdf (fetched if model2 is in the list)
    |____ hello_world.sdf (not fetched)

    Args:
        dir_path (str): queried directory path
        model_list (list[str]): list of model names

    Raises:
        DirNotFoundError: if directory is not found

    Returns:
        dict[str, str]: dictionary of {model: sdf filename}
    """

    generation_folder_dirpath = Path(dir_path)
    if not generation_folder_dirpath.is_dir():
        logging.exception('You have prepared your molecule yet!')
        raise DirNotFoundError(f'{generation_folder_dirpath} is not found!'
                                'Please make sure you run p3_prepare_molecule.py')

    valid_molecule_file_dict = {}
    for model in model_list:
        valid_molecule_file_dict[model] = find_file_name_through_regex(character=model,
                                                                       file_format='.sdf',
                                                                       dirname=Path(generation_folder_dirpath))

    return valid_molecule_file_dict


def fetch_model_file_from_model_dir(dir_path: str, model_list: list[str]) -> dict[str, str]:
    """Fetch SDF files from a folder containing the model name as subfolder (case insensitive).
    For example,

    folderA
    |---- model1
    |       |---- hello_world_model1.sdf (fetched if model1 is in the list)
    |---- model2
    |       |---- hello_world_model2.sdf (fetched if model2 is in the list)
    |       |---- hello_world_model1.sdf (not fetched)
    |
    |---- hello_world_model1.sdf (not fetched)
    |____ hello_world.sdf (not fetched)

    Args:
        dir_path (str): queried directory path
        model_list (list[str]): list of model names

    Raises:
        DirNotFoundError: if directory is not found

    Returns:
        dict[str, str]: dictionary of {model: sdf filename}
    """
    generation_folder_dirpath = Path(dir_path)
    if not generation_folder_dirpath.is_dir():
        logging.exception('You have prepared your molecule yet!')
        raise DirNotFoundError(f'{generation_folder_dirpath} is not found!'
                                'Please make sure you run p3_prepare_molecule.py')

    valid_molecule_file_dict = {}
    for model in model_list:

        model_dir = Path(generation_folder_dirpath) / model
        valid_molecule_file_dict[model] = find_file_name_through_regex(character=model,
                                                                       file_format='.sdf',
                                                                       dirname=Path(model_dir))

    return valid_molecule_file_dict


def fetch_model_plif_file_from_model_dir(dir_path: str, model_list: list[str]) -> dict[str, str]:
    """Fetch SDF files from a folder containing the model name + plif as subfolder (case insensitive).
    For example,

    folderA
    |---- model1
    |       |---- hello_world_model2.sdf (not fetched)
    |       |---- plif
    |               |---- hello_world_model1.sdf (fetched if model1 is in the list)
    |
    |---- hello_world_model1.sdf (not fetched)
    |____ hello_world.sdf (not fetched)

    Args:
        dir_path (str): queried directory path
        model_list (list[str]): list of model names

    Raises:
        DirNotFoundError: if directory is not found

    Returns:
        dict[str, str]: dictionary of {model: sdf filename}
    """

    generation_folder_dirpath = Path(dir_path)
    if not generation_folder_dirpath.is_dir():
        logging.exception('You have prepared your molecule yet!')
        raise DirNotFoundError(f'{generation_folder_dirpath} is not found!'
                               'Please make sure you run p3_prepare_molecule.py')

    valid_molecule_file_dict = {}
    for model in model_list:

        model_dir = Path(generation_folder_dirpath) / model / 'plif'
        valid_molecule_file_dict[model] = find_file_name_through_regex(character=model,
                                                                       file_format='.sdf',
                                                                       dirname=Path(model_dir))

    return valid_molecule_file_dict


def create_case_insensitive_regex(pattern: str) -> str:
    """create the pattern where insensitive case can be used.

    Args:
        pattern (str): character to be regex'ed

    Returns:
        str: insensitive case as string
    """
    return f"{''.join([ '[' + char.upper() + char.lower() + ']' for char in pattern])}"


def find_file_name_through_regex(character: str, file_format: str, dirname: str) -> str:
    """Find file name through regex. There should only be one result.
    This script will also help to match whole word, since some model or protein name is too short.

    Args:
        character (str): character (mostly model names) to be regex'ed
        file_format (str): file format (typically .sdf)
        dirname (str): directory to be searched

    Raises:
        ValueError: if there is more than 1 matching file in the folder

    Returns:
        str: fetched filename
    """
    assert file_format[0] == '.'

    glob_matching_fname = f'*{create_case_insensitive_regex(character)}*{file_format}'
    regex_pattern_non_alphanumeric_left_and_right = re.compile(f'(?<![A-Za-z0-9]){create_case_insensitive_regex(character)}(?![A-Za-z0-9])')

    matched_fname_list = [fname for fname in dirname.glob(glob_matching_fname) \
                          if regex_pattern_non_alphanumeric_left_and_right.search(fname.name)]
    if len(matched_fname_list) > 1:
        logging.exception(f'We found more than 1 matching file for {character} characters:'
                          f'{matched_fname_list}. Please ensure only 1 is detected')
        raise ValueError(f'We found more than 1 matching file for {character} characters:'
                         f'{matched_fname_list}. Please ensure only 1 is detected')
    if len(matched_fname_list) == 0:
        logging.warning(f'There are no matched file for {character} characters in {dirname}.'
                         'Make sure this is intended')
        return ''

    return str(matched_fname_list[0])


def find_file_name_through_regex_multiple_match(character: str, file_format: str, dirname: str) -> list[str]:
    """Find file name through regex. There can be multiple results
    This script will also help to match whole word, since some model or protein name is too short.

    Args:
        character (str): character (mostly model names) to be regex'ed
        file_format (str): file format (typically .sdf)
        dirname (str): directory to be searched

    Returns:
        str: fetched filename
    """
    assert file_format[0] == '.'

    glob_matching_fname = f'*{create_case_insensitive_regex(character)}*{file_format}'
    regex_pattern_non_alphanumeric_left_and_right = re.compile(f'(?<![A-Za-z0-9]){create_case_insensitive_regex(character)}(?![A-Za-z0-9])')

    matched_fname_list = [fname for fname in dirname.glob(glob_matching_fname) \
                          if regex_pattern_non_alphanumeric_left_and_right.search(fname.name)]
    return matched_fname_list
