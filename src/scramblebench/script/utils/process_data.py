from pathlib import Path

import logging
import re
from scramblebench.script.config_preparation import config_constant, config_parameter
from scramblebench.script.utils.error_handler import DirNotFoundError

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
            yaml_file = yaml_file.strip()
            yaml_filepath = Path(yaml_file).resolve()

            if not yaml_filepath.is_file():
                raise FileNotFoundError(f'The file {yaml_file} is not found. Please check your directory')
    
            if yaml_filepath.suffix.lower() not in ['.yaml', '.yml']:
                raise ValueError(f'Incorrect file {yaml_file}. We only support yaml file (i.e., .yaml or .yml). Please use the p1_generate_config.py to prepare your config file.')

            yaml_file_list.append(yaml_file)

        return yaml_file_list


def fetch_model_file_from_dir(dir_path, model_list) -> dict[str, str]:

    generation_folder_dirpath = Path(dir_path)
    if not generation_folder_dirpath.is_dir():
        logging.exception('You have prepared your molecule yet!')
        raise DirNotFoundError(f'{generation_folder_dirpath} is not found! Please make sure you run p3_prepare_molecule.py')
    
    valid_molecule_file_dict = {}
    for model in model_list:
        valid_molecule_file_dict[model] = find_file_name_through_regex(character=model, file_format='.sdf', dirname=Path(generation_folder_dirpath))
    
    return valid_molecule_file_dict


def fetch_model_file_from_model_dir(dir_path, model_list) -> dict[str, str]:

    generation_folder_dirpath = Path(dir_path)
    if not generation_folder_dirpath.is_dir():
        logging.exception('You have prepared your molecule yet!')
        raise DirNotFoundError(f'{generation_folder_dirpath} is not found! Please make sure you run p3_prepare_molecule.py')
    
    valid_molecule_file_dict = {}
    for model in model_list:

        model_dir = Path(generation_folder_dirpath) / model
        valid_molecule_file_dict[model] = find_file_name_through_regex(character=model, file_format='.sdf', dirname=Path(model_dir))
    
    return valid_molecule_file_dict


def fetch_model_plif_file_from_model_dir(dir_path, model_list) -> dict[str, str]:

    generation_folder_dirpath = Path(dir_path)
    if not generation_folder_dirpath.is_dir():
        logging.exception('You have prepared your molecule yet!')
        raise DirNotFoundError(f'{generation_folder_dirpath} is not found! Please make sure you run p3_prepare_molecule.py')
    
    valid_molecule_file_dict = {}
    for model in model_list:

        model_dir = Path(generation_folder_dirpath) / model / 'plif'
        valid_molecule_file_dict[model] = find_file_name_through_regex(character=model, file_format='.sdf', dirname=Path(model_dir))
    
    return valid_molecule_file_dict


def create_case_insensitive_regex(pattern: str) -> str:
    return f"{''.join([ '[' + char.upper() + char.lower() + ']' for char in pattern])}"


def find_file_name_through_regex(character, file_format, dirname):
    assert file_format[0] == '.'

    glob_matching_fname = f'*{create_case_insensitive_regex(character)}*{file_format}'
    regex_pattern_non_alphanumeric_left_and_right = re.compile(f'(?<![A-Za-z0-9]){create_case_insensitive_regex(character)}(?![A-Za-z0-9])')

    matched_fname_list = [fname for fname in dirname.glob(glob_matching_fname) if regex_pattern_non_alphanumeric_left_and_right.search(fname.name)]
    if len(matched_fname_list) > 1:
        logging.exception(f'We found more than 1 matching file for {character} characters: {matched_fname_list}. Please ensure only 1 is detected')
        raise ValueError(f'We found more than 1 matching file for {character} characters: {matched_fname_list}. Please ensure only 1 is detected')
    elif len(matched_fname_list) == 0:
        logging.warning(f'There are no matched file for {character} characters in {dirname}. Make sure this is intended')
    else:
        return str(matched_fname_list[0])