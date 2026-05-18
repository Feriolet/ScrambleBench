from typing import Any, Callable
from pathlib import Path

import yaml
import argparse
import logging
import sys
import subprocess
import os
import re
from scramblebench.script.config_preparation import config_constant, config_post_generation, config_parameter
from scramblebench.script.utils.error_handler import DirNotFoundError
from scramblebench.script.utils.process_data import read_input, fetch_model_folder_name, find_file_name_through_regex
import rdkit
from rdkit import Chem
import numpy as np
from copy import deepcopy
from collections import defaultdict
import json
logger = logging.getLogger(__name__)


# def fetch_valid_generated_molecule_file(config_data) -> dict[str, str]:
#     model_for_generation_list = fetch_model_folder_name(config_data=config_data)

#     GENERATION_SUMMARY_FOLDER = 'summary'
#     generation_folder_dirpath = Path(config_data[config_constant.POST_GENERATION_KEY]['input']) / GENERATION_SUMMARY_FOLDER
#     if not generation_folder_dirpath.is_dir():
#         logging.exception('You have not run any generation yet!')
#         raise DirNotFoundError(f'{generation_folder_dirpath} is not found! Please make sure you run p2_execute_generation.py')
    
#     valid_molecule_file_dict = {}
#     for model in model_for_generation_list:

#         matched_fname_list = find_file_name_through_regex(character=model, file_format='.sdf', dirname=generation_folder_dirpath)
#         if len(matched_fname_list) > 1:
#             logging.exception(f'We found more than 1 matching file for {model} model: {matched_fname_list}. Please ensure only 1 is detected')
#             raise ValueError(f'We found more than 1 matching file for {model} model: {matched_fname_list}. Please ensure only 1 is detected')
#         elif len(matched_fname_list) == 0:
#             logging.warning(f'There are no matched file for {model} model in {generation_folder_dirpath}. Make sure this is intended')

#         valid_molecule_file_dict[model] = str(matched_fname_list[0])
    
#     return valid_molecule_file_dict


def fetch_valid_generated_molecule_file(dir_path, model_list) -> dict[str, str]:

    generation_folder_dirpath = Path(dir_path)
    if not generation_folder_dirpath.is_dir():
        logging.exception('You have prepared your molecule yet!')
        raise DirNotFoundError(f'{generation_folder_dirpath} is not found! Please make sure you run p3_prepare_molecule.py')
    
    valid_molecule_file_dict = {}
    for model in model_list:

        matched_fname_list = find_file_name_through_regex(character=model, file_format='.sdf', dirname=Path(generation_folder_dirpath))
        if len(matched_fname_list) > 1:
            logging.exception(f'We found more than 1 matching file for {model} model: {matched_fname_list}. Please ensure only 1 is detected')
            raise ValueError(f'We found more than 1 matching file for {model} model: {matched_fname_list}. Please ensure only 1 is detected')
        elif len(matched_fname_list) == 0:
            logging.warning(f'There are no matched file for {model} model in {generation_folder_dirpath}. Make sure this is intended')
        else:
            valid_molecule_file_dict[model] = str(matched_fname_list[0])
    
    return valid_molecule_file_dict

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


def neutralize_atoms(mol):
    pattern = Chem.MolFromSmarts("[+1!h0!$([*]~[-1,-2,-3,-4]),-1!$([*]~[+1,+2,+3,+4])]")
    at_matches = mol.GetSubstructMatches(pattern)
    at_matches_list = [y[0] for y in at_matches]
    if len(at_matches_list) > 0:
        for at_idx in at_matches_list:
            atom = mol.GetAtomWithIdx(at_idx)
            chg = atom.GetFormalCharge()
            hcount = atom.GetTotalNumHs()
            atom.SetFormalCharge(0)
            atom.SetNumExplicitHs(hcount - chg)
            atom.UpdatePropertyCache()
    return mol


def compute_uniqueness_percentage(mol_l) -> float:
    total_mol = len(mol_l)
    temp_mol_l = [deepcopy(mol) for mol in mol_l]
    temp_mol_l = validate_mol_list(temp_mol_l)

    smi_l = [Chem.MolToSmiles(neutralize_atoms(mol)) for mol in temp_mol_l]
    
    return len(set(smi_l)) / total_mol


def compute_validity2d_percentage(mol_l) -> float:
    total_mol = len(list(mol_l))
    temp_mol_l = [deepcopy(mol) for mol in mol_l]
    temp_mol_l = validate_mol_list(temp_mol_l)

    return len(temp_mol_l) / total_mol


def process_mol(mol_l: list[Chem.Mol], 
                model_name: str,
                config_data: dict[str, Any]) -> list[Chem.Mol]:
    """Check if the generated ligand is below or equal to num_sample parameter (default 100).
    Invalid ligands will be removed, and announced in the print output.

    If Pocket2Mol and Lingo3DMol ligand > num_sample, pick the last num_sample molecules.
    If Chem42 ligand > num_sample, pick random molecule.

    Side note: I am not sure why Lingo3DMol always have dupes, maybe
    there's an error with the saving file at the last line or smth.

    Args:
        sdf_file (str): generated sdf filename. In this context, this should be the 'summary' folder
        config_data (dict[Any]): The config data parsed by yaml.
    """

    num_sample = config_data[config_constant.PARAMETER_KEY]['num_sample']
    repeat_parameter_dict = config_data[config_constant.PARAMETER_KEY]['repeat_parameter']
    model_pick_last_str = config_data[config_constant.POST_GENERATION_KEY].get('pick_last')
    model_pick_random_str = config_data[config_constant.POST_GENERATION_KEY].get('pick_random')

    if len(mol_l) > num_sample:
        logging.info(f'{model_name} model has {len(mol_l)} which exceed {num_sample}')
        if model_pick_last_str and model_name.lower() in [model.lower() for model in model_pick_last_str.split(',')]:
            logging.info(f'trimming {model_name} by picking last {num_sample} molecules')
            mol_l = mol_l[len(mol_l)-num_sample:]

        elif model_pick_random_str and model_name.lower() in [model.lower() for model in model_pick_random_str.split(',')]:
            logging.info(f'trimming {model_name} by picking randomly')
            mol_l = np.random.choice(mol_l, num_sample, replace=False).tolist()
        else:
            logging.warning(f'The model {model_name} has ligand more than {num_sample} and you have \
                            not specified how you want to filter it through the model_pick_last_l or model_pick_random_l. \
                            Picking random ligand instead')
            mol_l = np.random.choice(mol_l, num_sample, replace=False).tolist()
    
    assert len(mol_l) <= num_sample

    protein_name = config_data[config_constant.PARAMETER_KEY]['protein_name']

    for index, m in enumerate(mol_l):
        m.SetProp('GenAI_Model', model_name)
        m.SetProp('Protein', protein_name)

        additional_name_string = ''
        if repeat_parameter_dict:
            for key, key_config_path in repeat_parameter_dict.items():
                value = deep_get(config_data, key_config_path.split(','))
                additional_name_string += f'{key}_{value}'
                m.SetProp(key, str(value))

        if additional_name_string:
            m.SetProp('_Name', f'{model_name}_{index}_{additional_name_string}')
        else:
            m.SetProp('_Name', f'{model_name}_{index}')

    return mol_l


def validate_mol_list(mol_list: list[Chem.Mol]) -> list[Chem.Mol]:
    # only filter valid and unique molecule
    validated_mol_l = []
    validated_smi_l = []
    for mol in mol_list:
        if not mol:
            continue

        smi = Chem.MolToSmiles(neutralize_atoms(mol))

        if smi != '' and mol.GetNumAtoms() > 1 and smi not in validated_smi_l:
            validated_mol_l.append(mol)
            validated_smi_l.append(smi)
    
    return validated_mol_l


def prepare_molecule(config_data):

    post_generation_config_data = config_data[config_constant.POST_GENERATION_KEY]
    parameter_data = config_parameter.ParameterConfig(config_data=config_data)
    postgen_data = config_post_generation.PostGenerationConfig(config_data)
    valid_molecule_file_dict = fetch_valid_generated_molecule_file(dir_path=str(Path(postgen_data.input_value) / 'summary'),
                                                                   model_list=parameter_data.model_list_value)
    post_generation_output_root_dirpath = post_generation_config_data['output']
    json_dir = Path(post_generation_output_root_dirpath).parent
    json_content = defaultdict(dict)
    for model, fname in valid_molecule_file_dict.items():
        logging.info(f'log {fname} {model}')
        file_name = Path(fname).stem + '_prepared.sdf'
        output_dir = Path(post_generation_output_root_dirpath) / model
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / file_name
        mol_l = Chem.SDMolSupplier(fname)
        uniqueness = compute_uniqueness_percentage(mol_l)
        json_content[model]['total_mol'] = len(mol_l)
        json_content[model]['uniqueness'] = uniqueness
        validity = compute_validity2d_percentage(mol_l=mol_l)
        json_content[model]['validity2d'] = validity
        print(f'{uniqueness=}, {validity=}')
        mol_l = process_mol(validate_mol_list(mol_l), model, config_data=config_data)

        with Chem.SDWriter(output_file) as w:
            for mol in mol_l:
                w.write(mol)

        logging.info(f'Prepared molecule saved in {output_file}')

    with open(Path(json_dir) / 'data.json', 'w') as file:
        json.dump(json_content, file, indent=4)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Prepare generated molecules for downstream analysis")

    parser.add_argument("-i", "--input", help="config yaml input file or txt file containing yaml filepath", required=True, type=str)
    args = parser.parse_args()

    logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='%(asctime)s - %(module)s: - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
    
    logging.info('Running p3_prepare_molecule.py')
    logging.info('Reading the config filename :)')

    yaml_file_list = read_input(args.input)

    for yaml_file in yaml_file_list:
        data_input = yaml.safe_load(open(yaml_file, 'r'))
        prepare_molecule(config_data=data_input)