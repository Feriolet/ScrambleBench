
import yaml
import argparse
import logging
import sys
import numpy as np
import json
import rdkit

from rdkit import Chem
from typing import Any, Callable
from pathlib import Path
from copy import deepcopy
from collections import defaultdict

from scramblebench.script.config_preparation import config_post_generation, config_parameter
from scramblebench.script.utils.process_data import read_input, fetch_model_file_from_dir
from scramblebench.script.utils.process_mol import validate_mol_list, compute_generation_performance


logger = logging.getLogger(__name__)


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


def add_property_to_mol(mol_l, model_name, parameter_data: config_parameter.ParameterConfig):

    protein_name = parameter_data.protein_value
    batch_parameter = parameter_data.batch_parameter_dict
    for index, m in enumerate(mol_l):
        m.SetProp('GenAI_Model', model_name)
        m.SetProp('Protein', protein_name)

        additional_name_string = ''
        if batch_parameter:
            for key, value in batch_parameter.items():
                m.SetProp(key, str(value))

        if additional_name_string:
            m.SetProp('_Name', f'{model_name}_{index}_{additional_name_string}')
        else:
            m.SetProp('_Name', f'{model_name}_{index}')

    return mol_l


def prepare_mol(mol_l: list[Chem.Mol], 
                model_name: str,
                config_data: dict[str, Any]) -> list[Chem.Mol]:
    """Check if the generated ligand is below or equal to num_sample parameter (default 100).
    Invalid ligands will be removed, and announced in the print output.

    Args:
        sdf_file (str): generated sdf filename. In this context, this should be the 'summary' folder
        config_data (dict[Any]): The config data parsed by yaml.
    """

    postgeneration_data = config_post_generation.PostGenerationConfig(config_data=config_data)
    parameter_data = config_parameter.ParameterConfig(config_data=config_data)

    model_pick_last = postgeneration_data.pick_last_value or []
    model_pick_random = postgeneration_data.pick_random_value or []

    if model_pick_last:
        model_pick_last = [model.lower() for model in model_pick_last]
    if model_pick_random:
        model_pick_random = [model.lower() for model in model_pick_random]

    num_sample = parameter_data.num_sample_value
    if len(mol_l) > num_sample:
        logging.info(f'{model_name} model has {len(mol_l)} which exceed {num_sample}')

        if model_name.lower() in model_pick_last:
            logging.info(f'trimming {model_name} by picking last {num_sample} molecules')
            mol_l = mol_l[len(mol_l)-num_sample:]

        elif model_name.lower() in model_pick_random:
            logging.info(f'trimming {model_name} by picking randomly')
            mol_l = np.random.choice(mol_l, num_sample, replace=False).tolist()
            
        else:
            logging.warning(f'The model {model_name} has ligand more than {num_sample} and you have \
                            not specified how you want to filter it through the model_pick_last_l or model_pick_random_l. \
                            Picking random ligand instead')
            mol_l = np.random.choice(mol_l, num_sample, replace=False).tolist()
    
    return add_property_to_mol(mol_l, model_name=model_name, parameter_data=parameter_data)


def prepare_molecule(config_data):

    post_generation_data = config_post_generation.PostGenerationConfig(config_data)
    postgen_input_dir = f'{post_generation_data.input_value}/summary'
    postgen_output_root_dirpath = Path(post_generation_data.output_value)

    model_list = config_parameter.ParameterConfig(config_data=config_data).model_list_value
    valid_molecule_file_dict = fetch_model_file_from_dir(dir_path=postgen_input_dir, model_list=model_list)

    json_dir = postgen_output_root_dirpath.parent
    json_content = defaultdict(dict)
    for model, fname in valid_molecule_file_dict.items():

        mol_l = Chem.SDMolSupplier(fname, removeHs=False)

        json_content[model] = compute_generation_performance(mol_l=mol_l)

        mol_l = prepare_mol(validate_mol_list(mol_l), model, config_data=config_data)

        output_dir = postgen_output_root_dirpath / model
        output_dir.mkdir(parents=True, exist_ok=True)
        output_file = output_dir / f'{Path(fname).stem}_prepared.sdf'

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