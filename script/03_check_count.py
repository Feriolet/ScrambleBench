import argparse
import numpy as np
import os
import yaml

from glob import glob
from multiprocessing import Pool
from pathlib import Path
from rdkit import Chem
from typing import Any

def check_count(sdf_file: str, 
                config_data: dict[str, Any], 
                output_dirname: str = 'cheminformatics_input_prepared',
                model_name_l: list[str] = ['Pocket2Mol', 'PocketFlow', 'Lingo3DMol', 'PMDM', 'DiffSBDD', 'Chem42'],
                model_pick_last_l: list[str]=['Lingo3DMol', 'Pocket2Mol'],
                model_pick_random_l: list[str]=['Chem42']) -> None:
    """Check if the generated ligand is below or equal to num_sample parameter (default 100).
    Invalid ligands will be removed, and announced in the print output.

    If Pocket2Mol and Lingo3DMol ligand > num_sample, pick the last num_sample molecules.
    If Chem42 ligand > num_sample, pick random molecule.

    Side note: I am not sure why Lingo3DMol always have dupes, maybe
    there's an error with the saving file at the last line or smth.

    Args:
        sdf_file (str): generated sdf filename. In this context, this should be the 'summary' folder
        config_data (dict[Any]): The config data parsed by yaml.
        output_dirname (str, optional): the name of the output after filtering the ligands. Defaults to 'cheminformatics_input_prepared'.
        model_name_l (list[str], optional): List of model name used to generate ligand. Defaults to ['Pocket2Mol', 'PocketFlow', 'Lingo3DMol', 'PMDM', 'DiffSBDD', 'Chem42'].
        model_pick_last_l (list[str], optional): List of model name used to pick last num_sample ligands as filter. Defaults to ['Lingo3DMol', 'Pocket2Mol'].
        model_pick_random_l (list[str], optional): List of model name used to pick random num_sample ligands as filter. Defaults to ['Chem42'].
    """

    sdf_directory_split = sdf_file.rsplit('/',2)
    Path(os.path.join(sdf_directory_split[0], output_dirname)).mkdir(parents=True, exist_ok=True)
    output_file = os.path.join(sdf_directory_split[0], output_dirname, sdf_directory_split[-1])
    mol_l = [m for m in Chem.SDMolSupplier(sdf_file) if m is not None]
    num_sample = int(config_data['parameter']['num_sample'])

    if len(mol_l) > num_sample:
        if any([model_name.lower() in sdf_file.lower() for model_name in model_pick_last_l]):
            mol_l = [x for x in mol_l]
            mol_l = mol_l[len(mol_l)-num_sample:]
        elif any([model_name.lower() in sdf_file.lower() for model_name in model_pick_random_l]):
            mol_l = np.random.choice(mol_l, num_sample, replace=False).tolist()
        else:
            raise TypeError(f'The file {sdf_file} has ligand more than {num_sample} and you have \
                            not specified how you want to filter it through the model_pick_last_l or model_pick_random_l')

    assert len(mol_l) <= num_sample

    for model in model_name_l:
        if model in sdf_file:
            model_name = model
            break

    with Chem.SDWriter(output_file) as w:
        for index, m in enumerate(mol_l):
            m.SetProp('_Name', f'{model_name}_num{num_sample}_{index}')
            m.SetProp('GenAI_Model', model_name)
            m.SetProp('num_sample', num_sample)
            m.SetProp('Protein', config_data['parameter']['protein_title'])
            w.write(m)

def main():

    parser = argparse.ArgumentParser(description="Standardise generated compound output to create at most num_sample (default 100) ligand per generated compound.\n\n" \
                                     "This step is done after generating all of the Gen AI model. For most cases, you should use " \
                                     "the 'summary' folder for each GenAI model performed in step 01_run_all.sh. " \
                                     "If each model generates > num_sample ligand, this script will either filter it through picking the last compound or picking randomly," \
                                     "\n\n" \
                                     "The model name should be listed in the filename", formatter_class=argparse.RawTextHelpFormatter )

    parser.add_argument("--input", "-i", type=str, required=True,
                        help="could be an individual sdf file or folder containing sdf file")
    parser.add_argument("--config_path", "-c", type=int, required=False, default=100,
                        help="path of the config file generated from step 00_prepare_config_file.py")
    parser.add_argument("--output_dirname", "-o", type=str, required=False, default='cheminformatics_input_prepared',
                        help="Name out directory as output.")
    parser.add_argument('--pick_last', nargs='+', type=str, required=False, help='Model name to filter excessive ligands by picking the last ligands. Default is Lingo3DMol and Pocket2Mol')
    parser.add_argument('--pick_random', nargs='+', type=str, required=False, help='Model name to filter excessive ligands by picking ligands randomly. Default is Chem42')
    parser.add_argument('--model_list', nargs='+', type=str, required=False, help='List all models used to generate ligands within the specified directory. Default is: ' \
                        'Pocket2Mol, Lingo3DMol, PocketFlow, DiffSBDD, PMDM, Chem42')
    
    args = parser.parse_args()

    with open(args.config_file, 'r') as config_f:
        config_data = yaml.safe_load(config_f)

    if args.pick_last:
        model_pick_last_l = args.pick_last
    else:
        model_pick_last_l = ['Pocket2Mol', 'Lingo3DMol']
    
    if args.pick_random:
        model_pick_random_l = args.pick_random
    else:
        model_pick_random_l = ['Chem42']
    
    if args.model_list:
        model_name_l = args.model_list
    else:
        model_name_l = ['Pocket2Mol', 'PocketFlow', 'Lingo3DMol', 'PMDM', 'DiffSBDD', 'Chem42']
    
    for model_pick_random in model_pick_random_l:
        if model_pick_random not in model_name_l:
            print(f'[WARNING]: You specified {model_pick_random} in args.pick_random, but this is not in args.model_name')
        for model_pick_last in model_pick_last_l:
            if model_pick_random.lower() == model_pick_last.lower():
                raise ValueError(f'You specified {model_pick_last} in both args.pick_last and args.pick_random. Terminating...')
            if model_pick_last not in model_name_l:
                print(f'[WARNING]: You specified {model_pick_last} in args.pick_last, but this is not in args.model_name')
    
    if os.path.isdir(args.input):
        print(f'This script detected input as directory')
        print(f'searching for sdf file...')
        sdf_filenames = glob(f'{args.input}/*.sdf')

        if len(sdf_filenames) == 0:
            raise TypeError(f'no sdf filename is detected, terminating...')
        else:
            for sdf_file in sdf_filenames:
                check_count(sdf_file, 
                            config_data=config_data,
                            output_dirname=args.output_dirname,
                            model_name_l=model_name_l,
                            model_pick_last_l=model_pick_last_l,
                            model_pick_random_l=model_pick_random_l)

    elif os.path.isfile(args.input):
        print('This script detected input as a file')
        check_count(args.input,
                    config_data=config_data,
                    output_dirname=args.output_dirname,
                    model_name_l=model_name_l,
                    model_pick_last_l=model_pick_last_l,
                    model_pick_random_l=model_pick_random_l)
    else:
        raise TypeError('Input filename must be an existing directory of file! Exiting...')


if __name__ == '__main__':
    main()