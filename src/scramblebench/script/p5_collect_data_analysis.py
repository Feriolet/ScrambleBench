from typing import Any, Callable
from pathlib import Path

import yaml
import argparse
import logging
import sys
import subprocess
import json
from scramblebench.script.config_preparation import config_constant, config_parameter, config_diversity, config_genbench3d, config_redocking
from scramblebench.script.utils.error_handler import DirNotFoundError
from scramblebench.script.utils.process_data import read_input, find_file_name_through_regex

import os
import pandas as pd
import tempfile
from copy import deepcopy
from collections import defaultdict
import rdkit
from rdkit import Chem
from enum import Enum, IntEnum
from rdkit.Chem import rdRascalMCES

class GenBenchDockingMethod(Enum):
    VINA_MININPLACE = "Vina score"
    VINA_INPLACE = "Minimized Vina score"
    GLIDE_MININPLACE = "Glide score"
    GLIDE_INPLACE = "Minimized Glide score"


def fetch_valid_prepared_molecule_file(dir_path, model_list) -> dict[str, str]:

    generation_folder_dirpath = Path(dir_path)
    if not generation_folder_dirpath.is_dir():
        logging.exception('You have prepared your molecule yet!')
        raise DirNotFoundError(f'{generation_folder_dirpath} is not found! Please make sure you run p3_prepare_molecule.py')
    
    valid_molecule_file_dict = {}
    for model in model_list:

        matched_fname_list = find_file_name_through_regex(character=model, file_format='.sdf', dirname=Path(generation_folder_dirpath) / model)
        if len(matched_fname_list) > 1:
            logging.exception(f'We found more than 1 matching file for {model} model: {matched_fname_list}. Please ensure only 1 is detected')
            raise ValueError(f'We found more than 1 matching file for {model} model: {matched_fname_list}. Please ensure only 1 is detected')
        elif len(matched_fname_list) == 0:
            logging.warning(f'There are no matched file for {model} model in {generation_folder_dirpath}. Make sure this is intended')
        else:
            valid_molecule_file_dict[model] = str(matched_fname_list[0])
    
    return valid_molecule_file_dict


def fetch_valid_genbench_json_file(genbench_data: config_genbench3d.GenBench3DConfig, model) -> dict[str, str]:

    
    output_dirpath = Path(genbench_data.output_value)
    if not output_dirpath.is_dir():
        logging.exception('You have prepared your molecule yet!')
        raise DirNotFoundError(f'{output_dirpath} is not found! Please make sure you run p3_prepare_molecule.py')
    
    valid_molecule_file_dict = {}

    complex_minimisation = ['unminimised']
    if genbench_data.do_complex_forcefield_minimisation_value:
        complex_minimisation += ['minimised']

    for minimisation in complex_minimisation:

        matched_fname_list = find_file_name_through_regex(character=model, file_format='.json', dirname=Path(output_dirpath) / 'json_output')
        matched_fname_list = [matched_fname for matched_fname in matched_fname_list if f'_{minimisation}' in str(matched_fname)]
        if len(matched_fname_list) > 1:
            logging.exception(f'We found more than 1 matching file for {model} model: {matched_fname_list}. Please ensure only 1 is detected')
            raise ValueError(f'We found more than 1 matching file for {model} model: {matched_fname_list}. Please ensure only 1 is detected')
        elif len(matched_fname_list) == 0:
            logging.warning(f'There are no matched file for {model} model in {output_dirpath}. Make sure this is intended')
        else:
            valid_molecule_file_dict[minimisation] = str(matched_fname_list[0])
    
    return valid_molecule_file_dict


def collect_genbench3d_data(analysis_data, parameter_class):
    genbench3d_data = config_genbench3d.GenBench3DConfig(analysis_data)

    valid_molecule_file_dict = fetch_valid_prepared_molecule_file(dir_path=genbench3d_data.input_value,
                                                                    model_list=parameter_class.model_list_value)    

    genbench_dict = defaultdict(list)
    for model, input_sdf_fname in valid_molecule_file_dict.items():
        input_mol = Chem.SDMolSupplier(input_sdf_fname)
        name_l = [mol.GetProp('_Name') for mol in input_mol]
        genbench_dict['mol_id'] +=  name_l
        genbench_dict['Model'] += [model] * len(name_l)
        valid_genbench_json_dict = fetch_valid_genbench_json_file(genbench_data=genbench3d_data,
                                                    model=model)
                                                    
        for minimisation, json_fname in valid_genbench_json_dict.items():
            with open(json_fname, 'r') as json_f:
                genbench3d_json_data = json.load(json_f)

            genbench_docking_column = [GenBenchDockingMethod.VINA_INPLACE]
            if genbench3d_data.do_docking_forcefield_minimisation_value:
                genbench_docking_column += [GenBenchDockingMethod.VINA_MININPLACE]
            if genbench3d_data.schrodinger_dir_value:
                genbench_docking_column += [GenBenchDockingMethod.GLIDE_INPLACE]
            if genbench3d_data.schrodinger_dir_value and genbench3d_data.do_docking_forcefield_minimisation_value:
                genbench_docking_column += [GenBenchDockingMethod.GLIDE_MININPLACE]

            for genbench_metric in genbench_docking_column:
                if minimisation == 'minimised':
                    genbench_dict[f'FF_minimised_{genbench_metric.value}'] += genbench3d_json_data[genbench_metric.value]
                elif minimisation == 'unminimised':
                    genbench_dict[f'FF_unminimised_{genbench_metric.value}'] += genbench3d_json_data[genbench_metric.value]


    return pd.DataFrame.from_dict(genbench_dict)
        

def collect_redocking_glide_data(docking_data, schrodinger_dir, parameter_class):

    valid_molecule_file_dict = fetch_valid_prepared_molecule_file(dir_path=docking_data.output_value,
                                                                    model_list=parameter_class.model_list_value)    
    
    glide_dict = defaultdict(list)
    for model, glide_fname in valid_molecule_file_dict.items():
        if not Path(glide_fname).suffix == '.sdf':
            raise ValueError(f'please convert your glide fname from {Path(glide_fname.suffix)} to .sdf')
        
        
        with tempfile.TemporaryDirectory() as tempfile_dir:
            best_score_sdf_fname = Path(tempfile_dir) / f'{Path(glide_fname).stem}_best_score.sdf'
            cmd = [f'{schrodinger_dir}/utilities/glide_sort', '-best_by_title', 
                '-o', str(best_score_sdf_fname), glide_fname]
            
            subprocess.run(cmd, text=True)

            mol_l = [mol for mol in Chem.SDMolSupplier(best_score_sdf_fname) if mol]
            glide_dict['mol_id'] += [mol.GetProp('_Name') for mol in mol_l]
            glide_dict['Model'] += [model] * len(mol_l)
            glide_dict['glide_redocking_score'] += [mol.GetProp('r_i_docking_score') for mol in mol_l]

    return pd.DataFrame.from_dict(glide_dict)


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


def calculate_rms(mol1, mol2):

    mol1 = neutralize_atoms(Chem.RemoveHs(mol1))
    mol2 = neutralize_atoms(Chem.RemoveHs(mol2))

    try:
        return Chem.rdMolAlign.CalcRMS(mol1, mol2)
    except RuntimeError:
        opts = rdRascalMCES.RascalOptions()
        opts.completeAromaticRings = False
        opts.timeout = 30
        opts.similarityThreshold = 0.1
        opts.maxBondMatchPairs = 2500
        res = rdRascalMCES.FindMCES(mol1, mol2, opts)
        if not res:
            opts.ignoreBondOrders = True
            res = rdRascalMCES.FindMCES(mol1, mol2, opts)
        matches = res[0].atomMatches()
        return Chem.rdMolAlign.CalcRMS(mol1, mol2, map=[matches])


def collect_redocking_easydock_data(docking_data, parameter_class):

    valid_easydock_input_output_dict = defaultdict(dict)

    for model, input_fname in fetch_valid_prepared_molecule_file(dir_path=docking_data.input_value,
                                                                    model_list=parameter_class.model_list_value).items():
        valid_easydock_input_output_dict[model]['input'] = input_fname

    for model, output_fname in fetch_valid_prepared_molecule_file(dir_path=docking_data.output_value,
                                                                    model_list=parameter_class.model_list_value).items():
        valid_easydock_input_output_dict[model]['output'] = output_fname

    easydock_dict = defaultdict(list)
    for model, easydock_fnames in valid_easydock_input_output_dict.items():

        input_fname = easydock_fnames['input']
        output_fname = easydock_fnames['output']

        if not Path(input_fname).suffix == '.sdf':
            raise ValueError(f'please convert your glide fname from {Path(input_fname.suffix)} to .sdf')
        if not Path(output_fname).suffix == '.sdf':
            raise ValueError(f'please convert your glide fname from {Path(output_fname.suffix)} to .sdf')
        
        mol_dict = defaultdict(dict)
        input_mol_l = [mol for mol in Chem.SDMolSupplier(input_fname) if mol]
        for mol in input_mol_l:
            mol_dict[mol.GetProp("_Name")]['input'] = mol

        output_mol_l = [mol for mol in Chem.SDMolSupplier(output_fname) if mol]
        output_mol_name_l = [mol.GetProp('_Name') for mol in output_mol_l]
        for mol in output_mol_l:
            mol_dict[mol.GetProp("_Name")]['output'] = mol

        easydock_dict['easydock_redocking_rmsd'] += [calculate_rms(mol_dict[mol_name]['input'], mol_dict[mol_name]['output']) for mol_name in output_mol_name_l]
        easydock_dict['mol_id'] += [mol.GetProp('_Name') for mol in output_mol_l]
        easydock_dict['Model'] += [model] * len(output_mol_l)
        easydock_dict['easydock_redocking_score'] += [mol.GetProp('docking_score') for mol in output_mol_l]

    return pd.DataFrame.from_dict(easydock_dict)


def collect_redocking_data(analysis_data, parameter_class):
    redocking_data = config_redocking.RedockingConfig(analysis_data)

    redocking_df = pd.DataFrame()
    for docking_data in redocking_data.docking_value.valid_key_list:
        if isinstance(docking_data, config_redocking.GlideConfig):
            glide_df = collect_redocking_glide_data(docking_data, docking_data.dir_value, parameter_class)
            if redocking_df.empty:
                redocking_df = glide_df
            else:
                redocking_df = pd.merge(redocking_df, glide_df, on=['mol_id', 'Model'], how='outer')
        elif isinstance(docking_data, config_redocking.EasyDockConfig):
            easydock_df = collect_redocking_easydock_data(docking_data, parameter_class)
            if redocking_df.empty:
                redocking_df = easydock_df
            else:
                redocking_df = pd.merge(redocking_df, easydock_df, on=['mol_id', 'Model'], how='outer')        

    return redocking_df



def combine_analysis_df_with_parameter(analysis_df, parameter_data: config_parameter.ParameterConfig):

    parameter_order = []
    analysis_columns = analysis_df.columns

    analysis_df[parameter_data.protein_name] = parameter_data.protein_value
    parameter_order.append(parameter_data.protein_name)
    for key, val in parameter_data.repeat_parameter_dict.items():
        parameter_order.append(key)
        analysis_df[key] = val

    analysis_df = analysis_df.reindex(columns=parameter_order + list(analysis_columns))
    return analysis_df


def collect_analysis_metric(config_data):

    analysis_data = config_data[config_constant.ANALYSIS_KEY]
    parameter_data = config_parameter.ParameterConfig(config_data=config_data)

    analysis_df = pd.DataFrame()
    if config_constant.ANALYSIS_GENBENCH3D_KEY in analysis_data:
        genbench_df = collect_genbench3d_data(analysis_data=analysis_data, parameter_class = parameter_data)
        if analysis_df.empty:
            analysis_df = genbench_df

    if config_constant.ANALYSIS_REDOCKING_KEY in analysis_data:
        if config_constant.ANALYSIS_REDOCKING_DOCKING_KEY in analysis_data[config_constant.ANALYSIS_REDOCKING_KEY]:
            redocking_df = collect_redocking_data(analysis_data=analysis_data, parameter_class = parameter_data)
            if analysis_df.empty:
                analysis_df = redocking_df
            else:
                analysis_df = pd.merge(analysis_df, redocking_df, on=['mol_id', 'Model'], how='outer')

    analysis_df = combine_analysis_df_with_parameter(analysis_df, parameter_data=parameter_data)
    
    analysis_df.to_csv(Path(config_genbench3d.GenBench3DConfig(analysis_data).input_value).parent / 'all.csv')

    return analysis_df

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Collect data of analysis done in p4_analyse.py")

    parser.add_argument("-i", "--input", help="config yaml input file or txt file containing yaml filepath", required=True, type=str)
    args = parser.parse_args()

    logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='%(asctime)s - %(module)s: - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
    
    logging.info('Running p5_collect_data_analysis.py')
    logging.info('Reading the config filename :)')

    yaml_file_list = read_input(args.input)

    output_analysis_fname = Path(args.input).parent / 'compiled_result.csv'
    output_df = pd.DataFrame()

    for yaml_file in yaml_file_list:
        checkpoint_json = Path(yaml_file).parent
        data_input = yaml.safe_load(open(yaml_file, 'r'))

        analysis_df = collect_analysis_metric(data_input) 
        if output_df.empty:
            output_df = analysis_df
        else:
            output_df = pd.concat([output_df, analysis_df])
    
    output_df.to_csv(output_analysis_fname)