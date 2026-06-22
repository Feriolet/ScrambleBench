from typing import Any, Callable
from pathlib import Path

import yaml
import argparse
import logging
import sys
import subprocess
import json
from scramblebench.script.config_preparation import config_constant, config_parameter, config_diversity, config_genbench3d, config_redocking, config_post_generation, config_virtual_hit
from scramblebench.script.utils.error_handler import DirNotFoundError
from scramblebench.script.utils.process_data import read_input, find_file_name_through_regex, fetch_model_file_from_model_dir, fetch_model_plif_file_from_model_dir
from scramblebench.script.utils.process_mol import calculate_rms, calculate_physicochemical_properties, PhysicoChemicalProperties
import os
import pandas as pd
import tempfile
from copy import deepcopy
from collections import defaultdict
import rdkit
from rdkit import Chem
from enum import Enum, IntEnum
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

class GenBenchDockingMethod(Enum):
    VINA_MININPLACE = "Vina score"
    VINA_INPLACE = "Minimized Vina score"
    GLIDE_MININPLACE = "Glide score"
    GLIDE_INPLACE = "Minimized Glide score"

GENBENCH3D_VALIDITY_3D_STRING = 'Validity3D'

class GenBenchValidity3D(Enum):
    VALIDITY3D_BOND_LENGTH = 'Geometric mean bond q-value (Validity3D)'
    VALIDITY3D_BOND_ANGLE = 'Geometric mean angle q-value (Validity3D)'
    VALIDITY3D_BOND_TORSION = 'Geometric mean torsion q-value (Validity3D)'

def fetch_valid_genbench3d_json_file(genbench3d_data: config_genbench3d.GenBench3DConfig, model) -> dict[str, str]:

    output_dirpath = Path(genbench3d_data.output_value)
    if not output_dirpath.is_dir():
        logging.exception('You have prepared your molecule yet!')
        raise DirNotFoundError(f'{output_dirpath} is not found! Please make sure you run p3_prepare_molecule.py')
    
    valid_molecule_file_dict = {}

    complex_minimisation = ['unminimised']
    if genbench3d_data.do_complex_forcefield_minimisation_value:
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

    valid_molecule_file_dict = fetch_model_file_from_model_dir(dir_path=genbench3d_data.input_value,
                                                                    model_list=parameter_class.model_list_value)    

    genbench3d_dict = defaultdict(list)
    single_validity_dict = defaultdict(dict)
    for model, input_sdf_fname in valid_molecule_file_dict.items():
        input_mol = Chem.SDMolSupplier(input_sdf_fname, removeHs=False)
        name_l = [mol.GetProp('_Name') for mol in input_mol]
        genbench3d_dict['mol_id'] +=  name_l
        genbench3d_dict['Model'] += [model] * len(name_l)
        valid_genbench3d_json_dict = fetch_valid_genbench3d_json_file(genbench3d_data=genbench3d_data,
                                                    model=model)
                                                    
        for minimisation, json_fname in valid_genbench3d_json_dict.items():
            with open(json_fname, 'r') as json_f:
                genbench3d_json_data = json.load(json_f)

            genbench3d_docking_column = [GenBenchDockingMethod.VINA_INPLACE]
            if genbench3d_data.do_docking_forcefield_minimisation_value:
                genbench3d_docking_column += [GenBenchDockingMethod.VINA_MININPLACE]
            if genbench3d_data.schrodinger_dir_value:
                genbench3d_docking_column += [GenBenchDockingMethod.GLIDE_INPLACE]
            if genbench3d_data.schrodinger_dir_value and genbench3d_data.do_docking_forcefield_minimisation_value:
                genbench3d_docking_column += [GenBenchDockingMethod.GLIDE_MININPLACE]

            for genbench3d_metric in genbench3d_docking_column:
                if minimisation == 'minimised':
                    genbench3d_dict[f'FF_minimised_{genbench3d_metric.value}'] += genbench3d_json_data[genbench3d_metric.value]
                elif minimisation == 'unminimised':
                    genbench3d_dict[f'FF_unminimised_{genbench3d_metric.value}'] += genbench3d_json_data[genbench3d_metric.value]

            for validity in GenBenchValidity3D:
                if minimisation == 'minimised':
                    genbench3d_dict[f'FF_minimised_{validity.name}'] += genbench3d_json_data[validity.value]
                    single_validity_dict[model][f'minimised_Validity3D'] = genbench3d_json_data['Validity3D']
                elif minimisation == 'unminimised':
                    genbench3d_dict[f'FF_unminimised_{validity.name}'] += genbench3d_json_data[validity.value]
                    single_validity_dict[model][f'unminimised_Validity3D'] = genbench3d_json_data['Validity3D']


    return pd.DataFrame.from_dict(genbench3d_dict), single_validity_dict
        

def replace_best_score_with_plif(input_mol_list, plif_mol_list):

    new_mol_list = []

    input_mol_names = [mol.GetProp("_Name") for mol in input_mol_list]
    plif_mol_names = [mol.GetProp('_Name') for mol in plif_mol_list]

    for index, mol_name in enumerate(input_mol_names):
        if mol_name not in plif_mol_names:
            new_mol_list.append(input_mol_list[index])
        else:
            new_mol_list.append(plif_mol_list[plif_mol_names.index(mol_name)])

    return new_mol_list


def collect_redocking_glide_data(docking_data, schrodinger_dir, parameter_class):

    valid_glide_input_output_dict = defaultdict(dict)

    for model, input_fname in fetch_model_file_from_model_dir(dir_path=docking_data.input_value,
                                                                    model_list=parameter_class.model_list_value).items():
        valid_glide_input_output_dict[model]['input'] = input_fname

    for model, output_fname in fetch_model_file_from_model_dir(dir_path=docking_data.output_value,
                                                                    model_list=parameter_class.model_list_value).items():
        valid_glide_input_output_dict[model]['output'] = output_fname

    if docking_data.plif_value:
        for model, plif_fname in fetch_model_plif_file_from_model_dir(dir_path=docking_data.output_value,
                                                                        model_list=parameter_class.model_list_value).items():
            valid_glide_input_output_dict[model]['plif'] = plif_fname

    glide_dict = defaultdict(list)
    for model, glide_fnames in valid_glide_input_output_dict.items():

        input_fname = glide_fnames['input']
        output_fname = glide_fnames['output']

        if not Path(input_fname).suffix == '.sdf':
            raise ValueError(f'please convert your glide fname from {Path(input_fname.suffix)} to .sdf')
        if not Path(output_fname).suffix == '.sdf':
            raise ValueError(f'please convert your glide fname from {Path(output_fname.suffix)} to .sdf')
        
        mol_dict = defaultdict(dict)
        input_mol_l = [mol for mol in Chem.SDMolSupplier(input_fname, removeHs=False) if mol]
        for mol in input_mol_l:
            mol_dict[mol.GetProp("_Name")]['input'] = mol

        with tempfile.TemporaryDirectory() as tempfile_dir:
            best_score_sdf_fname = Path(tempfile_dir) / f'{Path(output_fname).stem}_best_score.sdf'
            cmd = [f'{schrodinger_dir}/utilities/glide_sort', '-best_by_title', 
                '-o', str(best_score_sdf_fname), output_fname]
            
            subprocess.run(cmd, text=True)

            output_mol_l = [mol for mol in Chem.SDMolSupplier(best_score_sdf_fname, removeHs=False) if mol]

            if docking_data.plif_value:
                try:
                    plif_mol_l =  [mol for mol in Chem.SDMolSupplier(glide_fnames['plif'], removeHs=False) if mol]
                    output_mol_l = replace_best_score_with_plif(output_mol_l, plif_mol_l)
                except OSError:
                    logging.warning(f'PLIF file exists {glide_fnames["plif"]} but is empty. Skipping integration of PLIF to Glide')
                except KeyError:
                    pass
                
        output_mol_name_l = [mol.GetProp('_Name') for mol in output_mol_l]
        for mol in output_mol_l:
            mol_dict[mol.GetProp("_Name")]['output'] = mol


        glide_dict['glide_redocking_rmsd'] += [calculate_rms(mol_dict[mol_name]['input'], mol_dict[mol_name]['output']) for mol_name in output_mol_name_l]
        glide_dict['mol_id'] += [mol.GetProp('_Name') for mol in output_mol_l]
        if docking_data.plif_value:
            glide_dict['glide_plif'] += [float(mol.GetProp('r_phase_PhaseScreenScore')) if mol.HasProp('r_phase_PhaseScreenScore') else 0 for mol in output_mol_l ]
        glide_dict['Model'] += [model] * len(output_mol_l)
        glide_dict['glide_redocking_score'] += [float(mol.GetProp('r_i_docking_score')) for mol in output_mol_l]

    return pd.DataFrame.from_dict(glide_dict)


def collect_redocking_easydock_data(docking_data, parameter_class):

    valid_easydock_input_output_dict = defaultdict(dict)

    for model, input_fname in fetch_model_file_from_model_dir(dir_path=docking_data.input_value,
                                                                    model_list=parameter_class.model_list_value).items():
        valid_easydock_input_output_dict[model]['input'] = input_fname

    for model, output_fname in fetch_model_file_from_model_dir(dir_path=docking_data.output_value,
                                                                    model_list=parameter_class.model_list_value).items():
        valid_easydock_input_output_dict[model]['output'] = output_fname

    if docking_data.plif_value:
        for model, plif_fname in fetch_model_plif_file_from_model_dir(dir_path=docking_data.output_value,
                                                                        model_list=parameter_class.model_list_value).items():
            valid_easydock_input_output_dict[model]['plif'] = plif_fname

    easydock_dict = defaultdict(list)
    for model, easydock_fnames in valid_easydock_input_output_dict.items():

        input_fname = easydock_fnames['input']
        output_fname = easydock_fnames['output']

        if not Path(input_fname).suffix == '.sdf':
            raise ValueError(f'please convert your glide fname from {Path(input_fname.suffix)} to .sdf')
        if not Path(output_fname).suffix == '.sdf':
            raise ValueError(f'please convert your glide fname from {Path(output_fname.suffix)} to .sdf')
        
        mol_dict = defaultdict(dict)
        input_mol_l = [mol for mol in Chem.SDMolSupplier(input_fname, removeHs=False) if mol]
        for mol in input_mol_l:
            mol_dict[mol.GetProp("_Name")]['input'] = mol

        output_mol_l = [mol for mol in Chem.SDMolSupplier(output_fname, removeHs=False) if mol]


        output_mol_name_l = [mol.GetProp('_Name') for mol in output_mol_l]
        for mol in output_mol_l:
            mol_dict[mol.GetProp("_Name")]['output'] = mol

        if docking_data.plif_value:
            try:
                plif_mol_l =  [mol for mol in Chem.SDMolSupplier(easydock_fnames['plif'], removeHs=False) if mol]
                output_mol_l = replace_best_score_with_plif(output_mol_l, plif_mol_l)
            except OSError:
                logging.warning(f'PLIF file exists {easydock_fnames["plif"]} but is empty. Skipping integration of PLIF to easydock')
            except KeyError:
                pass

        easydock_dict['easydock_redocking_rmsd'] += [calculate_rms(mol_dict[mol_name]['input'], mol_dict[mol_name]['output']) for mol_name in output_mol_name_l]
        easydock_dict['mol_id'] += [mol.GetProp('_Name') for mol in output_mol_l]
        if docking_data.plif_value:
            easydock_dict['easydock_plif'] += [float(mol.GetProp('plif_similarity')) if mol.HasProp('plif_similarity') else 0 for mol in output_mol_l ]

        easydock_dict['Model'] += [model] * len(output_mol_l)
        easydock_dict['easydock_redocking_score'] += [float(mol.GetProp('docking_score')) for mol in output_mol_l]

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


def collect_physicochemical_data(postgeneration_data: config_post_generation.PostGenerationConfig, parameter_class):

    physicochemical_dict = defaultdict(list)
    for model, input_fname in fetch_model_file_from_model_dir(dir_path=postgeneration_data.output_value,
                                                              model_list=parameter_class.model_list_value).items():
        if not Path(input_fname).suffix == '.sdf':
            raise ValueError(f'please convert your glide fname from {Path(input_fname.suffix)} to .sdf')


        input_mol_l = [mol for mol in Chem.SDMolSupplier(input_fname, removeHs=False) if mol]
        output_mol_l = [calculate_physicochemical_properties(mol) for mol in input_mol_l]

        physicochemical_dict['mol_id'] += [mol.GetProp('_Name') for mol in output_mol_l]
        physicochemical_dict['Model'] += [model] * len(output_mol_l)
        physicochemical_dict['SMILES'] += [Chem.MolToSmiles(mol) for mol in output_mol_l]
        physicochemical_dict[PhysicoChemicalProperties.MW.value] += [float(mol.GetProp(PhysicoChemicalProperties.MW.value)) for mol in output_mol_l]
        physicochemical_dict[PhysicoChemicalProperties.QED.value] += [float(mol.GetProp(PhysicoChemicalProperties.QED.value)) for mol in output_mol_l]
        physicochemical_dict[PhysicoChemicalProperties.LOGP.value] += [float(mol.GetProp(PhysicoChemicalProperties.LOGP.value)) for mol in output_mol_l]
        physicochemical_dict[PhysicoChemicalProperties.SA_SCORE.value] += [float(mol.GetProp(PhysicoChemicalProperties.SA_SCORE.value)) for mol in output_mol_l]
        
    return pd.DataFrame.from_dict(physicochemical_dict)


def combine_analysis_df_with_parameter(analysis_df, parameter_data):

    parameter_order = []
    analysis_columns = analysis_df.columns

    analysis_df[parameter_data.protein_name] = parameter_data.protein_value
    parameter_order.append(parameter_data.protein_name)
    for key, val in parameter_data.batch_parameter_dict.items():
        parameter_order.append(key)
        analysis_df[key] = val

    analysis_df = analysis_df.reindex(columns=parameter_order + list(analysis_columns))
    return analysis_df


def deep_merge(dict1, dict2):
    """Recursively merges dict2 into dict1 without mutating the inputs."""
    result = dict1.copy()  # Create a shallow copy to prevent modifying dict1
    for key, value in dict2.items():
        if key in result and isinstance(result[key], dict) and isinstance(value, dict):
            result[key] = deep_merge(result[key], value)
        else:
            result[key] = value
    return result


def collect_model_general_performance(json_fname):

    with open(json_fname) as json_f:

        return json.load(json_f)


def collect_diversity_data(analysis_data, parameter_class: config_parameter.ParameterConfig):

    diversity_data = config_diversity.DiversityConfig(analysis_data)
    diversity_dict = defaultdict(dict)
    for model in parameter_class.model_list_value:

        diversity_fname = Path(Path(diversity_data.output_value) / model / 'diversity_output.json')
        if diversity_fname.is_file():
            with open(diversity_fname) as diversity_f:
                diversity_output = json.load(diversity_f)
            diversity_dict[model]['diversity'] =diversity_output[model]
    

    return diversity_dict


def filter_catalogue(smi, catalog):
    mol = Chem.MolFromSmiles(smi)
    if catalog.HasMatch(mol):
        return 0
    
    return 1

def collect_virtual_hit(virtual_hit_data, analysis_df):

    virtual_hit_dict = defaultdict(dict)
    for model, df in analysis_df.groupby(['Model']):
        mol_length = len(df)

        has_been_filtered = False
        if virtual_hit_data.query_value:
            df = df.query(virtual_hit_data.query_value)

            virtual_hit_count = len(df)
            has_been_filtered = True

        if virtual_hit_data.filter_value:
            params_chembl = FilterCatalogParams()

            if virtual_hit_data.filter_value.lower() == 'pains':
                params_chembl.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
            else:
                raise ValueError(f'{virtual_hit_data.filter_value} is not supported')
            
            catalog_chembl = FilterCatalog(params_chembl)

            df['PAINS'] = df['SMILES'].apply(filter_catalogue, catalog=catalog_chembl)
            df = df[df['PAINS'] == 1]

            virtual_hit_count = len(df)
            has_been_filtered = True

        if not has_been_filtered:
            return None
        
        virtual_hit_dict[model[0]]['virtual_hit'] = {'rate': f'{(virtual_hit_count / mol_length):.2f}',
                                                      'name': df['mol_id'].tolist(),
                                                      'filter': virtual_hit_data.filter_value,
                                                      'query': virtual_hit_data.query_value}
    
    return virtual_hit_dict


def collect_analysis_metric(config_data):

    analysis_data = config_data[config_constant.ANALYSIS_KEY]
    parameter_data = config_parameter.ParameterConfig(config_data=config_data)
    postgeneration_data = config_post_generation.PostGenerationConfig(config_data=config_data)

    
    config_dir = Path(postgeneration_data.input_value).parent
    analysis_df = pd.DataFrame()

    model_dict = collect_model_general_performance(config_dir / 'data.json')

    if config_constant.ANALYSIS_GENBENCH3D_KEY in analysis_data:
        genbench3d_df, validity3d_dict = collect_genbench3d_data(analysis_data=analysis_data, parameter_class = parameter_data)
        if analysis_df.empty:
            analysis_df = genbench3d_df
        
        model_dict = deep_merge(model_dict, validity3d_dict)


    if config_constant.ANALYSIS_REDOCKING_KEY in analysis_data:
        if config_constant.ANALYSIS_REDOCKING_DOCKING_KEY in analysis_data[config_constant.ANALYSIS_REDOCKING_KEY]:
            redocking_df = collect_redocking_data(analysis_data=analysis_data, parameter_class = parameter_data)
            if analysis_df.empty:
                analysis_df = redocking_df
            else:
                analysis_df = pd.merge(analysis_df, redocking_df, on=['mol_id', 'Model'], how='outer')


    physicochemical_df = collect_physicochemical_data(postgeneration_data=postgeneration_data, parameter_class=parameter_data)
    if analysis_df.empty:
        analysis_df = physicochemical_df
    else:
        analysis_df = pd.merge(analysis_df, physicochemical_df, on=['mol_id', 'Model'], how='outer')


    if config_constant.ANALYSIS_DIVERSITY_KEY in analysis_data:

        diversity_dict = collect_diversity_data(analysis_data=analysis_data, parameter_class = parameter_data)
        model_dict = deep_merge(model_dict, diversity_dict)

    analysis_df = combine_analysis_df_with_parameter(analysis_df, parameter_data=parameter_data)
    
    if config_constant.ANALYSIS_VIRTUAL_HIT_KEY in data_input[config_constant.ANALYSIS_KEY].keys():
        virtual_hit_dict = collect_virtual_hit(config_virtual_hit.VirtualHitConfig(analysis_data), analysis_df)
        if virtual_hit_dict:
            model_dict = deep_merge(model_dict, virtual_hit_dict)

    analysis_df.to_csv(config_dir / 'summary.csv')

    with open(config_dir / 'summary.json', 'w') as json_f:
        json.dump(model_dict, json_f, indent=4)

    return analysis_df, model_dict


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
    output_json_fname = Path(args.input).parent / 'compiled_result.json'
    output_df = pd.DataFrame()

    for yaml_file in yaml_file_list:
        checkpoint_json = Path(yaml_file).parent
        data_input = yaml.safe_load(open(yaml_file, 'r'))

        analysis_df, model_dict = collect_analysis_metric(data_input) 



        if output_df.empty:
            output_df = analysis_df
        else:
            output_df = pd.concat([output_df, analysis_df])
    

    output_df.to_csv(output_analysis_fname)