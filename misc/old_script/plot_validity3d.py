import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import json
from glob import glob
from rdkit import Chem
import ptitprince as pt
import seaborn as sns
from pathlib import Path
import re
from collections import defaultdict

COLORBLIND_PALETTE = ["#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#999999"]


def load_json_file(ligand_directory, extra_query_string=''):
    file_list = glob(f'{ligand_directory}/json_output/*{extra_query_string}*.json')
    ordered_file_list = []
    for ordered_tools in ['Pocket2Mol', 'Lingo3DMol', 'PocketFlow', 'PMDM', 'DiffSBDD', 'Chem42']:
        for lig_json_file in file_list:
            if ordered_tools in lig_json_file:
                ordered_file_list.append(lig_json_file)
                break
            
    #manually create a dict one by one in tools naming
    generated_ligand_dict = {}
    for lig_json_file in ordered_file_list:
        lig_dict = json.load(open(lig_json_file))
        if 'PMDM' in lig_json_file:
            generated_ligand_dict['PMDM'] = lig_dict
        elif 'DiffSBDD' in lig_json_file:
            generated_ligand_dict['DiffSBDD'] = lig_dict
        elif 'Pocket2Mol' in lig_json_file:
            generated_ligand_dict['Pocket2Mol'] = lig_dict
        elif 'PocketFlow' in lig_json_file:
            generated_ligand_dict['PocketFlow'] = lig_dict
        elif 'Lingo3DMol' in lig_json_file:
            generated_ligand_dict['Lingo3DMol'] = lig_dict
        elif 'Chem42' in lig_json_file:
            generated_ligand_dict['Chem42'] = lig_dict
        else:
            ValueError('No tools implemented or wrong files')

    return generated_ligand_dict


def plot_double_validity3d(generated_ligand_dict, generated_ligand_dict2, output_file, fig, title='', protein_name=[1,2]):
    generative_model1 = generated_ligand_dict.keys()
    generative_model2 = generated_ligand_dict2.keys()

    validity_3d = {}

    validity_3d[protein_name[0]] = [[round(generated_ligand_dict[tools]['Validity3D'],2) for tools in generative_model1], generative_model1]
    validity_3d[protein_name[1]] = [[round(generated_ligand_dict2[tools]['Validity3D'],2) for tools in generative_model2], generative_model2]

    validity_3d_df = pd.DataFrame()
    for protein in validity_3d.keys():
        validity_3d_value = validity_3d[protein][0]

        validity_3d_df1 = pd.DataFrame(validity_3d_value, columns=['Validity3D'])
        validity_3d_df1['generative_model'] = validity_3d[protein][1]
        validity_3d_df1['protein_name'] = protein

        if validity_3d_df.empty:
            validity_3d_df = validity_3d_df1
        else:
            validity_3d_df = pd.concat([validity_3d_df, validity_3d_df1])


    ax = fig.subplots(1)
    sns.barplot(x="generative_model", y="Validity3D", hue='protein_name',
                 data=validity_3d_df, palette=COLORBLIND_PALETTE)
    ax.set_xlabel('Generative Model',  weight = 'bold', fontsize=20)
    ax.set_ylabel('Validity3D',  weight = 'bold', fontsize=20)
    ax.tick_params(axis='both', which='major', labelsize=20)
    ax.set_title(title, fontsize=20)
    handles, labels = ax.get_legend_handles_labels()
    ax.get_legend().remove()
    ax.set_ylim(0,1)
    ax.legend(handles[:6], labels[:6], ncol=1, fontsize=20, framealpha=0.5)

    return fig


def create_case_insensitive_regex(pattern: str) -> str:
    return f"{''.join([ '[' + char.upper() + char.lower() + ']' for char in pattern])}"

root_dir = '/opt/veincent/GenAI_manuscript/output_test_multiple_numsample'
protein_family_dict = {'GPCR':['GPCR_5HT2C', 'GPCR_DRD2'],
                       'Kinase': ['Kinase_CDK2', 'Kinase_GSK3b'],
                       'Hydrolase': ['Hydrolase_3CLpro', 'Hydrolase_NA']}
num_sample = 500
genbench_dirname = 'cheminformatics_analysis'

forcefield_minimisation_list = ['_unminimised', '_minimised']

fig = plt.figure(layout='constrained', figsize=(25, 25))
subfigs = fig.subfigures(3, 2, wspace=0.07)

idx_x = 0
idx_y = 0
for protein_family in protein_family_dict.keys():
    idx_y = 0
    for force_field_minimisation in forcefield_minimisation_list:
        json_list = []
        protein_list = protein_family_dict[protein_family]

        for protein_name in protein_list:
            root_dirpath: Path = Path(root_dir)
            if not root_dirpath.exists():
                raise ValueError('Genbench root path is missing')
            if not root_dirpath.is_dir():
                raise TypeError(f'Genbench root path {root_dirpath} is not a directory')

            glob_matching_dirname = f'*{create_case_insensitive_regex(protein_name)}*{num_sample}*'
            regex_pattern_non_alphanumeric_left_and_right = re.compile(f'(?<![A-Za-z0-9]){create_case_insensitive_regex(protein_name)}(?![A-Za-z0-9])')

            matched_dirname_list = [fname for fname in root_dirpath.glob(glob_matching_dirname) if regex_pattern_non_alphanumeric_left_and_right.search(fname.name)]

            docking_dict: dict[tuple[str, str], dict[str, object]] = defaultdict(dict)
            for protein_dirpath in matched_dirname_list:

                if not protein_dirpath.is_dir():
                    raise ValueError(f'{protein_dirpath} is not a directory')
                
                genbench_folder = protein_dirpath.joinpath(genbench_dirname)
                if not genbench_folder.is_dir():
                    raise ValueError(f'{genbench_folder} is not a directory')
                
                json_data = load_json_file(str(genbench_folder), force_field_minimisation)

                json_list.append(json_data)

        subfigs[(idx_x, idx_y)] = plot_double_validity3d(json_list[0], 
                                                         json_list[1], 
                                                         f'validity3d_{protein_family}{force_field_minimisation}.png',  
                                                         subfigs[(idx_x, idx_y)],  
                                                         f'Validity3D of Generative Model on {protein_family} protein target (MMFF94 {force_field_minimisation.replace("_", "")})', 
                                                         protein_name=protein_list)
        idx_y += 1
    idx_x +=1

plt.savefig('validity_3d_hole.png')