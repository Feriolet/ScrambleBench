import matplotlib.pyplot as plt
import numpy as np
import os
import pandas as pd
import seaborn as sns

from tsp_diversity import diversity_all
from glob import glob
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from time import time
from typing import Union, Optional, Any

COLOR_PALLETE = ["#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#999999"]

def average_tanimoto_distance_ecfp(mol_list: list[Chem.Mol]) -> list[np.ndarray]:
    """Calculate Average Pairwise Tanimoto Distance.
    This metric is calculated by the following formula:

    Given a list M with mol m,
    Tanimoto Similarity       = TanimotoSimilarity(MorganFP(m_i), MorganFP(m_j))
    Tanimoto Distance         = 1 - Tanimoto Similarity
    Average Tanimoto Distance = Average(Tanimoto Distance)

    Args:
        mol_list (list[Chem.Mol]): List of Mol object from rdkit

    Returns:
        np.ndarray: numpy array of average tanimoto distance with length M
    """
    mfp_fpg = rdFingerprintGenerator.GetMorganGenerator(radius=2)
    average_tanimoto_dist = []
    mol_fp_list =  [mfp_fpg.GetFingerprint(mol) for mol in mol_list]
    for mol_fp_i in mol_fp_list:
        tanimoto_dist = []
        for mol_fp_j in mol_fp_list:
            if mol_fp_i == mol_fp_j:
                continue
            else:
                tanimoto_dist.append(1 - DataStructs.TanimotoSimilarity(mol_fp_i, mol_fp_j))

        tanimoto_dist_array = np.array(tanimoto_dist)
        average_tanimoto_dist.append(np.average(tanimoto_dist_array))

    return average_tanimoto_dist

def combine_all_generated_ligand(input_dir: str, 
                                 genai_models: list[str], 
                                 additional_property: Optional[list[list[Any]]]=None) -> list[Chem.Mol]:
    """Combine generated ligand from different model contained in the same folder into one list. 
    The ligands will be identified through property in the Chem.Mol object using SetProp() method. 

    Args:
        input_dir (str): The directory path containing multiple or one sdf file
        genai_models (list[str]): List of GenAI model. The str should match the one in the file.
        additional_property (Optional[list[list[Any]]], optional): _description_. Defaults to None.

    Returns:
        list[Chem.Mol]: _description_
    """

    cheminformatic_input_fnames = glob(f'{input_dir}/*.sdf')
    all_generated_ligand_list = []

    for model in genai_models:
        detected_sdf_file = 'Unknown'
        # Find which file belongs to the corresponding model. If no output, skipped
        for generated_input_fname in cheminformatic_input_fnames:
            if model in generated_input_fname:
                detected_sdf_file = generated_input_fname
                break

        if detected_sdf_file == 'Unknown':
            continue

        sdmol = Chem.SDMolSupplier(detected_sdf_file)
        mol_list = []

        for mol in sdmol:
            if mol:
                if mol.GetNumAtoms() > 0:
                    if additional_property:
                        for prop in additional_property:
                            assert len(prop) == 2
                            mol.SetProp(prop[0], str(prop[1]))
                    mol_list.append(mol)
       
        all_generated_ligand_list += mol_list

    return all_generated_ligand_list

if __name__ == "__main__":

    # Please use the python_tsp conda environment to run this
    # This script reads this kind of folder structure
    # input_dirname_dirname_dirname
    # |
    # |---protein_date_numsample
    # |     |
    # |     |---cheminformatics_input_prepared (this is your input_dirname variable)
    # |           |
    # |           |---output1.sdf
    # |           |---output2.sdf
    # |
    # |---protein_date_numsample
    #       |
    #       |---cheminformatics_input_prepared (this is your input_dirname variable)
    #             |
    #             |---output1.sdf
    #             |---output2.sdf    
    # Your input parameter
    genai_models = ['Pocket2Mol', 'PocketFlow', 'Lingo3DMol', 'DiffSBDD', 'PMDM', 'Chem42']
    protein_list = ['GPCR_5HT2C','GPCR_DRD2', 'Kinase_CDK2', 'Kinase_GSK3b', 'Hydrolase_3CLpro', 'Hydrolase_NA']
    input_dirname_dirname_dirname = '/opt/veincent/GenAI_manuscript/output_test_multiple_numsample'
    input_dirname = 'cheminformatics_input_prepared'
    output_suffix = 'SixProteinBenchmark_DiversityAnalysis_average_Tanimoto'
    is_calculate_hamdiv_mces = False
    is_calculate_hamdiv_ecfp = False
    is_calculate_generic_bm = False
    is_calculate_average_tanimoto_distance_ecfp = True
    
    # Start script
    # Tanimoto_array_type = list[np.ndarray]
    # Mixed_data_type = Union[int, str, float]
    # diversity_library_dict_as_list: list[list[Union[Tanimoto_array_type, Mixed_data_type]]] = []
    # diversity_library_df_column = ['Protein', 'Model', 'num_sample']
    # if is_calculate_hamdiv_mces:
    #     diversity_library_df_column += ['Hamiltonian Tanimoto Diversity based on MCES', 'Hamiltonian Tanimoto Diversity MCES Time'] 
    # if is_calculate_generic_bm:
    #     diversity_library_df_column += ['Number of Generic BM', 'Generic BM Time']
    # if is_calculate_average_tanimoto_distance_ecfp:
    #     diversity_library_df_column += ['ECFP Average Tanimoto Distance']
    # if is_calculate_hamdiv_ecfp:
    #     diversity_library_df_column += ['Hamiltonian Tanimoto Diversity based on ECFP', 'Hamiltonian Tanimoto Diversity ECFP Time']
    

    # # First loop read the top most directory (protein_date_numsample)
    # # Second loop read through each cheminformatics directory
    # # Third loop read each sdf file in the cheminformatics directory
    # for protein in protein_list:
    #     output_step03_dirname_list = glob(f'{input_dirname_dirname_dirname}/{protein}*/{input_dirname}')
    #     for output_step03_dirname in output_step03_dirname_list:
    #         mol_list = combine_all_generated_ligand(output_step03_dirname, genai_models)
    #         num_sample = os.path.dirname(output_step03_dirname).split('_')[-1]

    #         for model in genai_models:
    #             library_data: list[Union[Tanimoto_array_type, Mixed_data_type]] = []
    #             library_data.extend([protein, model, int(num_sample)])

    #             is_calculate_hamdiv_mces_model = is_calculate_hamdiv_mces
    #             is_calculate_hamdiv_ecfp_model = is_calculate_hamdiv_ecfp
    #             is_calculate_generic_bm_model = is_calculate_generic_bm
    #             is_calculate_average_tanimoto_distance_ecfp_model = is_calculate_average_tanimoto_distance_ecfp

    #             specific_genai_model_mol_list = [mol for mol in mol_list if mol.GetProp('GenAI_Model') == model]
    #             if len(specific_genai_model_mol_list) == 0:
    #                 is_calculate_hamdiv_mces_model = False
    #                 is_calculate_hamdiv_ecfp_model = False
    #                 is_calculate_generic_bm_model = False
    #                 is_calculate_average_tanimoto_distance_ecfp_model = False

    #             print(f'running {model} for {protein} (sample size {num_sample})...')
    #             if is_calculate_hamdiv_mces_model:
    #                 start_tanimoto_time = time()
    #                 tanimoto_hamdiv_mces = diversity_all(mols=specific_genai_model_mol_list, mode='HamDiv', hamdiv_method='MCES', ncpu=200)
    #                 tanimoto_hamdiv_mces_total_time = time() - start_tanimoto_time
    #                 print(f'calculation done in {tanimoto_hamdiv_mces_total_time}s ! HamDiv = {tanimoto_hamdiv_mces}')
    #                 library_data.extend([tanimoto_hamdiv_mces, tanimoto_hamdiv_mces_total_time])
    #             elif is_calculate_hamdiv_mces:
    #                 tanimoto_hamdiv_mces, tanimoto_hamdiv_mces_total_time = 0, 0
    #                 library_data.extend([tanimoto_hamdiv_mces, tanimoto_hamdiv_mces_total_time])

    #             if is_calculate_generic_bm_model:
    #                 start_generic_bm_time = time() 
    #                 generic_bm = diversity_all(mols=specific_genai_model_mol_list, mode='Generic_BM')
    #                 generic_bm_total_time = time() - start_generic_bm_time
    #                 print(f'calculation done in {generic_bm_total_time}s ! Generic BM = {generic_bm}')
    #                 library_data.extend([generic_bm, generic_bm_total_time])
    #             elif is_calculate_generic_bm:
    #                 generic_bm, generic_bm_total_time = 0, 0
    #                 library_data.extend([generic_bm, generic_bm_total_time])

    #             if is_calculate_average_tanimoto_distance_ecfp_model:
    #                 start_tanimoto_time = time()
    #                 tanimoto_ecfp = average_tanimoto_distance_ecfp(specific_genai_model_mol_list)
    #                 tanimoto_total_time = time() - start_tanimoto_time
    #                 print(f'calculation done in {tanimoto_total_time}s !')
    #                 library_data.extend([tanimoto_ecfp])
    #             elif is_calculate_average_tanimoto_distance_ecfp:
    #                 tanimoto_ecfp = []
    #                 library_data.extend([tanimoto_ecfp])

    #             if is_calculate_hamdiv_ecfp_model:
    #                 start_tanimoto_time = time()
    #                 tanimoto_hamdiv_ecfp = diversity_all(mols=specific_genai_model_mol_list, mode='HamDiv', hamdiv_method='ECFP', ncpu=200)
    #                 tanimoto_hamdiv_ecfp_total_time = time() - start_tanimoto_time
    #                 print(f'calculation done in {tanimoto_hamdiv_ecfp_total_time}s ! HamDiv = {tanimoto_hamdiv_ecfp}')
    #                 library_data.extend([tanimoto_hamdiv_ecfp, tanimoto_hamdiv_ecfp_total_time])
    #             elif is_calculate_hamdiv_ecfp:
    #                 tanimoto_hamdiv_ecfp, tanimoto_hamdiv_ecfp_total_time = 0, 0
    #                 library_data.extend([tanimoto_hamdiv_ecfp, tanimoto_hamdiv_ecfp_total_time])

    #             diversity_library_dict_as_list.append(library_data)


    # diversity_df = pd.DataFrame(diversity_library_dict_as_list, 
    #                             columns=diversity_library_df_column)
    
    # if is_calculate_average_tanimoto_distance_ecfp:
    #     diversity_df = diversity_df.explode('ECFP Average Tanimoto Distance').reset_index()    
    
    # # # Save the data so you don't have to rerun it hahahaha
    # diversity_df.to_csv(f'{output_suffix}.csv', index=False)
        
    diversity_df = pd.read_csv('/opt/veincent/GenAI_manuscript/scramblebench/script/diversity/SixProteinBenchmark_DiversityAnalysis_average_Tanimoto.csv')

    
    if is_calculate_average_tanimoto_distance_ecfp:
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.flatten()

        for i, protein in enumerate(protein_list):
            filtered_data = diversity_df[diversity_df['Protein'] == protein]
            sns.violinplot(data=filtered_data, x='num_sample', 
                ax=axes[i],
                y='ECFP Average Tanimoto Distance', 
                hue='Model', 
                palette= COLOR_PALLETE)
            axes[i].set_title(protein)
            axes[i].set_ylim(bottom=0)

        plt.savefig(f'{output_suffix}_TanimotoECFP.png')

    if is_calculate_hamdiv_ecfp:
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.flatten()
        for i, protein in enumerate(protein_list):
            filtered_data = diversity_df[diversity_df['Protein'] == protein]
            sns.lineplot(data=filtered_data, x='num_sample', 
                ax=axes[i],
                y='Hamiltonian Tanimoto Diversity based on ECFP', 
                hue='Model', 
                style='Model',
                markers=True,
                markersize=10,
                palette= COLOR_PALLETE)
            axes[i].set_title(protein)
            axes[i].set_ylim(bottom=0)

        plt.savefig(f'{output_suffix}_HamiltonianTanimotoDistanceECFP.png')


    if is_calculate_hamdiv_mces:
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.flatten()

        for i, protein in enumerate(protein_list):
            filtered_data = diversity_df[diversity_df['Protein'] == protein]
            sns.lineplot(data=filtered_data, x='num_sample', 
                ax=axes[i],
                y='HamDiv MCES Tanimoto Distance', 
                hue='Model', 
                style='Model',
                markers=True,
                markersize=10,
                palette= COLOR_PALLETE)
            axes[i].set_title(protein)
            axes[i].set_ylim(bottom=0)

        plt.savefig(f'{output_suffix}_HamiltonianTanimotoDistanceMCES.png')

    if is_calculate_generic_bm:
        fig, axes = plt.subplots(2, 3, figsize=(15, 10))
        axes = axes.flatten()

        for i, protein in enumerate(protein_list):
            filtered_data = diversity_df[diversity_df['Protein'] == protein]
            sns.lineplot(data=filtered_data, x='num_sample', 
                ax=axes[i],
                y='Number of Generic BM', 
                hue='Model', 
                style='Model',
                markers=True,
                markersize=10,
                palette= COLOR_PALLETE)
            axes[i].set_title(protein)
            axes[i].set_ylim(bottom=0)

        plt.savefig(f'{output_suffix}_GenericBMScaffold.png')

    

