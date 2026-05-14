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
import argparse
from pathlib import Path
import json


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


def calculate_cpu():
    cpu_available = len(os.sched_getaffinity(0))
    CPU_BUFFER = 50
    return max(1, int(cpu_available - CPU_BUFFER))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate Diversity of generated compound")

    parser.add_argument("-i", "--input", help="sdf input", required=True, type=str)
    parser.add_argument('-o', '--output', help='json output', required=True, type=str)
    parser.add_argument("--distance", help="distance metric for diversity analysis", type=str)
    parser.add_argument("--diversity", help="diversity metric for diversity analysis", required=True, type=str)


    args = parser.parse_args()

    if not Path(args.input).is_file():
        raise ValueError(f'{args.input} does not exist')
    if Path(args.input).suffix.lower() != '.sdf':
        raise ValueError(f'{args.input} must be an sdf file')

    if Path(args.output).suffix.lower() != '.json':
        raise ValueError(f'{args.input} must be an json file')
    

    mol_l = Chem.SDMolSupplier(args.input)    

    distance_metric = args.distance
    diversity_metric = args.diversity

    if diversity_metric == 'hamdiv' and distance_metric in ['ecfp', 'mces']:
        start = time()
        result = diversity_all(mols=mol_l, mode='HamDiv', hamdiv_method=distance_metric.upper(), ncpu=calculate_cpu())
        diversity_time = time() - start
        method = f'HamDiv {distance_metric.upper()}'

    elif diversity_metric == 'average' and distance_metric == 'ecfp':
        start = time()
        result = average_tanimoto_distance_ecfp(mols=mol_l)
        diversity_time = time() - start  
        method = f'Average ECFP Tanimoto'
    else:
        start = time()
        result = diversity_all(mols=mol_l, mode=diversity_metric.upper())
        diversity_time = time() - start   
        method = diversity_metric     
    

    with open(args.output, 'w') as output_f:
        json.dump({'method': method,
                   'score': result,
                   'time': diversity_time}, output_f)

    # if is_calculate_average_tanimoto_distance_ecfp:
    #     fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    #     axes = axes.flatten()

    #     for i, protein in enumerate(protein_list):
    #         filtered_data = diversity_df[diversity_df['Protein'] == protein]
    #         sns.violinplot(data=filtered_data, x='num_sample', 
    #             ax=axes[i],
    #             y='ECFP Average Tanimoto Distance', 
    #             hue='Model', 
    #             palette= COLOR_PALLETE)
    #         axes[i].set_title(protein)
    #         axes[i].set_ylim(bottom=0)

    #     plt.savefig(f'{output_suffix}_TanimotoECFP.png')

    # if is_calculate_hamdiv_ecfp:
    #     fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    #     axes = axes.flatten()
    #     for i, protein in enumerate(protein_list):
    #         filtered_data = diversity_df[diversity_df['Protein'] == protein]
    #         sns.lineplot(data=filtered_data, x='num_sample', 
    #             ax=axes[i],
    #             y='Hamiltonian Tanimoto Diversity based on ECFP', 
    #             hue='Model', 
    #             style='Model',
    #             markers=True,
    #             markersize=10,
    #             palette= COLOR_PALLETE)
    #         axes[i].set_title(protein)
    #         axes[i].set_ylim(bottom=0)

    #     plt.savefig(f'{output_suffix}_HamiltonianTanimotoDistanceECFP.png')


    # if is_calculate_hamdiv_mces:
    #     fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    #     axes = axes.flatten()

    #     for i, protein in enumerate(protein_list):
    #         filtered_data = diversity_df[diversity_df['Protein'] == protein]
    #         sns.lineplot(data=filtered_data, x='num_sample', 
    #             ax=axes[i],
    #             y='HamDiv MCES Tanimoto Distance', 
    #             hue='Model', 
    #             style='Model',
    #             markers=True,
    #             markersize=10,
    #             palette= COLOR_PALLETE)
    #         axes[i].set_title(protein)
    #         axes[i].set_ylim(bottom=0)

    #     plt.savefig(f'{output_suffix}_HamiltonianTanimotoDistanceMCES.png')

    # if is_calculate_generic_bm:
    #     fig, axes = plt.subplots(2, 3, figsize=(15, 10))
    #     axes = axes.flatten()

    #     for i, protein in enumerate(protein_list):
    #         filtered_data = diversity_df[diversity_df['Protein'] == protein]
    #         sns.lineplot(data=filtered_data, x='num_sample', 
    #             ax=axes[i],
    #             y='Number of Generic BM', 
    #             hue='Model', 
    #             style='Model',
    #             markers=True,
    #             markersize=10,
    #             palette= COLOR_PALLETE)
    #         axes[i].set_title(protein)
    #         axes[i].set_ylim(bottom=0)

    #     plt.savefig(f'{output_suffix}_GenericBMScaffold.png')

    

