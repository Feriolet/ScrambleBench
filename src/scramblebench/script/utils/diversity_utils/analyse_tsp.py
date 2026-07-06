"""This util file handles the main interface for p4_analyse_diversity.py"""
import json
import os
import sys
import argparse
import logging

from time import time
from pathlib import Path

import numpy as np
import rdkit

from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from tsp_diversity import diversity_all


logger = logging.getLogger(__name__)

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

            tanimoto_dist.append(1 - DataStructs.TanimotoSimilarity(mol_fp_i, mol_fp_j))

        tanimoto_dist_array = np.array(tanimoto_dist)
        average_tanimoto_dist.append(np.average(tanimoto_dist_array))

    return average_tanimoto_dist


def calculate_cpu():
    """claculate number of available cpu

    Returns:
        int: available cpu
    """
    if hasattr(os, 'sched_getaffinity'):
        cpu_available = len(os.sched_getaffinity(0))
    else:
        cpu_available = os.cpu_count() or 1
    CPU_BUFFER = 50
    return max(1, int(cpu_available - CPU_BUFFER))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Calculate Diversity of generated compound")

    parser.add_argument("-i", "--input", help="sdf input", required=True, type=str)
    parser.add_argument('-o', '--output', help='json output', required=True, type=str)
    parser.add_argument("--distance", help="distance metric for diversity analysis", type=str)
    parser.add_argument("--diversity", help="diversity metric for diversity analysis", required=True, type=str)

    logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='%(asctime)s - %(module)s: - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
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
        result = average_tanimoto_distance_ecfp(mol_list=mol_l)
        diversity_time = time() - start
        method = 'Average ECFP Tanimoto'
    else:
        start = time()
        result = diversity_all(mols=mol_l, mode=diversity_metric.upper())
        diversity_time = time() - start
        method = diversity_metric


    with open(args.output, 'w', encoding='utf-8') as output_f:
        json.dump({'method': method,
                   'score': result,
                   'time': diversity_time}, output_f)
