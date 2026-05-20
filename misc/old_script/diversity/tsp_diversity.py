import numpy as np
from rdkit import Chem
from rdkit.DataStructs import TanimotoSimilarity
from rdkit.Chem import rdRascalMCES, rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from rdkit.Chem.Scaffolds.MurckoScaffold import GetScaffoldForMol, MakeScaffoldGeneric
import six
import sys
sys.modules['sklearn.externals.six'] = six

from tqdm import tqdm
from multiprocessing import Pool

import networkx as nx
from python_tsp.exact import solve_tsp_dynamic_programming
from python_tsp.heuristics import solve_tsp_local_search
from tsp_utils import identify_functional_groups, GetRingSystems
# Original file based on: https://github.com/HXYfighter/HamDiv


def compute_tanimoto_distance(mol1, mol2, mces):

    total_weight_mol1 = mol1.GetNumAtoms() + mol1.GetNumBonds()  
    total_weight_mol2 = mol2.GetNumAtoms() + mol2.GetNumBonds()
    total_weight_scaf = mces.GetNumAtoms() + mces.GetNumBonds()

    return 1 - (total_weight_scaf /(total_weight_mol1 + total_weight_mol2 - total_weight_scaf))


def calculate_single_rascal_similarity(mol_tuple):
    ref_mol, mol = mol_tuple
    opts = rdRascalMCES.RascalOptions()
    opts.similarityThreshold = 0.1
    opts.completeAromaticRings = False
    opts.maxBondMatchPairs = 2000
    #does not help compute time, sometimes make it slower for the TSP algo I think
    #opts.ringMatchesRingOnly = True
    opts.equivalentAtoms = "[F,Cl,Br,I] [S,O]"
    opts.timeout = 3600

    try:
            mces_res = rdRascalMCES.FindMCES(ref_mol, mol, opts)[0]
            return Chem.MolFromSmarts(mces_res.smartsString)
    
    except IndexError:
        return Chem.MolFromSmarts('')


def calculate_single_ecfp_similarity(mol_tuple):
        ref_mol, mol = mol_tuple
        mfp_fpg = rdFingerprintGenerator.GetMorganGenerator(radius=2)

        return TanimotoSimilarity(mfp_fpg.GetFingerprint(ref_mol), mfp_fpg.GetFingerprint(mol))


def dist_array(smiles = None, mols = None, hamdiv_method = 'ECFP', ncpu = 1):
    if mols == None:
        l = len(smiles)
        mols = [Chem.MolFromSmiles(s) for s in smiles]
    else:
        l = len(mols)
    
    '''
    You can replace the Tanimoto distances of ECFPs with other molecular distance metrics!
    '''

    if hamdiv_method == 'MCES':
        comparison_list = []
        for i in range(l):
            for j in range(i + 1, l):
                comparison_list.append((mols[i], mols[j]))

        dist_matrix = np.zeros((l,l))
        with Pool(ncpu) as p:
            mces_list = p.map(calculate_single_rascal_similarity, comparison_list, chunksize=1)

        mces_index = 0
        for i in range(l):
            for j in range(i + 1, l):
                dist_matrix[i, j] = compute_tanimoto_distance(mols[i], mols[j], mces_list[mces_index])
                mces_index += 1
                dist_matrix[j, i] = dist_matrix[i, j]
            
        return dist_matrix
    
    elif hamdiv_method == 'ECFP':
        comparison_list = []
        for i in range(l):
            for j in range(i + 1, l):
                comparison_list.append((mols[i], mols[j]))

        dist_matrix = np.zeros((l,l))
        with Pool(ncpu) as p:
            ecfp_list = p.map(calculate_single_ecfp_similarity, comparison_list, chunksize=10)

        ecfp_index = 0
        for i in range(l):
            for j in range(i + 1, l):
                dist_matrix[i, j] = 1 - ecfp_list[ecfp_index]
                ecfp_index += 1
                dist_matrix[j, i] = dist_matrix[i, j]
            
        return dist_matrix



def diversity_all(smiles = None, mols = None, dists = None, mode = "HamDiv", args = None, disable = False,
                  hamdiv_method = 'ECFP', ncpu = 1):
    if mode == "Richness":
        if smiles != None:
            return len(set(smiles))
        else:
            smiles = set()
            for mol in mols:
                smiles.add(Chem.MolToSmiles(mol))
            return len(smiles)
    elif mode == "FG":
        func_groups = set()
        for i in range(len(mols) if mols is not None else len(smiles)):
            smi = Chem.MolToSmiles(mols[i]) if smiles is None else smiles[i]
            grps = identify_functional_groups(smi)
            func_groups.update(grps)

        return(len(func_groups))

    elif mode == "RS":
        ring_sys = set()
        for i in range(len(mols) if mols is not None else len(smiles)):
            smi = Chem.MolToSmiles(mols[i]) if smiles is None else smiles[i]
            grps = GetRingSystems(smi)
            ring_sys.update(grps)

        return(len(ring_sys))

    elif mode == "BM":
        scaffolds = set()
        for i in range(len(mols) if mols is not None else len(smiles)):
            if mols is not None:
                scaf = MurckoScaffold.GetScaffoldForMol(mols[i])
            else:
                mol = Chem.MolFromSmiles(smiles[i])
                scaf = MurckoScaffold.GetScaffoldForMol(mol)
            scaf_smi = Chem.MolToSmiles(scaf)
            scaffolds.update([scaf_smi])

        return(len(scaffolds))

    elif mode == "Generic_BM":
        scaffolds = set()
        for i in range(len(mols) if mols is not None else len(smiles)):
            if mols[i] is not None and mols[i].GetNumAtoms() > 0:
                smi = Chem.MolToSmiles(mols[i])
                scaf = GetScaffoldForMol(mols[i])
                scaf = MakeScaffoldGeneric(scaf)
            else:
                mol = Chem.MolFromSmiles(smiles[i])
                scaf = GetScaffoldForMol(mol)
                scaf = MakeScaffoldGeneric(scaf)
            scaf_smi = Chem.MolToSmiles(scaf, canonical=True)
            scaffolds.update([scaf_smi])

        return(len(scaffolds))

    if type(dists) is np.ndarray:
        l = len(dists)
    elif mols == None:
        l = len(smiles)
        assert l >= 2
        dists = dist_array(smiles, hamdiv_method = hamdiv_method, ncpu= ncpu)
    else:
        l = len(mols)
        assert l >= 2
        dists = dist_array(smiles, mols, hamdiv_method = hamdiv_method, ncpu = ncpu)
    
    if mode == "IntDiv":
        if l == 1:
            return 0
        return np.sum(dists) / l / (l - 1)
    elif mode == "SumDiv":
        if l == 1:
            return 0
        return np.sum(dists)/ (l - 1)
    elif mode == "Diam":
        if l == 1:
            return 0
        d_max = 0
        for i in range(l):
            for j in range(i + 1, l):
                if d_max < dists[i, j]:
                    d_max = dists[i, j]
        return d_max
    elif mode == "SumDiam":
        if l == 1:
            return 0
        sum_d_max = 0
        for i in range(l):
            d_max_i = 0
            for j in range(l):
                if j != i and d_max_i < dists[i, j]:
                    d_max_i = dists[i, j]
            sum_d_max += d_max_i
        return sum_d_max
    elif mode == "Bot":
        if l == 1:
            return 0
        d_min = 1
        for i in range(l):
            for j in range(i + 1, l):
                if d_min > dists[i, j]:
                    d_min = dists[i, j]
        return d_min
    elif mode == "SumBot":
        if l == 1:
            return 0
        sum_d_min = 0
        for i in range(l):
            d_min_i = 1
            for j in range(l):
                if j != i and d_min_i > dists[i, j]:
                    d_min_i = dists[i, j]
            sum_d_min += d_min_i
        return sum_d_min
    elif mode == "DPP":
        return np.linalg.det(1 - dists)
    elif mode.split('-')[0] == 'NCircles':
        threshold = float(mode.split('-')[1])
        circs_sum = []
        for k in tqdm(range(1), disable = disable):
            circs = np.zeros(l)
            rs = np.arange(l)
            # random.shuffle(rs)
            for i in rs:
                circs_i = 1
                for j in range(l):
                    if j != i and circs[j] == 1 and dists[i, j] <= threshold:
                        circs_i = 0
                        break
                circs[i] = circs_i
            circs_sum.append(np.sum(circs))
        return np.max(np.array(circs_sum))            
    elif mode == "HamDiv":
        total = HamDiv(dists=dists, method='local_search', hamdiv_method=hamdiv_method)
        
        return total
    
    else:
        raise TypeError('Mode Undefined!')


def HamDiv(smiles = None, mols = None, dists=None, method="greedy_tsp", hamdiv_method='ECFP', ncpu = 1):
    l = dists.shape[0] if dists is not None else len(mols) if mols is not None else len(smiles)
    if l == 1:
        return 0

    dists = dist_array(smiles, method = hamdiv_method, ncpu = ncpu) if dists is None else dists
    
    remove = np.zeros(l)
    for i in range(l):
        for j in range(i + 1, l):
            if dists[i, j] == 0:
                remove[i] = 1
    remove = np.argwhere(remove == 1)
    dists = np.delete(dists, remove, axis = 0)
    dists = np.delete(dists, remove, axis = 1)
    
    G = nx.from_numpy_array(dists)
    
    if method == "exact_dp":
        tsp, total = solve_tsp_dynamic_programming(dists)
    elif method == "christofides":
        tsp = nx.approximation.christofides(G, weight='weight')
    elif method == "greedy_tsp":
        tsp = nx.approximation.greedy_tsp(G, weight='weight')
    elif method == "simulated_annealing_tsp":
        tsp = nx.approximation.simulated_annealing_tsp(G, init_cycle="greedy", weight='weight')
    elif method == "threshold_accepting_tsp":
        tsp = nx.approximation.threshold_accepting_tsp(G, init_cycle="greedy", weight='weight')
    elif method == "local_search":
        tsp, total = solve_tsp_local_search(dists, max_processing_time=300)
    else:
        Exception("Undefined method")
    
    if method not in ["exact_dp", "local_search"]:
        total = 0
        for i in range(1, len(tsp)):
            total += dists[tsp[i - 1], tsp[i]]

    return total

if __name__ == "__main__":
    pass