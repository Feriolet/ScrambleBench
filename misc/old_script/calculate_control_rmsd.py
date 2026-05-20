import importlib
from glob import glob
from pathlib import Path
import os
import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, QED, RDConfig, rdRascalMCES, rdMolAlign, AllChem
import re


generate_summary_function = importlib.import_module('06_generate_summary_refactored')

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

    if not mol1 and not mol2:
        return None
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


vina_control_dirname = '/opt/veincent/GenAI_manuscript/output_test_multiple_numsample/Docking/vina_control_23dec'
glide_control_dirname = '/opt/veincent/GenAI_manuscript/output_test_multiple_numsample/Docking/glide_control'
input_control_dirname = '/opt/veincent/GenAI_manuscript/protein_input'
protein_name = ['5ht2c', 'drd2', 'cdk2', 'gsk3b', '3clpro', 'na']

for protein in protein_name:
    input_dir = Path(input_control_dirname)

    glob_matching_fname = f'*{generate_summary_function.create_case_insensitive_regex(protein)}*ligand.sdf'
    regex_pattern_non_alphanumeric_left_and_right = re.compile(f'(?<![A-Za-z0-9]){generate_summary_function.create_case_insensitive_regex(protein)}(?![A-Za-z0-9])')

    input_sdf = [fname for fname in input_dir.rglob(glob_matching_fname) if regex_pattern_non_alphanumeric_left_and_right.search(fname.name)][0]

    input_mol = Chem.SDMolSupplier(input_sdf)[0]
    vina_dir = Path(vina_control_dirname)
    
    glob_matching_fname = f'*{generate_summary_function.create_case_insensitive_regex(protein)}*.sdf'
    regex_pattern_non_alphanumeric_left_and_right = re.compile(f'(?<![A-Za-z0-9]){generate_summary_function.create_case_insensitive_regex(protein)}(?![A-Za-z0-9])')

    vina_sdf = [fname for fname in vina_dir.glob(glob_matching_fname) if regex_pattern_non_alphanumeric_left_and_right.search(fname.name)][0]

    vina_mol = Chem.SDMolSupplier(vina_sdf)

    print(protein)
    for mol in vina_mol:
        print(calculate_rms(input_mol, mol), mol.GetProp('docking_score'))
    

    glide_dir = Path(glide_control_dirname)

    glob_matching_fname = f'*{generate_summary_function.create_case_insensitive_regex(protein)}*.sdfgz'
    regex_pattern_non_alphanumeric_left_and_right = re.compile(f'(?<![A-Za-z0-9]){generate_summary_function.create_case_insensitive_regex(protein)}(?![A-Za-z0-9])')

    glide_sdf = [fname for fname in glide_dir.rglob(glob_matching_fname) if regex_pattern_non_alphanumeric_left_and_right.search(fname.name)][0]
    glide_mol = generate_summary_function.process_sdfgz_fname(glide_sdf)


    for mol in glide_mol:
        print(calculate_rms(input_mol, mol), mol.GetProp('r_i_docking_score'))


# [-11.019, -12.511, 0.8755, 0.8442]
# [-10.504, -10.141, 0.7024, 0.6223]
# [-9.5779, -9.8834, 0.8814, 1.2211]
# [-6.8313, -8.6500, 0.9133, 1.4014]
# [-8.1221, -6.6220, 5.2166, 3.8787]
# [-8.4364, -7.0780, 0.2527, 0.3195]