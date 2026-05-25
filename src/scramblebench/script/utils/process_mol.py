import rdkit
from rdkit import Chem
import os
import sys
from rdkit.Chem import rdRascalMCES
from enum import Enum, IntEnum
from copy import deepcopy
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, QED, RDConfig, rdRascalMCES

sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

class PhysicoChemicalProperties(Enum):
    MW = 'MW'
    LOGP = 'logP'
    QED = 'QED'
    SA_SCORE = 'SAScore'

MOL_PROPERTY_CALCULATORS = {PhysicoChemicalProperties.MW.value: Descriptors.MolWt,
                            PhysicoChemicalProperties.LOGP.value: Crippen.MolLogP,
                            PhysicoChemicalProperties.SA_SCORE.value: sascorer.calculateScore,
                            PhysicoChemicalProperties.QED.value : QED.qed
                            }


def calculate_physicochemical_properties(mol: Chem.Mol) -> Chem.Mol:
    """Calculate the four physicochemical properties of list of mol.
    These properties are molecular weight, hydrophobicity, synthetic accessibility score and drug-likeness

    Args:
        mol (Chem.Mol): mol

    Returns:
        Chem.Mol: a ligand with physicochemical properties stored in SetProp of rdkit
    """
    for physicochemical_property, property_function in MOL_PROPERTY_CALCULATORS.items():
        mol.SetProp(physicochemical_property, str(property_function(mol)))
            
    return mol


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


def validate_mol_list(mol_list: list[Chem.Mol]) -> list[Chem.Mol]:
    # only filter valid and unique molecule
    validated_mol_l = []
    validated_smi_l = []
    for mol in mol_list:
        if not mol:
            continue

        smi = Chem.MolToSmiles(neutralize_atoms(mol))

        if smi != '' and mol.GetNumAtoms() > 1 and smi not in validated_smi_l:
            validated_mol_l.append(mol)
            validated_smi_l.append(smi)
    
    return validated_mol_l


def compute_uniqueness_percentage(mol_l) -> float:

    smi_l = [Chem.MolToSmiles(neutralize_atoms(mol)) for mol in validate_mol_list(mol_l)]
    return len(set(smi_l)) / len(list(mol_l))


def compute_validity2d_percentage(mol_l) -> float:

    return len(validate_mol_list(mol_l)) / len(list(mol_l))


def compute_generation_performance(mol_l):

    computed_mol_l = deepcopy(list(mol_l))

    performance_dict = {}
    performance_dict['total_mol'] = len(computed_mol_l)
    performance_dict['uniqueness'] = compute_uniqueness_percentage(computed_mol_l)
    performance_dict['validity2d'] = compute_validity2d_percentage(computed_mol_l)

    return performance_dict

