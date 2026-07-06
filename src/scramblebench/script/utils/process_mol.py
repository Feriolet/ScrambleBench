"""This script handles how to manipulate and fetch property of Mol objects."""
import os
import sys

from enum import Enum
from copy import deepcopy

import rdkit
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, QED, RDConfig, rdRascalMCES

sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer


class PhysicoChemicalProperties(Enum):
    """class entailing the physicochemical properties of Mol. This is usually coupled with the
    MOL_PROPERTY_CALCULATORS constant.

    Args:
        Enum (Enum): Enum class
    """
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


def neutralize_atoms(mol: Chem.Mol) -> Chem.Mol:
    """neutralize atoms. Function is copy-pasted from RDKit handbook.

    Args:
        mol (Chem.Mol): Mol object

    Returns:
        Chem.Mol: neutralized Mol
    """
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


def calculate_rms(mol1: Chem.Mol, mol2: Chem.Mol) -> float:
    """calculate the RMSD of two molecules

    Args:
        mol1 (Chem.Mol): first Mol
        mol2 (Chem.Mol): second Mol

    Returns:
        float: RMSD value
    """

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


def validate_mol_list(mol_list: list[Chem.Mol, None]) -> list[Chem.Mol]:
    """validate if the Mol in the list is actually a molecule. Empty string still can
    be converted to Chem.Mol (i.e., Chem.MolFromSmiles('') is valid), 
    so you need to count the atom numbers

    Args:
        mol_list (list[Chem.Mol]): list of Mol objects

    Returns:
        list[Chem.Mol]: list of valid Mol objects
    """
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


def compute_uniqueness_percentage(mol_l: list[Chem.Mol, None]) -> float:
    """calculate the proportion of molecules that are unique (i.e., have different SMILES)

    Args:
        mol_l (list[Chem.Mol, None]): list of Mol

    Returns:
        float: uniqueness percentage ranging from 0-1
    """

    smi_l = [Chem.MolToSmiles(neutralize_atoms(mol)) for mol in validate_mol_list(mol_l)]
    return len(set(smi_l)) / len(list(mol_l))


def compute_validity2d_percentage(mol_l: list[Chem.Mol, None]) -> float:
    """calculate the proportion of molecules that are valid (i.e., parseable by RDKit)

    Args:
        mol_l (list[Chem.Mol, None]): list of Mol

    Returns:
        float: uniqueness percentage ranging from 0-1
    """

    return len(validate_mol_list(mol_l)) / len(list(mol_l))


def compute_generation_performance(mol_l: list[Chem.Mol, None]) -> dict[str, float]:
    """calculate general model performance based on their generated ligands.

    Args:
        mol_l (list[Chem.Mol, None]): list of Mol objects

    Returns:
        dict[str, float]: dictionary detailing the total molecule, uniqueness and validity of mol
    """

    computed_mol_l = deepcopy(list(mol_l))

    performance_dict = {}
    performance_dict['total_mol'] = len(computed_mol_l)
    performance_dict['uniqueness'] = compute_uniqueness_percentage(computed_mol_l)
    performance_dict['validity2d'] = compute_validity2d_percentage(computed_mol_l)

    return performance_dict
