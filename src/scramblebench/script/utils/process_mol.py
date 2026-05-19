import rdkit
from rdkit import Chem
from rdkit.Chem import rdRascalMCES

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
