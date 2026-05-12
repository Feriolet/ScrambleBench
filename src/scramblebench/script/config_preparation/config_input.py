
from typing import Any, Optional
from pathlib import Path

from scramblebench.script.utils.error_handler import FileDataError, FileTypeError
from scramblebench.script.config_preparation import config_constant
import Bio
import rdkit
from Bio.PDB import PDBParser
from rdkit import Chem
import os
from scramblebench.script.utils.split_protein_ligand import split_pocket_ligand
from oddt.toolkits.extras.rdkit import fixer
import numpy as np
import logging
import re

logger = logging.getLogger(__name__)


class InputConfig:

    def __init__(self, config_data: dict[str, Any]):

        self.name = config_constant.INPUT_KEY
        input_data = config_data[self.name]

        self.inputstructure_dict = {input_keyname: InputStructure(structinput_data) for input_keyname, structinput_data in input_data.items()}

    def __iter__(self):
        """Returns the iterator object itself."""
        return iter(list(self.inputstructure_dict.items()))
    
    def update(self, input_keyname, input_key, value):
        try:
            self.inputstructure_dict[input_keyname] = self.inputstructure_dict[input_keyname].update(input_key, value)
        except KeyError:
            raise KeyError(f'no key called {input_keyname}')

        return self
    
    # see __iter__ method for how this is possible.
    def validate_config(self) -> None:
        for input in self:
            (input[-1].validate_config()) 

    def get_input_list(self):
        return [inputstruct.input_name_value for inputstruct in self.inputstructure_dict.values()]

    def write(self, cutoff:int=10) -> dict[str, Any]:
        input_data = {}
        for key, inputstruct in self.inputstructure_dict.items():
            input_data[key] = inputstruct.write(cutoff=cutoff)
            
        return {self.name: input_data}
    
class InputStructure:

    def __init__(self, input_dict: dict[str, Any]):

        self.complex_name = 'complex_path'
        self.pdb_name = 'pdb_path'
        self.sdf_name = 'sdf_path'
        self.protein_name = 'name'
        self.pocket_path_name = 'pocket_path'
        self.pocket_coord_name = 'pocket_coord'

        self.complex_value = input_dict[self.complex_name]
        self.pdb_value = input_dict[self.pdb_name]
        self.sdf_value = input_dict[self.sdf_name]
        self.protein_value = input_dict.get(self.protein_name)
        self.pocket_path_value = input_dict.get(self.pocket_path_name)
        self.pocket_coord_value = input_dict.get(self.pocket_coord_name)

    def update(self, key: str, value: str):
        if key == self.complex_name:
            self.complex_value = value
        elif key == self.pdb_name:
            self.pdb_value = value
        elif key == self.sdf_name:
            self.sdf_value = value
        else:
            raise TypeError(f'no key called {key}')

        return self

    def validate_config(self) -> None:
        logging.info('Validating Input Config.')
        logging.info(f'Validating Complex Structure. Complex filename: {self.complex_value}')

        check_complex_content(self.complex_value,
                            self.pdb_value,
                            self.sdf_value)
    
    def write(self, cutoff: int = 10) -> dict[str, Any]:
        lig_mol = Chem.SDMolSupplier(self.sdf_value)[0]

        if not self.pocket_path_value:
            self.pocket_path_value = split_pocket_ligand(self.complex_value, cutoff=cutoff)
        
        if not self.pocket_coord_value:
            self.pocket_coord_value = np.array(list(Chem.rdMolTransforms.ComputeCentroid(lig_mol.GetConformer(0), ignoreHs=True)))

        data = {self.complex_name: str(Path(self.complex_value).resolve()),
                self.pdb_name: str(Path(self.pdb_value).resolve()),
                self.sdf_name: str(Path(self.sdf_value).resolve()),
                self.pocket_path_name: str(Path(self.pocket_path_value).resolve()),
                self.pocket_coord_name: f" {','.join([str(np.round(coord, 2)) for coord in self.pocket_coord_value])}"} 
        
        if self.protein_value:
            data[self.protein_name] = self.protein_value

        return data


class InputDirConfig:

    def __init__(self, config_data: dict[str, Any]):

        self.name = config_constant.INPUT_DIR_KEY
        input_dir_data = config_data[self.name]
        self.dirpath_name = 'dirpath'
        self.dirpath_value = input_dir_data[self.dirpath_name]

        self.complex_name = 'complex_path'
        self.pdb_name = 'pdb_path'
        self.sdf_name = 'sdf_path'
        self.input_data = InputConfig(self.search_for_filepath())

        
    def search_for_filepath(self):
        
        pathdir_dirname = Path(self.dirpath_value)
        assert pathdir_dirname.resolve().is_dir()

        protein_dirname_list = pathdir_dirname.iterdir()
        protein_dict = {}
        for protein_dirname in protein_dirname_list:
            assert Path(protein_dirname).is_dir()

            protein_name = Path(protein_dirname).stem
            
            protein_dict[protein_name] = {}

            ### Find potential Protein PDB
            glob_matching_fname = f'*{create_case_insensitive_regex("protein")}*.pdb'
            regex_pattern_non_alphanumeric_left_and_right = re.compile(f'(?<![A-Za-z0-9]){create_case_insensitive_regex("protein")}(?![A-Za-z0-9])')

            matched_fname_list = [fname for fname in protein_dirname.glob(glob_matching_fname) if regex_pattern_non_alphanumeric_left_and_right.search(fname.name)]
            if len(matched_fname_list) == 0:
                logging.warning(f'No pdb file found for {protein_name}. Skipping...')
                continue
            elif len(matched_fname_list) > 1:
                logging.warning(f'Multiple pdb file found for {protein_name}. Skipping...')
                continue
            else:
                protein_dict[protein_name][self.pdb_name] = str(matched_fname_list[0])

            ### Find potential Ligand SDF
            glob_matching_fname = f'*{create_case_insensitive_regex("ligand")}*.sdf'
            regex_pattern_non_alphanumeric_left_and_right = re.compile(f'(?<![A-Za-z0-9]){create_case_insensitive_regex("ligand")}(?![A-Za-z0-9])')

            matched_fname_list = [fname for fname in protein_dirname.glob(glob_matching_fname) if regex_pattern_non_alphanumeric_left_and_right.search(fname.name)]
            if len(matched_fname_list) == 0:
                logging.warning(f'No sdf file found for {protein_name}. Skipping...')
                continue
            elif len(matched_fname_list) > 1:
                logging.warning(f'Multiple sdf file found for {protein_name}. Skipping...')
                continue
            else:
                protein_dict[protein_name][self.sdf_name] = str(matched_fname_list[0])

            ### Find potential complex PDB
            glob_matching_fname = f'*{create_case_insensitive_regex("complex")}*.pdb'
            regex_pattern_non_alphanumeric_left_and_right = re.compile(f'(?<![A-Za-z0-9]){create_case_insensitive_regex("complex")}(?![A-Za-z0-9])')

            matched_fname_list = [fname for fname in protein_dirname.glob(glob_matching_fname) if regex_pattern_non_alphanumeric_left_and_right.search(fname.name)]
            if len(matched_fname_list) == 0:
                logging.warning(f'No pdb file found for {protein_name}. Skipping...')
                continue
            elif len(matched_fname_list) > 1:
                logging.warning(f'Multiple pdb file found for {protein_name}. Skipping...')
                continue
            else:
                protein_dict[protein_name][self.complex_name] = str(matched_fname_list[0])           

        return {config_constant.INPUT_KEY: protein_dict}
    
    def validate_config(self):
        self.input_data.validate_config()
    
    def write(self, cutoff: int = 10) -> dict[str, Any]:
        return self.input_data.write(cutoff=cutoff)
    
def create_case_insensitive_regex(pattern: str) -> str:
    return f"{''.join([ '[' + char.upper() + char.lower() + ']' for char in pattern])}"


def check_pdb_path(pdb_fname: str) -> None:
    if not isinstance(pdb_fname, str):
        raise TypeError(f'{pdb_fname} is not a string!')

    pdb_path = Path(pdb_fname)
    if not pdb_path.is_file():
        raise FileNotFoundError(f'No file is found for {pdb_fname}')
    if not pdb_path.name.endswith('.pdb'):
        raise FileTypeError(f'{pdb_fname} does not end with PDB for protein')

def check_complex_path(complex_fname: str) -> None:

    if not isinstance(complex_fname, str):
        raise TypeError(f'{complex_fname} is not a string!')
    pdb_path = Path(complex_fname)
    if not pdb_path.is_file():
        raise FileNotFoundError(f'No file is found for {complex_fname}')
    if not pdb_path.name.endswith('.pdb'):
        raise FileTypeError(f'{complex_fname} does not end with PDB for complex')
    
def check_sdf_path(sdf_fname: str) -> None:
    if not isinstance(sdf_fname, str):
        raise TypeError(f'{sdf_fname} is not a string!')
    sdf_path = Path(sdf_fname)
    if not sdf_path.is_file():
        raise FileNotFoundError(f'No file is found for {sdf_fname}')
    if not sdf_path.name.endswith('.sdf'):
        raise FileTypeError(f'{sdf_fname} does not end with SDF for ligand')

def get_residues_from_pdb(pdb_file_path: str) -> Optional[tuple[str, list[Optional[str]]]]:

    # You can use a dict to convert three letter code to one letter code
    d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
    'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
    'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
    'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

    parser = PDBParser(QUIET=True)
    structure = parser.get_structure('struct', pdb_file_path)    

    # iterate each model, chain, and residue
    # printing out the sequence for each chain

    for model in structure:
        for chain in model:
            prot_seq = ''.join([d3to1[residue.resname] for residue in chain if residue.resname in d3to1])
            lig_seq = [residue.resname for residue in chain if residue.resname not in d3to1]

            return prot_seq, lig_seq

def check_pdb_contains_residue(pdb_fname: str) -> bool:
    
    with open(pdb_fname, 'r') as pdb_fn:
        pdb_data = pdb_fn.readlines()
        for pdb_line in pdb_data:
            if 'ATOM' in pdb_line:
                return True
    
    return False

def validate_complex_and_protein_content(complex_fname: str, protein_fname: str) -> None:
    
    check_complex_path(complex_fname)
    check_pdb_path(protein_fname)

    if not check_pdb_contains_residue(complex_fname):
        raise FileDataError(f'There is no protein in complex filename {complex_fname}')
    
    if not check_pdb_contains_residue(protein_fname):
        raise FileDataError(f'There is no protein in protein filename {protein_fname}')
    
    complex_prot_seq, complex_lig_seq = get_residues_from_pdb(complex_fname)
    if len(complex_lig_seq) < 1:
        raise FileDataError(f'There is no ligand in complex filename {complex_fname}')
    if len(complex_lig_seq) > 1:
        raise FileDataError(f'There is more than 1 ligand in complex filename {complex_fname}')

    protein_prot_seq, protein_lig_seq = get_residues_from_pdb(protein_fname)
    if len(protein_lig_seq) > 0:
        raise FileDataError(f'There is a ligand in protein filename {protein_fname}')
    
    if complex_prot_seq != protein_prot_seq:
        raise FileDataError(f'There is a mismatch in protein residue in {protein_fname} and {complex_fname}')


def parse_sdf(ligand_fname: str) -> list[Optional[Chem.Mol]] :

    mol_l = []
    try:
        supplier = Chem.SDMolSupplier(ligand_fname, removeHs=False)
        for mol in supplier:
            if isinstance(mol, Chem.Mol):
                if mol != '' and mol.GetNumAtoms() > 0:
                    mol_l.append(mol)
    except OSError:
        raise FileDataError(f'The content of {ligand_fname} is likely to be empty')
    return mol_l

def compute_ligand_centroid_distance(ligand1: Chem.Mol, ligand2: Chem.Mol) -> np.float_:
    ligand_centroid = np.array(list(Chem.rdMolTransforms.ComputeCentroid(ligand1.GetConformer(0), ignoreHs=True)))
    complex_centroid = np.array(list(Chem.rdMolTransforms.ComputeCentroid(ligand2.GetConformer(0), ignoreHs=True)))
    
    return np.linalg.norm(ligand_centroid - complex_centroid)

def validate_complex_and_ligand_content(complex_fname: str, ligand_fname: str) -> None:

    POCKET_RESIDUE_NUMBER_THRESHOLD = 5
    POCKET_LIGAND_DISTANCE_CUTOFF = 10
    TRANSLATED_LIGAND_DISTANCE_CUTOFF = 1

    check_sdf_path(ligand_fname)
    check_complex_path(complex_fname)

    ligand_mol_l = parse_sdf(ligand_fname)
    if len(ligand_mol_l) < 1:
        raise FileDataError(f'There is a no ligand in {ligand_fname}')
    if len(ligand_mol_l) > 1:
        raise FileDataError(f'There is more than 1 ligand in {ligand_fname}')

    complex_ = Chem.MolFromPDBFile(complex_fname, sanitize=False)
    extracted_pocket, extracted_ligand = fixer.ExtractPocketAndLigand(complex_, cutoff=POCKET_LIGAND_DISTANCE_CUTOFF)

    pocket_residue_number = len(set([atom.GetPDBResidueInfo().GetResidueNumber() for atom in extracted_pocket.GetAtoms()]))
    if pocket_residue_number < POCKET_RESIDUE_NUMBER_THRESHOLD:
        raise FileDataError(f'We only detected {pocket_residue_number} residue surrounding the ligand in {complex_fname}')

    ligand_lig = Chem.RemoveAllHs(ligand_mol_l[0])
    extracted_ligand = Chem.RemoveAllHs(extracted_ligand)

    # quick match, too lazy to neutralise the ligand.   
    if Chem.MolToSmiles(ligand_lig, canonical=True) != Chem.MolToSmiles(extracted_ligand, canonical=True):
        raise FileDataError(f'We detected different ligand from {ligand_fname} and {complex_fname}. If the molecule is the same, kindly check its protonation/isomer state.')

    if compute_ligand_centroid_distance(ligand_lig, extracted_ligand) > TRANSLATED_LIGAND_DISTANCE_CUTOFF:
        raise FileDataError(f'The ligand from {ligand_fname} and {complex_fname} are likely to be translated by more than {TRANSLATED_LIGAND_DISTANCE_CUTOFF} A')

def check_complex_content(complex_fname: str, protein_fname: str, ligand_fname: str) -> None:

    validate_complex_and_protein_content(complex_fname=complex_fname, protein_fname=protein_fname)
    validate_complex_and_ligand_content(complex_fname=complex_fname, ligand_fname=ligand_fname)


