"""This file handles the input key of user's YAML config. 
Mainly used for p1_generate_config.py and other downstream path to
access the subkey field easily.
"""

import logging
import subprocess

from pathlib import Path
from collections.abc import Iterable, ItemsView
from typing import Any, Optional
from typing_extensions import Self

import numpy as np
import rdkit

from Bio.PDB import PDBParser
from rdkit import Chem
from rdkit.Chem.rdMolTransforms import ComputeCentroid
from oddt.toolkits.extras.rdkit import fixer

from scramblebench.script.config_preparation.split_protein_ligand import split_pocket_ligand
from scramblebench.script.utils.error_handler import FileDataError, FileTypeError
from scramblebench.script.config_preparation import config_constant
from scramblebench.script.utils.process_data import find_file_name_through_regex

logger = logging.getLogger(__name__)


class InputConfig:
    """class of InputConfig used for downstream analysis
    """

    def __init__(self, input_dict: dict[str, Any], key: str=config_constant.INPUT_KEY):
        """initialize class

        Args:
            input_dict (dict[str, Any]): user's prepared YAML config with input key as dictionary
            key (str, optional): the key that contains input data. Defaults to config_constant.INPUT_KEY.
        """

        self.name = key
        self.complex_name = 'complex_path'
        self.pdb_name = 'pdb_path'
        self.sdf_name = 'sdf_path'
        self.protein_name = 'name'
        self.pocket_path_name = 'pocket_path'
        self.pocket_coord_name = 'pocket_coord'

        if key:
            input_dict = input_dict[self.name]

        self.complex_value = input_dict[self.complex_name]
        self.pdb_value = input_dict[self.pdb_name]
        self.sdf_value = input_dict[self.sdf_name]
        self.protein_value = input_dict.get(self.protein_name)
        self.pocket_path_value = input_dict.get(self.pocket_path_name)
        self.pocket_coord_value = input_dict.get(self.pocket_coord_name)

    def update(self, key: str, value: str) -> Self:
        """update value of key of the class

        Args:
            key (str): valid InputConfig key
            value (str): new value that will be updated

        Raises:
            TypeError: invalid key of InputConfig

        Returns:
            Self: InputConfig
        """
        if key == self.complex_name:
            self.complex_value = value
        elif key == self.pdb_name:
            self.pdb_value = value
        elif key == self.sdf_name:
            self.sdf_value = value
        elif key == self.protein_name:
            self.protein_value = value
        else:
            raise TypeError(f'no key called {key}')

        return self


    def validate_config(self) -> None:
        """validate whether the user InputConfig is valid.
        Namely, they will check if:
        1) Complex PDB matches with Protein PDB and Ligand SDF
        2) Only 1 complex PDB, 1 protein PDB, and 1 ligand SDF.
        3) ligand SDF only has 1 ligand.
        4) file is not empty
        5) ligand is within the protein in the complex PDB.
        """

        logging.info('Validating Input Config.')
        logging.info(f'Validating Complex Structure. Complex filename: {self.complex_value}')

        check_complex_content(self.complex_value,
                            self.pdb_value,
                            self.sdf_value)


    def get_pdb_sequence(self) -> str:
        """get sequence of protein

        Returns:
            str: protein sequence
        """

        prot_seq, _ = get_residues_from_pdb(self.complex_value)

        return prot_seq


    def get_pocket_residue_info(self, chain: str) -> tuple[list[int], list[str]]:
        """fetch pocket residues of the protein

        Args:
            chain (str): PDB chain

        Returns:
            tuple[list[int], list[str]]: list of protein residues number and protein residue name
        """

        prot_seq_num = [residue.id[1] for residue in chain]
        prot_seq_name =  [residue.resname for residue in chain]

        return prot_seq_name, prot_seq_num


    def fetch_residue_is_pocket(self, cutoff: int=10) -> list[list[int, str, bool]]:
        """fetch whether each residue of protein is considered a pocket

        Args:
            cutoff (int, optional): distance cutoff of pocket residue from ligand. Defaults to 10.

        Returns:
            list[list[int, str, bool]]: a list of pocket information:
                                        (protein index, protein name, is_pocket (True/False))
        """

        d3to1 = {'CYS': 'C', 'ASP': 'D', 'SER': 'S', 'GLN': 'Q', 'LYS': 'K',
        'ILE': 'I', 'PRO': 'P', 'THR': 'T', 'PHE': 'F', 'ASN': 'N', 
        'GLY': 'G', 'HIS': 'H', 'LEU': 'L', 'ARG': 'R', 'TRP': 'W', 
        'ALA': 'A', 'VAL':'V', 'GLU': 'E', 'TYR': 'Y', 'MET': 'M'}

        if not self.pocket_path_value:
            self.pocket_path_value = split_pocket_ligand(self.complex_value, cutoff=cutoff)

        parser = PDBParser(QUIET=True)

        pocket_list = []
        for prot_model, pocket_model in zip(parser.get_structure('struct', self.pdb_value), \
                                            parser.get_structure('struct', self.pocket_path_value)):

            pocket_chain = [chain.id for chain in pocket_model]

            for chain in prot_model:
                chain_data = []

                if chain.id not in pocket_chain:
                    pocket_list.append([(residue.id[1], d3to1[residue.resname], False) \
                                        for residue in chain if residue.resname in d3to1])
                    continue

                for residue in chain:

                    pocket_name_list, pocket_num_list = self.get_pocket_residue_info(pocket_model[chain.id])

                    res_num = residue.id[1]

                    if residue.resname in d3to1:
                        if res_num in pocket_num_list and pocket_name_list[pocket_num_list.index(res_num)] == residue.resname:
                            chain_data.append((res_num, d3to1[residue.resname], True))

                        else:
                            chain_data.append((res_num, d3to1[residue.resname], False))

                pocket_list.append(chain_data)

        return pocket_list


    def get_sdf_mol(self) -> Chem.Mol:
        """get the ligand as Mol object

        Returns:
            Chem.Mol: reference ligand object
        """
        try:
            mol_l = []
            supplier = Chem.SDMolSupplier(self.sdf_value, removeHs=False)
            for mol in supplier:
                if isinstance(mol, Chem.Mol):
                    if mol != '' and mol.GetNumAtoms() > 0:
                        mol_l.append(mol)
        except OSError as exc:
            raise FileDataError(f'The content of {self.sdf_value} is likely to be empty') from exc

        if len(mol_l) < 1:
            raise FileDataError(f'There is a no ligand in {self.sdf_value}')
        if len(mol_l) > 1:
            raise FileDataError(f'There is more than 1 ligand in {self.sdf_value}')

        return mol_l[0]


    def replace_input_fnames(self, prefix_dir: str) -> None:
        """function to replace input filenames in p1_generate_config.py This is done as I
        do not want to mix intermediate files involving the input folders, so I decided that
        I want to copy the filenames to the output folder instead.

        Args:
            prefix_dir (str): directory path to copy the filenames to

        Raises:
            TypeError: if prefix_dir is not a string
        """

        if not isinstance(prefix_dir, (str, Path)):
            raise TypeError(f'{prefix_dir} is not a string or pathlib, but {type(prefix_dir)}')


        prefix_dir = Path(prefix_dir) / f'{config_constant.INPUT_KEY}_{self.protein_value}' / 'input_intermediates'
        prefix_dir.mkdir(parents=True, exist_ok=True)

        subprocess.run(['cp', self.pdb_value, self.sdf_value, self.complex_value, prefix_dir], check=True)

        self.pdb_value = str(prefix_dir / Path(self.pdb_value).name)
        self.sdf_value = str(prefix_dir / Path(self.sdf_value).name)
        self.complex_value = str(prefix_dir / Path(self.complex_value).name)


    def write(self, cutoff: int = 10, prefix_dir: str = False) -> dict[str, Any]:
        """Write the standardised config format for InputConfig

        Args:
            cutoff (int, optional): distance cutoff of pocket residues from ligand. Defaults to 10.

        Raises:
            ValueError: if the pocket coordinate is not a list, np.ndarray, or string

        Returns:
            dict[str, Any]: InputConfig as dictionary
        """
        lig_mol = Chem.SDMolSupplier(self.sdf_value)[0]

        if prefix_dir:
            self.replace_input_fnames(prefix_dir=prefix_dir)

        if not self.pocket_path_value:
            self.pocket_path_value = split_pocket_ligand(self.complex_value, cutoff=cutoff)

        if not self.pocket_coord_value:
            self.pocket_coord_value = np.array(list(ComputeCentroid(lig_mol.GetConformer(0), ignoreHs=True)))

        data = {self.complex_name: str(Path(self.complex_value).resolve()),
                self.pdb_name: str(Path(self.pdb_value).resolve()),
                self.sdf_name: str(Path(self.sdf_value).resolve()),
                self.pocket_path_name: str(Path(self.pocket_path_value).resolve())}

        if isinstance(self.pocket_coord_value, str) and len(self.pocket_coord_value.split(',')) == 3:
            data[self.pocket_coord_name] = self.pocket_coord_value
        elif isinstance(self.pocket_coord_value, (list, np.ndarray)):
            data[self.pocket_coord_name] = ','.join([str(np.round(coord, 2)) for coord in self.pocket_coord_value])
        else:
            raise ValueError(f'{self.pocket_coord_value}')
        if self.protein_value:
            data[self.protein_name] = self.protein_value

        # because write method is only used for p1, I don't need to add {'input' : data} like others
        return data


class InputNonDirpathConfig:
    """Input Config Class specifically for p1_generate_config.py argument
    """

    def __init__(self, config_data: dict[str, Any]):
        """initiate class

        Args:
            config_data (dict[str, Any]): user config YAML input as dictionary

        Raises:
            KeyError: if 'dirpath_input' not in the key
        """

        self.name = config_constant.INPUT_KEY

        if self.name not in config_data:
            raise KeyError(f'You forgot to put "{self.name}" key in your config file!')

        input_data = config_data[self.name]

        self.inputstructure_dict = {input_keyname: InputConfig(structinput_data, key=None) \
                                    for input_keyname, structinput_data in input_data.items()}


    def __iter__(self) -> Iterable[ItemsView[str, InputConfig]]:
        """custom iteration method to return the list of InputConfig as dictionary

        Returns:
            Iterable[ItemsView[str, InputConfig]]: list of user's (protein name, InputConfig)
        """
        return iter(list(self.inputstructure_dict.items()))

    def update(self, input_keyname: str, input_key: str, value: Any) -> Self:
        """update value of key of the class

        Args:
            input_keyname (str): protein name
            input_key (str): valid key name of InputConfig
            value (Any): new value that will be updated

        Raises:
            KeyError: if input_key is not valid

        Returns:
            Self: InputNOnDirpathConfig class
        """
        try:
            self.inputstructure_dict[input_keyname] = self.inputstructure_dict[input_keyname].update(input_key, value)
        except KeyError as exc:
            raise KeyError(f'no key called {input_keyname}') from exc

        return self


    # see __iter__ method for how this is possible.
    def validate_config(self) -> None:
        """validate each InputConfig
        """
        for input_key in self:
            (input_key[-1].validate_config())


    def get_input_list(self) -> list[str]:
        """get all protein names in the keys

        Returns:
            list[str]: list of protein names
        """
        return [inputstruct.input_name_value for inputstruct in self.inputstructure_dict.values()]


    def write(self, cutoff:int=10, prefix_dir=None) -> dict[str, Any]:
        """Return standardised dictionary without changing the input structure format

        Args:
            cutoff (int, optional): cutoff of pocket residues from ligand. Defaults to 10.

        Returns:
            dict[str, Any]: standardised dictionary for user's input config
        """

        input_data = {}
        for key, inputstruct in self.inputstructure_dict.items():
            input_data[key] = inputstruct.update(key=inputstruct.protein_name, value=key)\
                                         .write(cutoff=cutoff, prefix_dir=prefix_dir)

        return {self.name: input_data}


    def write_inputstruct(self, cutoff:int=10) -> dict[str, Any]:
        """return standardised dictionary for a single InputConfig

        Args:
            cutoff (int, optional): cutoff of pocket residues from ligand. Defaults to 10.

        Returns:
            dict[str, Any]: standardised dictionary for a single input config
        """

        assert len(list(self.inputstructure_dict.keys())) == 1

        for key, inputstruct in self.inputstructure_dict.items():
            input_data = inputstruct.write(cutoff=cutoff)
            input_data[inputstruct.protein_name] = key

            return {self.name: input_data}


class InputDirpathConfig:
    """Input Config Class specifically for p1_generate_config.py using --dirpath_input argument
    """

    def __init__(self, config_data: dict[str, Any]):
        """initiate class

        Args:
            config_data (dict[str, Any]): user config YAML input as dictionary
        """

        self.name = config_constant.INPUT_DIR_KEY
        input_dir_data = config_data[self.name]
        self.dirpath_name = 'dirpath'
        self.dirpath_value = input_dir_data[self.dirpath_name]

        self.complex_name = 'complex_path'
        self.pdb_name = 'pdb_path'
        self.sdf_name = 'sdf_path'
        self.input_data = InputConfig(self.search_for_filepath())


    def search_for_filepath(self) -> dict[str, dict]:
        """use glob to search for the user's protein and ligand files. 
        There MUST be ONLY ONE file present for each type.

        Returns:
            dict[str, dict]: minimally required config data for InputConfig
        """

        pathdir_dirname = Path(self.dirpath_value)
        assert pathdir_dirname.resolve().is_dir()

        protein_dirname_list = pathdir_dirname.iterdir()
        protein_dict = {}
        for protein_dirname in protein_dirname_list:
            assert Path(protein_dirname).is_dir()
            protein_name = Path(protein_dirname).stem
            protein_dict[protein_name] = {}

            protein_pdb = find_file_name_through_regex(character='protein',
                                                       file_format='.pdb',
                                                       dirname=protein_dirname)
            if len(protein_pdb) == 0:
                continue

            ligand_sdf = find_file_name_through_regex(character='ligand',
                                                       file_format='.sdf',
                                                       dirname=protein_dirname)
            if len(ligand_sdf) == 0:
                continue

            complex_pdb = find_file_name_through_regex(character='complex',
                                                       file_format='.pdb',
                                                       dirname=protein_dirname)
            if len(complex_pdb) == 0:
                continue

            protein_dict[protein_name][self.pdb_name] = protein_pdb
            protein_dict[protein_name][self.sdf_name] = ligand_sdf
            protein_dict[protein_name][self.complex_name] = complex_pdb

        return {config_constant.INPUT_KEY: protein_dict}


    def validate_config(self) -> None:
        """validate each InputConfig
        """
        self.input_data.validate_config()


    def write(self, cutoff: int = 10) -> dict[str, Any]:
        """Return standardised dictionary for list of InputConfig

        Args:
            cutoff (int, optional): cutoff of pocket residues from ligand. Defaults to 10.

        Returns:
            dict[str, Any]: standardised dictionary for user's input config
        """
        return self.input_data.write(cutoff=cutoff)


def create_case_insensitive_regex(pattern: str) -> str:
    """create a template where regex allows case insensitive search

    Args:
        pattern (str): queried characters

    Returns:
        str: regex formatted queried characters as case insensitive
    """
    return f"{''.join([ '[' + char.upper() + char.lower() + ']' for char in pattern])}"


def check_pdb_path(pdb_fname: str) -> None:
    """check whether protein PDB file is valid

    Args:
        pdb_fname (str): protein PDB filename

    Raises:
        TypeError: if filename is not a string
        FileNotFoundError: if filename does not exist
        FileTypeError: if filename is not a PDB format
    """
    if not isinstance(pdb_fname, str):
        raise TypeError(f'{pdb_fname} is not a string!')

    pdb_path = Path(pdb_fname)
    if not pdb_path.is_file():
        raise FileNotFoundError(f'No file is found for {pdb_fname}')
    if not pdb_path.name.endswith('.pdb'):
        raise FileTypeError(f'{pdb_fname} does not end with PDB for protein')

def check_complex_path(complex_fname: str) -> None:
    """check whether complex PDB file is valid

    Args:
        complex_fname (str): complex PDB filename

    Raises:
        TypeError: if filename is not a string
        FileNotFoundError: if filename does not exist
        FileTypeError: if filename is not a PDB format
    """

    if not isinstance(complex_fname, str):
        raise TypeError(f'{complex_fname} is not a string!')
    pdb_path = Path(complex_fname)
    if not pdb_path.is_file():
        raise FileNotFoundError(f'No file is found for {complex_fname}')
    if not pdb_path.name.endswith('.pdb'):
        raise FileTypeError(f'{complex_fname} does not end with PDB for complex')


def check_sdf_path(sdf_fname: str) -> None:
    """check whether ligand SDF file is valid

    Args:
        sdf_fname (str): ligand SDF filename

    Raises:
        TypeError: if filename is not a string
        FileNotFoundError: if filename does not exist
        FileTypeError: if filename is not an SDF format
    """
    if not isinstance(sdf_fname, str):
        raise TypeError(f'{sdf_fname} is not a string!')
    sdf_path = Path(sdf_fname)
    if not sdf_path.is_file():
        raise FileNotFoundError(f'No file is found for {sdf_fname}')
    if not sdf_path.name.endswith('.sdf'):
        raise FileTypeError(f'{sdf_fname} does not end with SDF for ligand')


def get_residues_from_pdb(pdb_file_path: str) -> Optional[tuple[str, list[Optional[str]]]]:
    """get protein residue information from PDB

    Args:
        pdb_file_path (str): PDB filename

    Returns:
        Optional[tuple[str, list[Optional[str]]]]: returns a protein sequence (str) 
        and ligand residue name (if available)
    """

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
    """check if PDB file has protein residues labelled by ATOM

    Args:
        pdb_fname (str): PDB filename

    Returns:
        bool: True/False
    """

    with open(pdb_fname, 'r', encoding='utf-8') as pdb_fn:
        pdb_data = pdb_fn.readlines()
        for pdb_line in pdb_data:
            if 'ATOM' in pdb_line:
                return True

    return False


def validate_complex_and_protein_content(complex_fname: str, protein_fname: str) -> None:
    """check the content of protein PDB and complex PDB

    Args:
        complex_fname (str): complex PDB filename
        protein_fname (str): protein PDB filename

    Raises:
        FileDataError: if protein PDB does not have residue
        FileDataError: if complex PDB does not have residue
        FileDataError: if there is no ligand in complex PDB
        FileDataError: if there are multiple ligands in complex PDB
        FileDataError: if there is a ligand in protein PDB
        FileDataError: if the protein sequence of protein and complex PDB does not match
    """

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


def parse_sdf(ligand_fname: str) -> list[Optional[Chem.Mol]]:
    """read SDF filename to return a Mol object

    Args:
        ligand_fname (str): ligand SDF filename

    Raises:
        FileDataError: if the SDF file is empty

    Returns:
        list[Optional[Chem.Mol]]: ligand reference mol object
    """
    mol_l = []
    try:
        supplier = Chem.SDMolSupplier(ligand_fname, removeHs=False)
        for mol in supplier:
            if isinstance(mol, Chem.Mol):
                if mol != '' and mol.GetNumAtoms() > 0:
                    mol_l.append(mol)
    except OSError as exc:
        raise FileDataError(f'The content of {ligand_fname} is likely to be empty') from exc
    return mol_l


def compute_ligand_centroid_distance(ligand1: Chem.Mol, ligand2: Chem.Mol) -> np.float_:
    """compute the centroid distance of the ligand coordinate

    Args:
        ligand1 (Chem.Mol): reference ligand
        ligand2 (Chem.Mol): queried ligand

    Returns:
        np.float_: distance between the first and second ligand's centroid
    """
    ligand_centroid = np.array(list(ComputeCentroid(ligand1.GetConformer(0), ignoreHs=True)))
    complex_centroid = np.array(list(ComputeCentroid(ligand2.GetConformer(0), ignoreHs=True)))

    return np.linalg.norm(ligand_centroid - complex_centroid)


def validate_complex_and_ligand_content(complex_fname: str, ligand_fname: str) -> None:
    """check the content of complex PDB and ligand SDF

    Args:
        complex_fname (str): complex PDB filename
        ligand_fname (str): ligand SDF filename

    Raises:
        FileDataError: if pocket residue is not sufficient
        FileDataError: if ligand from complex PDB and ligand SDF is not the same
        FileDataError: if the position of ligand from complex PDB and ligand SDF is different
    """

    POCKET_RESIDUE_NUMBER_THRESHOLD = 5
    POCKET_LIGAND_DISTANCE_CUTOFF = 10
    TRANSLATED_LIGAND_DISTANCE_CUTOFF = 1

    check_sdf_path(ligand_fname)
    check_complex_path(complex_fname)

    ligand_mol_l = parse_sdf(ligand_fname)


    complex_ = Chem.MolFromPDBFile(complex_fname, sanitize=False)
    extracted_pocket, extracted_ligand = fixer.ExtractPocketAndLigand(complex_, cutoff=POCKET_LIGAND_DISTANCE_CUTOFF)

    pocket_residue_number = len({atom.GetPDBResidueInfo().GetResidueNumber() \
                                     for atom in extracted_pocket.GetAtoms()})
    if pocket_residue_number < POCKET_RESIDUE_NUMBER_THRESHOLD:
        raise FileDataError(f'We only detected {pocket_residue_number} pocket residue in {complex_fname}')

    ligand_lig = Chem.RemoveAllHs(ligand_mol_l[0])
    extracted_ligand = Chem.RemoveAllHs(extracted_ligand)

    # quick match, too lazy to neutralise the ligand.
    if Chem.MolToSmiles(ligand_lig, canonical=True) != Chem.MolToSmiles(extracted_ligand, canonical=True):
        raise FileDataError(f'We detected different ligand from {ligand_fname} and {complex_fname}. \
                            If the molecule is the same, kindly check its protonation/isomer state.')

    if compute_ligand_centroid_distance(ligand_lig, extracted_ligand) > TRANSLATED_LIGAND_DISTANCE_CUTOFF:
        raise FileDataError(f'The ligand from {ligand_fname} and {complex_fname} \
                            are likely to be translated by more than {TRANSLATED_LIGAND_DISTANCE_CUTOFF} A')


def check_complex_content(complex_fname: str, protein_fname: str, ligand_fname: str) -> None:
    """check the content of complex PDB with ligand SDF and protein PDB

    Args:
        complex_fname (str): complex PDB filename
        protein_fname (str): protein PDB filename
        ligand_fname (str): ligand SDF filename
    """

    validate_complex_and_protein_content(complex_fname=complex_fname, protein_fname=protein_fname)
    validate_complex_and_ligand_content(complex_fname=complex_fname, ligand_fname=ligand_fname)
