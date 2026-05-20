import numpy as np
import pandas as pd

import gzip
import json
import os
import re
import sqlite3
import sys


from collections import defaultdict
from enum import Enum, IntEnum
from dataclasses import dataclass
from itertools import chain
from pathlib import Path
from rdkit import Chem
from rdkit.Chem import Descriptors, Crippen, QED, RDConfig, rdRascalMCES, rdMolAlign, AllChem
from typing import Optional, Iterator
from meeko import (MoleculePreparation,
                   PDBQTMolecule,
                   PDBQTWriterLegacy,
                   RDKitMolCreate)

sys.path.append(os.path.join(RDConfig.RDContribDir, 'SA_Score'))
import sascorer

# require
# 1. Raw sdf folder (cheminformatics_input)
# 2. Vina and Glide folder from Genbench
# 3. Validity3D from GenBench
# 4. Phase score
# 5. Vina and Glide folder redocked

# WARNING: Because this analysis require the different isomer to be differentiated,
# glide is required for easeness in integration.
# This is because the ligand is prepared through glide LigPrep, so they have the E/Z isomer
# running easydock remove this isomer when converted using RDKit. So, make sure you have a tag
# to differentiate each isomer if you wish to differentiate them.
# I may try to make a flexible workflow to include the summary if easydock or glide is not run
# in the future, but it is what it is for now.

# Also, for 3CLpro, the <original name> property is wrong for half of the control, because
# one of the control does not have a smiles, so their name is off by one.


MOL_PROPERTY_CALCULATORS = {'MW': Descriptors.MolWt,
                            'logP': Crippen.MolLogP,
                            'SAScore': sascorer.calculateScore,
                            'QED' : QED.qed
                            }


class GenBenchDockingMethod(Enum):
    VINA_MININPLACE = "Vina score"
    VINA_INPLACE = "Minimized Vina score"
    GLIDE_MININPLACE = "Glide score"
    GLIDE_INPLACE = "Minimized Glide score"


class LigandSource(Enum):
    GENBENCH_FOLDER = "genbench_folder"
    VINA_REDOCK_FOLDER = "vina_redock_folder"
    GLIDE_REDOCK_FOLDER = "glide_redock_folder"
    GLIDE_PHASE_FOLDER = 'glide_phase_folder'
    RAW_INPUT_FOLDER = "raw_input_folder"


class FileCount(IntEnum):
    EMPTY = 0
    ONE_FILE = 1


class GenBenchForceFieldMinimisation(Enum):
    UNMINIMISED = 'unminimised'
    MINIMISED = 'minimised'


class GenBenchDockingMinimisation(Enum):
    UNMINIMISED = 'unminimised'
    MINIMISED = 'minimised'


class PhysicoChemicalProperties(Enum):
    MW = 'MW'
    LOGP = 'logP'
    QED = 'QED'
    SA_SCORE = 'SAScore'

 
@dataclass(kw_only=True)
class Raw3DInputLigandRecord:
    ligand_id: str
    smi: str
    mol_block: str
    molecular_weight: float
    logP: float
    SA_score: float
    QED: float
    protein_name: str
    genai_model: str
    num_sample: int


@dataclass(kw_only=True)
class EasydockLigandRecord:
    ligand_id: str
    stereo_id: int
    smi: str
    docked_mol_block: str
    docking_score: float
    protein_name: str


@dataclass(kw_only=True)
class GlideRedockedLigandRecord:
    ligand_id: str
    stereo_id: int
    smi: str
    docked_mol_block: str
    docking_score: float
    protein_name: str
    genai_model: str
    num_sample: int


@dataclass(kw_only=True)
class GlidePhaseLigandRecord:
    ligand_id: str
    stereo_id: int
    smi: str
    matched_ligand_site: str
    phase_score: float
    protein_name: str
    genai_model: str
    num_sample: int


@dataclass(kw_only=True)
class GenbenchGlideLigandRecord:
    ligand_id: str
    docked_mol_block: str
    docking_score: float
    docking_method: GenBenchDockingMethod
    forcefield_minimisation: GenBenchForceFieldMinimisation
    protein_name: str
    genai_model: str
    num_sample: int


@dataclass(kw_only=True)
class GenbenchVinaLigandRecord:
    ligand_id: str
    docked_pdbqt_mol_block: str
    docking_score: float
    docking_method: GenBenchDockingMethod
    forcefield_minimisation: GenBenchDockingMinimisation
    protein_name: str
    genai_model: str
    num_sample: int


#Merged_database = dict[ligand_ID, dict[Folder, FolderData]]
Merged_schema = dict[tuple[str, int, str], dict[LigandSource, dict[str, object]]]


def create_case_insensitive_regex(pattern: str) -> str:
    return f"{''.join([ '[' + char.upper() + char.lower() + ']' for char in pattern])}"


def process_sdf_fname(sdf_fname: Path) -> list[Chem.Mol | None]:
    mol_l :list[Chem.Mol | None] = []
    for mol in Chem.SDMolSupplier(sdf_fname, removeHs=False):
        if not isinstance(mol, Chem.Mol):
            mol_l.append(None)
        elif isinstance(mol, Chem.Mol) and mol.GetNumAtoms() <= 0:
            mol_l.append(None)
        else:
            mol_l.append(mol)
        
    return mol_l


def process_sdfgz_fname(sdfgz_fname: str | Path) -> list[Chem.Mol | None]:
    mol_l :list[Chem.Mol | None] = []
    with gzip.open(sdfgz_fname, 'rb') as f:
        for mol in Chem.ForwardSDMolSupplier(f, removeHs=False):
            if not isinstance(mol, Chem.Mol):
                mol_l.append(None)
            elif isinstance(mol, Chem.Mol) and mol.GetNumAtoms() <= 0:
                mol_l.append(None)
            else:
                mol_l.append(mol)
            
        return mol_l
    

def process_easydock_vina_folder(easydock_dir: str, protein_name: str) -> Iterator[tuple[str, int, EasydockLigandRecord]]:
    easydock_fetch_column_name = ['id', 'stereo_id', 'smi', 'mol_block', 'docking_score']
 
    easydock_dirpath = Path(easydock_dir)
    if not easydock_dirpath.exists():
        raise ValueError('Easydock path is missing')
    if not easydock_dirpath.is_dir():
        raise TypeError(f'Easydock Path {easydock_dirpath} is not a directory')
    
    glob_matching_fname = f'*{create_case_insensitive_regex(protein_name)}*.db'
    regex_pattern_non_alphanumeric_left_and_right = re.compile(f'(?<![A-Za-z0-9]){create_case_insensitive_regex(protein_name)}(?![A-Za-z0-9])')

    matched_fname_list = [fname for fname in easydock_dirpath.glob(glob_matching_fname) if regex_pattern_non_alphanumeric_left_and_right.search(fname.name)]
    matched_fname_count = 0
    for matched_fname_count, db_fname in enumerate(matched_fname_list, start=1):
        with sqlite3.connect(db_fname, timeout=90) as conn:
            conn.row_factory = sqlite3.Row
            cur = conn.cursor()
            data = cur.execute(f"""SELECT {",".join(easydock_fetch_column_name)} FROM mols""")

            for ligand_row in data.fetchall():
                yield ligand_row['id'], int(ligand_row['stereo_id']), EasydockLigandRecord(ligand_id=ligand_row['id'],
                                                             stereo_id=int(ligand_row['stereo_id']),
                                                             smi=ligand_row['smi'],
                                                             docked_mol_block=ligand_row['mol_block'],
                                                             docking_score=ligand_row['docking_score'],
                                                             protein_name=protein_name)

    if matched_fname_count == FileCount.EMPTY:
        print(f'No db is found for {protein_name}')
    elif matched_fname_count > FileCount.ONE_FILE:
        print(f'We detected multiple db files for {protein_name}. Please terminate if you do not intend to do this')


def process_glide_redocked_folder(glide_dir: str, protein_name: str) -> Iterator[Optional[tuple[str, int, GlideRedockedLigandRecord]]]:

    glide_dirpath = Path(glide_dir)
    if not glide_dirpath.exists():
        raise ValueError('Glide path is missing')
    if not glide_dirpath.is_dir():
        raise TypeError(f'Glide Path {glide_dirpath} is not a directory')
    
    glob_matching_fname = f'*{create_case_insensitive_regex(protein_name)}*SP*.sdf*'
    regex_pattern_non_alphanumeric_left_and_right = re.compile(f'(?<![A-Za-z0-9]){create_case_insensitive_regex(protein_name)}(?![A-Za-z0-9])')

    matched_fname_list = [fname for fname in glide_dirpath.glob(glob_matching_fname) if regex_pattern_non_alphanumeric_left_and_right.search(fname.name)]
    matched_fname_count = 0
    for matched_fname_count, glide_fname in enumerate(matched_fname_list, start=1):
        if glide_fname.name.endswith('sdfgz'):
            mol_l = process_sdfgz_fname(glide_fname)
        elif glide_fname.name.endswith('sdf'):
            mol_l = process_sdf_fname(glide_fname)

        for mol in mol_l:
            if isinstance(mol, Chem.Mol):
                #s_lp_Variant format is id-stereoid with stereo_id starting from 1
                ligand_id = str(mol.GetProp('_Name'))
                ligand_stereo_id = int(mol.GetProp('s_lp_Variant').split('-')[-1]) -1
                if not mol.HasProp('num_sample'):
                    mol.SetProp('num_sample', '0')

                yield ligand_id, ligand_stereo_id, GlideRedockedLigandRecord(ligand_id=ligand_id,
                                                                             stereo_id=ligand_stereo_id,
                                                                             smi=Chem.MolToSmiles(mol, canonical=True),
                                                                             docked_mol_block=Chem.MolToMolBlock(mol),
                                                                             docking_score=mol.GetProp('r_i_docking_score'),
                                                                             protein_name=protein_name,
                                                                             genai_model=mol.GetProp('GenAI_Model'),
                                                                             num_sample=int(mol.GetProp('num_sample')))
    if matched_fname_count == FileCount.EMPTY:
        print(f'No sdf/sdfgz is found for {protein_name} for glide analysis')
    elif matched_fname_count > FileCount.ONE_FILE:
        print(f'We detected multiple sdf/sdfgz files for {protein_name} for glide analysis. Please terminate if you do not intend to do this')


def process_glide_phase_folder(glide_phase_dir: str, protein_name: str) -> Iterator[Optional[tuple[str, int, GlidePhaseLigandRecord]]]:

    glide_phase_dirpath = Path(glide_phase_dir)
    if not glide_phase_dirpath.exists():
        raise ValueError('Glide phase path is missing')
    if not glide_phase_dirpath.is_dir():
        raise TypeError(f'Glide phase path {glide_phase_dirpath} is not a directory')
    
    glob_matching_fname = f'*{create_case_insensitive_regex(protein_name)}*.sdf*'
    regex_pattern_non_alphanumeric_left_and_right = re.compile(f'(?<![A-Za-z0-9]){create_case_insensitive_regex(protein_name)}(?![A-Za-z0-9])')

    matched_fname_list = [fname for fname in glide_phase_dirpath.glob(glob_matching_fname) if regex_pattern_non_alphanumeric_left_and_right.search(fname.name)]
    matched_fname_count = 0
    for matched_fname_count, glide_fname in enumerate(matched_fname_list, start=1):
        if glide_fname.name.endswith('sdfgz'):
            mol_l = process_sdfgz_fname(glide_fname)
        elif glide_fname.name.endswith('sdf'):
            mol_l = process_sdf_fname(glide_fname)

        for mol in mol_l:
            if isinstance(mol, Chem.Mol):
                #s_lp_Variant format is id-stereoid with stereo_id starting from 1
                ligand_id = str(mol.GetProp('_Name'))
                ligand_stereo_id = int(mol.GetProp('s_lp_Variant').split('-')[-1]) -1
                if not mol.HasProp('num_sample'):
                    mol.SetProp('num_sample', '0')

                yield ligand_id, ligand_stereo_id, GlidePhaseLigandRecord(ligand_id=ligand_id,
                                                                             stereo_id=ligand_stereo_id,
                                                                             smi=Chem.MolToSmiles(mol, canonical=True),
                                                                             matched_ligand_site=mol.GetProp('s_phase_Matched_Ligand_Sites'),
                                                                             phase_score=mol.GetProp('r_phase_PhaseScreenScore'),
                                                                             protein_name=protein_name,
                                                                             genai_model=mol.GetProp('GenAI_Model'),
                                                                             num_sample=int(mol.GetProp('num_sample')))
    if matched_fname_count == FileCount.EMPTY:
        print(f'No sdf/sdfgz is found for {protein_name} for glide phase analysis')
    elif matched_fname_count > FileCount.ONE_FILE:
        print(f'We detected multiple sdf/sdfgz files for {protein_name} for glide phase analysis. Please terminate if you do not intend to do this')


def process_genbench_vina_score(forcefield_minimisation_dir, docking_list, num_sample, docking_method):
    docking_list_array = np.array(docking_list)
    docking_list_without_nan = docking_list_array[~np.isnan(docking_list_array)]
    if docking_method == GenBenchDockingMethod.VINA_MININPLACE:
        docking_minimisation = GenBenchDockingMinimisation.MINIMISED.value
    elif docking_method == GenBenchDockingMethod.VINA_INPLACE:
        docking_minimisation = GenBenchDockingMinimisation.UNMINIMISED.value

    for num_id in range(num_sample):
        for pdbqt_fname in forcefield_minimisation_dir.glob(f'*num{num_sample}_{num_id}_{docking_minimisation}*.pdbqt'):
            if pdbqt_fname:
                ligand_id = pdbqt_fname.name.split(f'_{docking_minimisation}')[0]
                docking_score = float(docking_list_without_nan[0])
                docking_list_without_nan = np.delete(docking_list_without_nan, [0])
                pdbqt_block = str(pdbqt_fname.read_text())
                
                yield ligand_id, {
                    'ligand_id': ligand_id,
                    'docking_score': docking_score,
                    'docked_mol_block': pdbqt_block,
                    'docking_method': docking_method
                }
    

def process_genbench_vina_folder(genbench_dir: Path, protein_name: str):
    vina_dir = genbench_dir.joinpath('Vina')
    json_dirpath = genbench_dir.joinpath('json_output')
    num_sample = int(genbench_dir.parent.name.split('_')[-1])

    for vina_dir_content in vina_dir.iterdir():
        genai_model_name = vina_dir_content.stem
        for forcefield_minimisation_dir in vina_dir_content.iterdir():
            if '_minimised' in forcefield_minimisation_dir.name.lower():
                forcefield_minimisation_status = GenBenchForceFieldMinimisation.MINIMISED
            elif '_unminimised' in forcefield_minimisation_dir.name.lower():
                forcefield_minimisation_status = GenBenchForceFieldMinimisation.UNMINIMISED
            else:
                raise ValueError(f'unknown forcefield minimisation in {protein_name}. Dirpath: {forcefield_minimisation_dir.name}')
        
            glob_matching_fname =f'*{create_case_insensitive_regex(protein_name)}*{forcefield_minimisation_status.value}*.json'

            regex_pattern_non_alphanumeric_left_and_right = re.compile(f'(?<![A-Za-z0-9]){create_case_insensitive_regex(genai_model_name)}(?![A-Za-z0-9])(?=.*(?<![A-Za-z0-9]){forcefield_minimisation_status.value}(?![A-Za-z0-9]))')

            matched_fname_list = [fname for fname in json_dirpath.glob(glob_matching_fname) if regex_pattern_non_alphanumeric_left_and_right.search(fname.name)]
            matched_json_fname_count = 0
            for matched_json_fname_count, json_fname in enumerate(matched_fname_list, start=1):
                json_data = json.load(open(json_fname))

                vina_data = defaultdict(list)
                if GenBenchDockingMethod.VINA_INPLACE.value in json_data:
                    vina_data[GenBenchDockingMethod.VINA_INPLACE].append(json_data[GenBenchDockingMethod.VINA_INPLACE.value])
                if GenBenchDockingMethod.VINA_MININPLACE.value in json_data:
                    vina_data[GenBenchDockingMethod.VINA_MININPLACE].append(json_data[GenBenchDockingMethod.VINA_MININPLACE.value])
                
                for docking_method, docking_list in vina_data.items():
                    for ligand_id, vina_data in process_genbench_vina_score(forcefield_minimisation_dir, docking_list, num_sample, docking_method):
                        yield ligand_id, protein_name, GenbenchVinaLigandRecord(
                            ligand_id=ligand_id,
                            docked_pdbqt_mol_block=vina_data['docked_mol_block'],
                            docking_score=vina_data['docking_score'],
                            forcefield_minimisation=forcefield_minimisation_status,
                            docking_method=vina_data['docking_method'],
                            protein_name=protein_name,
                            genai_model=genai_model_name,
                            num_sample=num_sample
                        )
            if matched_json_fname_count == FileCount.EMPTY:
                print(f'We did not find any json file for {protein_name}, FF minimisation {forcefield_minimisation_status.value} for {genai_model_name} model, numsample {num_sample}')
            if matched_json_fname_count > FileCount.ONE_FILE:
                print(f'We found more than 1 json file for {protein_name}, FF minimisation {forcefield_minimisation_status.value} for {genai_model_name} model, numsample {num_sample}')


def process_genbench_glide_score(forcefield_minimisation_dir, docking_method):

    if docking_method == GenBenchDockingMethod.GLIDE_MININPLACE:
        docking_minimisation = GenBenchDockingMinimisation.MINIMISED.value
    elif docking_method == GenBenchDockingMethod.GLIDE_INPLACE:
        docking_minimisation = GenBenchDockingMinimisation.UNMINIMISED.value
    else:
        raise ValueError(f'no docking method supported for {docking_method} for Glide Genbench analysis')

    glob_matching_fname =f'glide_scoring*{docking_minimisation}*.sdf*'

    regex_pattern_non_alphanumeric_left_and_right = re.compile(f'(?<![A-Za-z0-9]){docking_minimisation}(?![A-Za-z0-9])')

    matched_fname_list = [fname for fname in forcefield_minimisation_dir.glob(glob_matching_fname) if regex_pattern_non_alphanumeric_left_and_right.search(fname.name)]
    matched_fname_count = 0
    for matched_fname_count, glide_fname in enumerate(matched_fname_list, start=1):
        if glide_fname.name.endswith('sdfgz'):
            mol_l = process_sdfgz_fname(glide_fname)
        elif glide_fname.name.endswith('sdf'):
            mol_l = process_sdf_fname(glide_fname)
        else:
            raise TypeError(f'please insert sdfgz or sdf file. We detected {glide_fname} instead')

        for mol in mol_l[1:]:
            if isinstance(mol, Chem.Mol):
                #s_lp_Variant format is id-stereoid with stereo_id starting from 1
                ligand_id = str(mol.GetProp('_Name'))
                yield ligand_id, {'ligand_id': ligand_id,
                    'docking_score': mol.GetProp('r_i_docking_score'),
                    'docked_mol_block': Chem.MolToMolBlock(mol),
                    'docking_method': docking_method
                }

    if matched_fname_count == FileCount.EMPTY:
        print(f'No sdf/sdfgz is found for {protein_name} for glide genbench3d analysis')
    elif matched_fname_count > FileCount.ONE_FILE:
        print(f'We detected multiple sdf/sdfgz files for {protein_name} for glide genbench3d analysis. Please terminate if you do not intend to do this')


def process_genbench_glide_folder(genbench_dir: Path, protein_name: str):
    glide_dir = genbench_dir.joinpath('Glide')
    json_dirpath = genbench_dir.joinpath('json_output')
    num_sample = int(genbench_dir.parent.name.split('_')[-1])

    for glide_dir_content in glide_dir.iterdir():
        genai_model_name = glide_dir_content.stem
        for forcefield_minimisation_dir in glide_dir_content.iterdir():
            if '_minimised' in forcefield_minimisation_dir.name.lower():
                forcefield_minimisation_status = GenBenchForceFieldMinimisation.MINIMISED
            elif '_unminimised' in forcefield_minimisation_dir.name.lower():
                forcefield_minimisation_status = GenBenchForceFieldMinimisation.UNMINIMISED
            else:
                raise ValueError(f'unknown forcefield minimisation in {protein_name}. Dirpath: {forcefield_minimisation_dir.name}')
        
            glob_matching_fname =f'*{create_case_insensitive_regex(protein_name)}*{forcefield_minimisation_status.value}*.json'

            regex_pattern_non_alphanumeric_left_and_right = re.compile(f'(?<![A-Za-z0-9]){create_case_insensitive_regex(genai_model_name)}(?![A-Za-z0-9])(?=.*(?<![A-Za-z0-9]){forcefield_minimisation_status.value}(?![A-Za-z0-9]))')

            matched_fname_list = [fname for fname in json_dirpath.glob(glob_matching_fname) if regex_pattern_non_alphanumeric_left_and_right.search(fname.name)]
            matched_json_fname_count = 0
            for matched_json_fname_count, json_fname in enumerate(matched_fname_list, start=1):
                json_data = json.load(open(json_fname))

                glide_docking_method_list = []
                if GenBenchDockingMethod.GLIDE_INPLACE.value in json_data:
                    glide_docking_method_list.append(GenBenchDockingMethod.GLIDE_INPLACE)
                if GenBenchDockingMethod.GLIDE_MININPLACE.value in json_data:
                    glide_docking_method_list.append(GenBenchDockingMethod.GLIDE_MININPLACE)
                
                for docking_method in glide_docking_method_list:
                    for ligand_id, glide_data in process_genbench_glide_score(forcefield_minimisation_dir, docking_method):
                        yield ligand_id, protein_name, GenbenchGlideLigandRecord(
                            ligand_id=ligand_id,
                            docked_mol_block=glide_data['docked_mol_block'],
                            docking_score=glide_data['docking_score'],
                            forcefield_minimisation=forcefield_minimisation_status,
                            docking_method=glide_data['docking_method'],
                            protein_name=protein_name,
                            genai_model=genai_model_name,
                            num_sample=num_sample
                        )
            if matched_json_fname_count == FileCount.EMPTY:
                print(f'We did not find any json file for {protein_name}, FF minimisation {forcefield_minimisation_status.value} for {genai_model_name} model, numsample {num_sample}')
            if matched_json_fname_count > FileCount.ONE_FILE:
                print(f'We found more than 1 json file for {protein_name}, FF minimisation {forcefield_minimisation_status.value} for {genai_model_name} model, numsample {num_sample}')


def process_genbench_folder(root_dir: str, protein_name: str, genbench_dirname: str):

    root_dirpath: Path = Path(root_dir)
    if not root_dirpath.exists():
        raise ValueError('Genbench root path is missing')
    if not root_dirpath.is_dir():
        raise TypeError(f'Genbench root path {root_dirpath} is not a directory')

    glob_matching_dirname = f'*{create_case_insensitive_regex(protein_name)}*'
    regex_pattern_non_alphanumeric_left_and_right = re.compile(f'(?<![A-Za-z0-9]){create_case_insensitive_regex(protein_name)}(?![A-Za-z0-9])')

    matched_dirname_list = [fname for fname in root_dirpath.glob(glob_matching_dirname) if regex_pattern_non_alphanumeric_left_and_right.search(fname.name)]
    docking_dict: dict[tuple[str, str], dict[str, object]] = defaultdict(dict)
    for protein_dirpath in matched_dirname_list:
        if protein_dirpath.is_dir():
            genbench_folder = protein_dirpath.joinpath(genbench_dirname)
            if genbench_folder.is_dir():
                for ligand, protein, docking_data in chain(process_genbench_glide_folder(genbench_folder, protein_name), process_genbench_vina_folder(genbench_folder, protein_name)):
                    docking_keyword = ''
                    if docking_data.forcefield_minimisation == GenBenchForceFieldMinimisation.MINIMISED:
                        docking_keyword += 'FF_minimised_'
                    else:
                        docking_keyword += 'FF_unminimised_'
                    docking_keyword += docking_data.docking_method.name.lower()

                    if 'vina' in docking_data.docking_method.name.lower():
                        docking_dict[(ligand, protein)].update({'ligand_id': docking_data.ligand_id,
                                                        f'{docking_keyword}_docked_pdbqt_mol_block': docking_data.docked_pdbqt_mol_block,
                                                        f'{docking_keyword}_docking_score': docking_data.docking_score,
                                                        'protein_name': protein,
                                                        'genai_model': docking_data.genai_model,
                                                        'num_sample': docking_data.num_sample})
                    elif 'glide' in docking_data.docking_method.name.lower():
                        docking_dict[(ligand, protein)].update({'ligand_id': docking_data.ligand_id,
                                                        f'{docking_keyword}_docked_mol_block': docking_data.docked_mol_block,
                                                        f'{docking_keyword}_docking_score': docking_data.docking_score,
                                                        'protein_name': protein,
                                                        'genai_model': docking_data.genai_model,
                                                        'num_sample': docking_data.num_sample})                        
                
    for key, value in docking_dict.items():
        ligand_id, protein_name = key
        yield ligand_id, protein_name, value


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


def process_raw_3d_sdf_folder(root_dir: str, protein_name: str, raw_3d_sdf_dirname: str):
    root_dirpath: Path = Path(root_dir)
    if not root_dirpath.exists():
        raise ValueError('Raw 3D SDF root path is missing')
    if not root_dirpath.is_dir():
        raise TypeError(f'Raw 3D SDF root path {root_dirpath} is not a directory')

    glob_matching_dirname = f'*{create_case_insensitive_regex(protein_name)}*'
    regex_pattern_non_alphanumeric_left_and_right = re.compile(f'(?<![A-Za-z0-9]){create_case_insensitive_regex(protein_name)}(?![A-Za-z0-9])')

    matched_dirname_list = [fname for fname in root_dirpath.glob(glob_matching_dirname) if regex_pattern_non_alphanumeric_left_and_right.search(fname.name)]
    for protein_dirpath in matched_dirname_list:
        if protein_dirpath.is_dir():
            raw_3d_sdf_folder = protein_dirpath.joinpath(raw_3d_sdf_dirname)
            if raw_3d_sdf_folder.is_dir():
                glob_matching_fname = f'*{create_case_insensitive_regex(protein_name)}*.sdf*'
                regex_pattern_non_alphanumeric_left_and_right = re.compile(f'(?<![A-Za-z0-9]){create_case_insensitive_regex(protein_name)}(?![A-Za-z0-9])')

                matched_fname_list = [fname for fname in raw_3d_sdf_folder.glob(glob_matching_fname) if regex_pattern_non_alphanumeric_left_and_right.search(fname.name)]
                matched_fname_count = 0
                for matched_fname_count, raw_input_fname in enumerate(matched_fname_list, start=1):
                    if raw_input_fname.name.endswith('sdfgz'):
                        mol_l = process_sdfgz_fname(raw_input_fname)
                    elif raw_input_fname.name.endswith('sdf'):
                        mol_l = process_sdf_fname(raw_input_fname)

                    for mol in mol_l:
                        if isinstance(mol, Chem.Mol):
                            #s_lp_Variant format is id-stereoid with stereo_id starting from 1
                            ligand_id = str(mol.GetProp('_Name'))
                            mol = calculate_physicochemical_properties(mol)

                            yield ligand_id, Raw3DInputLigandRecord(ligand_id=ligand_id,
                                                                    smi=Chem.MolToSmiles(mol, canonical=True),
                                                                    mol_block=Chem.MolToMolBlock(mol),
                                                                    molecular_weight=mol.GetProp(PhysicoChemicalProperties.MW.value),
                                                                    logP=mol.GetProp(PhysicoChemicalProperties.LOGP.value),
                                                                    SA_score=mol.GetProp(PhysicoChemicalProperties.SA_SCORE.value),
                                                                    QED=mol.GetProp(PhysicoChemicalProperties.QED.value),
                                                                    protein_name=protein_name,
                                                                    genai_model=mol.GetProp('GenAI_Model'),
                                                                    num_sample=int(mol.GetProp('num_sample')))
                if matched_fname_count == FileCount.EMPTY:
                    print(f'No sdf/sdfgz is found for {protein_name} for input raw 3d analysis')
    

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


def calculate_rms(mol_block_l1, mol_block_l2):

    rms_list = []
    for mb1, mb2 in zip(mol_block_l1, mol_block_l2):
        try:
            mol1  = Chem.MolFromMolBlock(mb1, removeHs=False)
            mol2 = Chem.MolFromMolBlock(mb2, removeHs=False)
        except TypeError:
            rms_list.append(None)
            continue
        if not mol1 and not mol2:
            rms_list.append(None)
            continue
        mol1 = neutralize_atoms(Chem.RemoveHs(mol1))
        mol2 = neutralize_atoms(Chem.RemoveHs(mol2))

        try:
            rms_list.append(Chem.rdMolAlign.CalcRMS(mol1, mol2))
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
            rms_list.append(Chem.rdMolAlign.CalcRMS(mol1, mol2, map=[matches]))

    return rms_list


# Taken from easydock repository
def assign_bonds_from_template(template_mol, mol):
    # explicit hydrogends are removed from carbon atoms (chiral hydrogens) to match pdbqt mol,
    # e.g. [NH3+][C@H](C)C(=O)[O-]
    template_mol_ = Chem.Mol(template_mol)
    template_mol_ = Chem.AddHs(template_mol_, explicitOnly=True,
                               onlyOnAtoms=[a.GetIdx() for a in template_mol_.GetAtoms() if
                                            a.GetAtomicNum() != 6])
    mol = AllChem.AssignBondOrdersFromTemplate(template_mol_, mol)
    Chem.SanitizeMol(mol)
    Chem.AssignStereochemistry(mol, cleanIt=True, force=True, flagPossibleStereoCenters=True)
    return mol


# Taken from easydock repository
def pdbqt2molblock(pdbqt_block, template_molblock):
    """
    The function takes PDBQT block with one or more poses and converts top pose to MDL MOL format.
    :param pdbqt_block: a single string with a single PDBQT block (a single pose)
    :param template_molblock: MolBlock of a reference structure to assign bond orders
    :param boron_replacement: indicate whether to try to return boron atoms instead af carbon ones
    :return: a single string with a MOL block, if conversion failed returns None
    """

    try:
        template_mol  = Chem.MolFromMolBlock(template_molblock, removeHs=False)
    except TypeError:
        return None

    if not isinstance(template_mol, Chem.Mol):
        return None
    
    mol_block = None
    mol_id = template_mol.GetProp('_Name')
    try:
        pdbqt_mol = PDBQTMolecule(pdbqt_block, is_dlg=False, skip_typing=True, poses_to_read=1)
        rdkitmol_list = RDKitMolCreate.from_pdbqt_mol(pdbqt_mol)
        rdkit_mol = rdkitmol_list[0]
    except Exception:
        print(f"{mol_id}, conversion PDB to MOL was failed.\n")
        return None

    try:
        mol = assign_bonds_from_template(template_mol, rdkit_mol)

        mol.SetProp("_Name", mol_id)
        mol_block = Chem.MolToMolBlock(mol)
    except Exception:
        print(f"{mol_id}, could not assign bond orders while parsing PDB. Trying to fix.\n")

    return mol_block


if __name__ == '__main__':
    # this is your main directory containing all of the raw results, even the one
    # outside the automated script, such as the easydock vina docking and glide docking

    root_pathdir = '/opt/veincent/GenAI_manuscript/output_test_multiple_numsample'
    raw_3d_sdf_foldername = 'cheminformatics_input_prepared'

    redock_pathdir = '/opt/veincent/GenAI_manuscript/output_test_multiple_numsample/Docking'
    chem42_gsk_redock_pathdir = '/opt/veincent/GenAI_manuscript/output_test_multiple_numsample/Docking/gsk3b_chem42'
    redock_vina_foldername = 'vina_output'
    redock_glide_foldername = 'glide_input'
    phase_screen_foldername = 'phase_score'

    protein_name_l = ['5ht2c', 'drd2', 'cdk2', 'gsk3b', '3clpro', 'na']

    genbench_foldername = 'cheminformatics_analysis'
    num_sample_list = [100, 200, 500]
    is_minimised = True
    is_unminimised = True
    GenAI_model_l = ['Pocket2Mol', 'DiffSBDD', 'PocketFlow', 'Chem42', 'PMDM', 'Lingo3DMol']
    
    
    merged_database: Merged_schema = defaultdict(dict)
    for protein_name in protein_name_l:
        for ligand_id, ligand_stereoid, easydock_data in process_easydock_vina_folder(os.path.join(redock_pathdir, redock_vina_foldername), protein_name):
            merged_database[(ligand_id, ligand_stereoid, protein_name)][LigandSource.VINA_REDOCK_FOLDER] = {"ligand_id": easydock_data.ligand_id,
                                                                            "stereo_id": easydock_data.stereo_id,
                                                                            "redocked_vina_mol_block": easydock_data.docked_mol_block,
                                                                            "redocked_vina_smi": easydock_data.smi,
                                                                            "protein_name": easydock_data.protein_name,
                                                                            "redocked_vina_docking_score": easydock_data.docking_score}
        # special case because I use chain B of GSK3b instead of chain A for Chem42, so the rmsd is acting weird
        if protein_name == 'gsk3b':
            for ligand_id, ligand_stereoid, easydock_data in process_easydock_vina_folder(os.path.join(chem42_gsk_redock_pathdir, redock_vina_foldername), protein_name):
                merged_database[(ligand_id, ligand_stereoid, protein_name)][LigandSource.VINA_REDOCK_FOLDER] = {"ligand_id": easydock_data.ligand_id,
                                                                                "stereo_id": easydock_data.stereo_id,
                                                                                "redocked_vina_mol_block": easydock_data.docked_mol_block,
                                                                                "redocked_vina_smi": easydock_data.smi,
                                                                                "protein_name": easydock_data.protein_name,
                                                                                "redocked_vina_docking_score": easydock_data.docking_score}
        for ligand_id, input_data in process_raw_3d_sdf_folder(os.path.join(root_pathdir), protein_name, raw_3d_sdf_foldername):
            for stereo_id in range(100):
                if (ligand_id, stereo_id, protein_name) in merged_database or stereo_id == 0:
                    merged_database[(ligand_id, stereo_id, protein_name)][LigandSource.RAW_INPUT_FOLDER] = {"ligand_id": input_data.ligand_id,
                                                                                                                 "stereo_id": stereo_id,
                                                                                                                 "input_smi": input_data.smi,
                                                                                                                 "input_mol_block": input_data.mol_block,
                                                                                                                 "MW": input_data.molecular_weight,
                                                                                                                 "log P": input_data.logP,
                                                                                                                 "SA_Score": input_data.SA_score,
                                                                                                                 "QED": input_data.QED,
                                                                                                                 "protein_name": input_data.protein_name,
                                                                                                                 "genai_model": input_data.genai_model,
                                                                                                                 "num_sample": input_data.num_sample}
                else:
                    break            

        #print(f' vina {count_nested_dicts(merged_database)=}')
        for ligand_id, ligand_stereoid, glide_redocked_data in process_glide_redocked_folder(os.path.join(redock_pathdir, redock_glide_foldername), protein_name):
            merged_database[(ligand_id, ligand_stereoid, protein_name)][LigandSource.GLIDE_REDOCK_FOLDER] = {"ligand_id": glide_redocked_data.ligand_id,
                                                                            "stereo_id": glide_redocked_data.stereo_id,
                                                                            "redocked_glide_mol_block": glide_redocked_data.docked_mol_block,
                                                                            "redocked_glide_smi": glide_redocked_data.smi,
                                                                            "protein_name": glide_redocked_data.protein_name,
                                                                            "redocked_glide_docking_score": glide_redocked_data.docking_score,
                                                                            "genai_model": glide_redocked_data.genai_model,
                                                                            "num_sample": glide_redocked_data.num_sample}
        
        if protein_name == 'gsk3b':
            for ligand_id, ligand_stereoid, glide_redocked_data in process_glide_redocked_folder(os.path.join(chem42_gsk_redock_pathdir, redock_glide_foldername), protein_name):
                merged_database[(ligand_id, ligand_stereoid, protein_name)][LigandSource.GLIDE_REDOCK_FOLDER] = {"ligand_id": glide_redocked_data.ligand_id,
                                                                                "stereo_id": glide_redocked_data.stereo_id,
                                                                                "redocked_glide_mol_block": glide_redocked_data.docked_mol_block,
                                                                                "redocked_glide_smi": glide_redocked_data.smi,
                                                                                "protein_name": glide_redocked_data.protein_name,
                                                                                "redocked_glide_docking_score": glide_redocked_data.docking_score,
                                                                                "genai_model": glide_redocked_data.genai_model,
                                                                                "num_sample": glide_redocked_data.num_sample}
            
        #print(f'glide {count_nested_dicts(merged_database)=}')
        for ligand_id, ligand_stereoid, glide_phase_data in process_glide_phase_folder(os.path.join(redock_pathdir, phase_screen_foldername), protein_name):
            merged_database[(ligand_id, ligand_stereoid, protein_name)][LigandSource.GLIDE_PHASE_FOLDER] = {"ligand_id": glide_phase_data.ligand_id,
                                                                            "stereo_id": glide_phase_data.stereo_id,
                                                                            "phase_matched_site": glide_phase_data.matched_ligand_site,
                                                                            "phase_glide_smi": glide_phase_data.smi,
                                                                            "protein_name": glide_phase_data.protein_name,
                                                                            "glide_phase_score": glide_phase_data.phase_score,
                                                                            "genai_model": glide_phase_data.genai_model,
                                                                            "num_sample": glide_phase_data.num_sample}


        for ligand_id, protein_name, genbench_data in process_genbench_folder(os.path.join(root_pathdir), protein_name, genbench_foldername):
            # check stereo_id up until 100
            for stereo_id in range(100):
                if (ligand_id, stereo_id, protein_name) in merged_database or stereo_id == 0:
                    merged_database[(ligand_id, stereo_id, protein_name)][LigandSource.GENBENCH_FOLDER] = genbench_data | {'stereo_id': stereo_id}
                else:
                    break

    rows = []

    for sample_id, sources in merged_database.items():
        row = {}
        for source_data in sources.values():
            row.update(source_data)
        rows.append(row)

    df = pd.DataFrame(rows)

    df['FF_minimised_vina_mininplace_docked_mol_block'] = [pdbqt2molblock(pdbqtblock, glidemolblock) for pdbqtblock, glidemolblock in zip(df['FF_minimised_vina_mininplace_docked_pdbqt_mol_block'], 
                                                                                                                                          df['FF_minimised_glide_mininplace_docked_mol_block'])]

    df['glide_redock-input_rms'] = calculate_rms(df['redocked_glide_mol_block'], df['input_mol_block'])
    df['glide_redock-FF_glide_mininplace_rms'] = calculate_rms(df['FF_minimised_glide_mininplace_docked_mol_block'], df['redocked_glide_mol_block'])
    df['FF_glide_mininplace_rms-input_rms'] = calculate_rms(df['FF_minimised_glide_mininplace_docked_mol_block'], df['input_mol_block'])

    df['vina_redock-input_rms'] = calculate_rms(df['redocked_vina_mol_block'], df['input_mol_block'])
    df['vina_redock-FF_vina_mininplace_rms'] = calculate_rms(df['FF_minimised_vina_mininplace_docked_mol_block'], df['redocked_vina_mol_block'])
    df['FF_vina_mininplace_rms-input_rms'] = calculate_rms(df['FF_minimised_vina_mininplace_docked_mol_block'], df['input_mol_block'])


    df.to_parquet('save_all_feature.parquet')
    df.drop(columns=['redocked_vina_mol_block',
                     'redocked_glide_mol_block',
                     'input_mol_block',
                     'FF_minimised_vina_mininplace_docked_mol_block',
                     'FF_minimised_glide_inplace_docked_mol_block', 
                     'FF_minimised_glide_mininplace_docked_mol_block',
                     'FF_unminimised_glide_inplace_docked_mol_block',
                     'FF_unminimised_glide_mininplace_docked_mol_block',
                     'FF_minimised_vina_inplace_docked_pdbqt_mol_block',
                     'FF_minimised_vina_mininplace_docked_pdbqt_mol_block',
                     'FF_unminimised_vina_inplace_docked_pdbqt_mol_block',
                     'FF_unminimised_vina_mininplace_docked_pdbqt_mol_block'], 
                     inplace=True)
    df.to_csv('please_remove_refactored.csv')
