"""This file util handles docking preparation for p4_analyse_redocking.py"""
import os
import subprocess
from pathlib import Path
from collections import defaultdict
from typing import Any

import rdkit
from rdkit import Chem

MAX_CPU_USED = 50

CPU_BUFFER = 5
CPU_FOR_VINA = 5

LIGPREP_LICENSE_PER_CPU = 1
GLIDE_SP_LICENSE_PER_CPU = 4

def calculate_easydock_cpu() -> int:
    """calculate available process for easydock. cpu for each vina process is set to 5
    based on my personal experience.

    Returns:
        int: available process
    """
    if hasattr(os, 'sched_getaffinity'):
        num_cpus = len(os.sched_getaffinity(0))
    else:
        num_cpus = os.cpu_count() or 1

    return max(1, int((num_cpus - CPU_BUFFER)/ CPU_FOR_VINA))



def fetch_glide_sp_cpu(schrodinger_dir: str) -> int:
    """calculate available cpu for Glide SP. Each process consumes 4? license token

    Args:
        schrodinger_dir (str): root directory for schrodinger

    Returns:
        int: number of cpu for Glide SP
    """

    cmd = [f'{schrodinger_dir}/run', 'lictool', 'status']

    schrodinger_license_strings = subprocess.run(cmd, capture_output=True, check=True, text=True).stdout
    glide_sp_license_available = int(schrodinger_license_strings.split('GLIDE_SP_DOCKING')[-1]
                                                                .split('\n')[0]
                                                                .split('licenses available')[0]
                                                                .split('Total of')[-1])

    if hasattr(os, 'sched_getaffinity'):
        num_cpus = len(os.sched_getaffinity(0))
    else:
        num_cpus = os.cpu_count() or 1

    available_cpu = num_cpus - CPU_BUFFER

    return max(1, min(int(glide_sp_license_available / GLIDE_SP_LICENSE_PER_CPU) - CPU_BUFFER,
                      MAX_CPU_USED, available_cpu))


def fetch_ligprep_cpu(schrodinger_dir: str) -> int:
    """calculate available cpu for LigPrep. Each process consumes 1 license token 

    Args:
        schrodinger_dir (str): root directory for schrodinger

    Returns:
        int: number of cpu for LigPrep
    """

    cmd = [f'{schrodinger_dir}/run', 'lictool', 'status']

    schrodinger_license_strings = subprocess.run(cmd, capture_output=True, check=True, text=True).stdout
    glide_sp_license_available = int(schrodinger_license_strings.split('LIGPREP_MAIN')[-1]
                                                                .split('\n')[0]
                                                                .split('licenses available')[0]
                                                                .split('Total of')[-1])



    if hasattr(os, 'sched_getaffinity'):
        num_cpus = len(os.sched_getaffinity(0))
    else:
        num_cpus = os.cpu_count() or 1

    available_cpu = num_cpus - CPU_BUFFER

    return max(1, min(int(glide_sp_license_available / LIGPREP_LICENSE_PER_CPU),
                      MAX_CPU_USED, available_cpu))


def prepare_ligprep_inp_file(ligand_fname: str, output_fname: str, write_fname: str):
    """prepare the .inp file for ligprep job. The input format is valid as of 2025-4 version.

    Args:
        ligand_fname (str): .sdf file for ligand name
        output_fname (str): .sdf file for ligprep output name
        write_fname (str): .inp output file
    """

    ligprep_inp_data = f'''
    INPUT_FILE_NAME   {Path(ligand_fname).name}
    MAX_ATOMS   500
    FORCE_FIELD   16
    PH   7.4
    PH_THRESHOLD   1.0
    EPIK   yes
    EPIKX   no
    EPIK_METAL_BINDING   no
    INCLUDE_ORIGINAL_STATE   no
    DETERMINE_CHIRALITIES   no
    IGNORE_CHIRALITIES   no
    NUM_STEREOISOMERS   32
    OUT_SD  {Path(output_fname).name}

    '''

    with open(write_fname, 'w', encoding='utf-8') as ligprep_fname:
        ligprep_fname.write(ligprep_inp_data)


def prepare_easydock_vina_parameter(config_data: dict[str, Any]) -> dict[str, Any]:
    """to prepare or validate the easydock.yml config file

    Args:
        config_data (dict[str, Any]): easydock .yml config file as dictionary

    Returns:
        dict[str, Any]: prepared easydock .yml config file as dictionary
    """

    config_data['exhaustiveness'] = config_data.get('exhaustiveness') or 32
    config_data['n_poses'] = config_data.get('n_poses') or 5
    config_data['seed'] = config_data.get('seed') or 32
    config_data['exhaustiveness'] = config_data.get('exhaustiveness') or 42

    if config_data.get('ncpu') is None:
        if hasattr(os, 'sched_getaffinity'):
            num_cpus = len(os.sched_getaffinity(0))
        else:
            num_cpus = os.cpu_count() or 1
        config_data['ncpu'] = max(1, min(num_cpus - 5, 5))

    return config_data


def prepare_easydock_grid(pocket_coordinates: list[float] | str, write_fname: str):
    """prepare the grid .txt file for easydock.

    Args:
        pocket_coordinates (list[float] | str): pocket coordinates of the protein
        write_fname (str): .txt output file
    """
    if isinstance(pocket_coordinates, str):
        pocket_coordinates = pocket_coordinates.split(',')

    assert isinstance(pocket_coordinates, list) and len(pocket_coordinates) == 3

    grid_txt = f'''
    center_x = {pocket_coordinates[0]}
    center_y = {pocket_coordinates[1]}
    center_z = {pocket_coordinates[2]}
    size_x = 25.0
    size_y = 25.0
    size_z = 25.0
    '''

    with open(write_fname, 'w', encoding='utf-8') as grid_f:
        grid_f.write(grid_txt)


def prepare_glide_inp_file(grid_fname: str, ligand_fname: str,
                           write_fname: str, intra_hbonds:bool=False):
    """prepare the .inp file for Glide SP

    Args:
        grid_fname (str): grid file
        ligand_fname (str): .sdf ligand file
        write_fname (str): .inp output file
        intra_hbonds (bool, optional): whether to reward intramolecular hydrogen bonds.
        Defaults to False.
    """
    glide_inp = f'''
        FORCEFIELD   OPLS_2005
        FORCEPLANAR   True
        GRIDFILE   {grid_fname}
        LIGANDFILE   {Path(ligand_fname).name}
        POSE_OUTTYPE   ligandlib_sd
        POSES_PER_LIG   5
        POSTDOCK_NPOSE   20
        POSTDOCKSTRAIN   True
        PRECISION   SP
        REWARD_INTRA_HBONDS {intra_hbonds}
        WRITE_VSDB   False
    '''

    with open(write_fname, 'w', encoding='utf-8') as glide_inp_f:
        glide_inp_f.write(glide_inp)


def prepare_easydock_ligand_input(input_fname: str, output_dir: str) -> str:
    """prepare the .sdf ligand for easydock docking. This is because sometimes ligands that
    are prepared externally can have multiple stereoisomer or protonation states, which will
    trigger an error to the easydock .DB in the past. Although this is fixed by only including
    one isomer and deleting the rest, the user may not want this. Hence, this is a workaround
    to append the molecule name (or id in easydock) depending on the number of isomer/protonation
    states to prevent the unique constraint of .db schema.

    Of course, adding a number to the user's molecule name means we need to remove them again in the
    final process. Please refer to the process_easydock_docking_output function.

    Args:
        input_fname (str): .sdf ligand files
        output_dir (str): output directory

    Returns:
        str: prepared .sdf ligand files
    """

    mol_l = Chem.SDMolSupplier(input_fname, removeHs=False)

    output_fname = Path(output_dir) / f'{Path(input_fname).stem}_ez_prepared.sdf'
    write_fname = Chem.SDWriter(output_fname)
    mol_dict = {}

    for mol in mol_l:
        mol_id = mol.GetProp('_Name')
        if mol_id in mol_dict:
            mol_dict[mol_id] += 1
        else:
            mol_dict[mol_id] = 0

        mol.SetProp('_Name', f'{mol_id}_{mol_dict[mol_id]}')
        write_fname.write(mol)

    return output_fname


def process_easydock_docking_output(input_fname: str) -> str:
    """This is to remove the additional number suffix tags to duplicates/multiple molecule names
    from the prepare_easydock_ligand_input function. 

    Args:
        input_fname (str): input .sdf file

    Returns:
        str: processed .sdf file
    """

    mol_l = Chem.SDMolSupplier(input_fname, removeHs=False)
    output_fname = Path(input_fname).parent / f'{Path(input_fname).stem}_processed.sdf'
    write_fname = Chem.SDWriter(output_fname)
    mol_dict = defaultdict(dict)

    for mol in mol_l:
        docking_score = float(mol.GetProp('docking_score'))
        parent_id = mol.GetProp('_Name').rsplit('_', 1)[0]
        if parent_id not in mol_dict:
            mol_dict[parent_id]['mol'] = mol
            mol_dict[parent_id]['docking_score'] = docking_score

        else:
            if docking_score < mol_dict[parent_id]['docking_score']:
                mol_dict[parent_id]['mol'] = mol
                mol_dict[parent_id]['docking_score'] = docking_score

    for data in mol_dict.values():
        mol = data['mol']
        mol.SetProp('_Name', mol.GetProp('_Name').rsplit('_', 1)[0])
        write_fname.write(mol)

    return output_fname
