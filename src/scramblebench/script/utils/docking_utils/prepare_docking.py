import os
import subprocess
from pathlib import Path

def calculate_easydock_cpu():
    if hasattr(os, 'sched_getaffinity'):
        num_cpus = len(os.sched_getaffinity(0))
    else:
        num_cpus = os.cpu_count() or 1
    CPU_FOR_VINA = 5
    CPU_BUFFER = 5

    return max(1, int((num_cpus - CPU_BUFFER)/ CPU_FOR_VINA))



def fetch_glide_sp_cpu(schrodinger_dir):

    cmd = [f'{schrodinger_dir}/run', 'lictool', 'status']

    schrodinger_license_strings = subprocess.run(cmd, capture_output=True, text=True).stdout
    glide_sp_license_available = int(schrodinger_license_strings.split('GLIDE_SP_DOCKING')[-1].split('\n')[0].split('licenses available')[0].split('Total of')[-1])

    MAX_CPU_USED = 50
    CPU_BUFFER = 5
    GLIDE_SP_LICENSE_PER_CPU = 4

    if hasattr(os, 'sched_getaffinity'):
        num_cpus = len(os.sched_getaffinity(0))
    else:
        num_cpus = os.cpu_count() or 1
    
    AVAILABLE_CPU = num_cpus - CPU_BUFFER

    return max(1, min(int(glide_sp_license_available / GLIDE_SP_LICENSE_PER_CPU) - CPU_BUFFER, MAX_CPU_USED, AVAILABLE_CPU))


def fetch_ligprep_cpu(schrodinger_dir):

    cmd = [f'{schrodinger_dir}/run', 'lictool', 'status']

    schrodinger_license_strings = subprocess.run(cmd, capture_output=True, text=True).stdout
    glide_sp_license_available = int(schrodinger_license_strings.split('LIGPREP_MAIN')[-1].split('\n')[0].split('licenses available')[0].split('Total of')[-1])

    MAX_CPU_USED = 50
    GLIDE_SP_LICENSE_PER_CPU = 1
    CPU_BUFFER = 5

    if hasattr(os, 'sched_getaffinity'):
        num_cpus = len(os.sched_getaffinity(0))
    else:
        num_cpus = os.cpu_count() or 1

    AVAILABLE_CPU = num_cpus - CPU_BUFFER

    return max(1, min(int(glide_sp_license_available / GLIDE_SP_LICENSE_PER_CPU), MAX_CPU_USED, AVAILABLE_CPU))


def prepare_ligprep_inp_file(ligand_fname, output_fname, write_fname):

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

    with open(write_fname, 'w') as ligprep_fname:
        ligprep_fname.write(ligprep_inp_data)


def prepare_easydock_vina_parameter(config_data):

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


def prepare_easydock_grid(pocket_coordinates, write_fname):
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

    with open(write_fname, 'w') as grid_f:
        grid_f.write(grid_txt)


def prepare_glide_inp_file(grid_fname, ligand_fname, write_fname, intra_hbonds=False):
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

    with open(write_fname, 'w') as glide_inp_f:
        glide_inp_f.write(glide_inp)
