"""Calculate conformation validity of de novo ligand."""

import argparse
import logging
import sys
import subprocess

from pathlib import Path
from copy import deepcopy
from collections import defaultdict
from typing import Any

import yaml
from scramblebench.script.config_preparation import config_constant, config_genbench3d, config_input, config_parameter
from scramblebench.script.utils.process_data import read_input, fetch_model_file_from_model_dir
from scramblebench.script.utils.job_manager import CheckpointManager, CheckpointStatus

logger = logging.getLogger(__name__)


def check_schrodinger_status(schrodinger_dir: str) -> None:
    """check if schrodinger jsc is running. Will check for schrodinger2025-4 onwards.
    JSC cannot accept jobs if the local-server has not been started.

    Args:
        schrodinger_dir (str): schrodinger rootdir

    Raises:
        ValueError: if root directory is not schrodinger's
    """
    if 'schrodinger' not in Path(schrodinger_dir).name:
        raise ValueError('schrodinger dirname should end in schrodinger string \
                         (e.g., /path/to/schrodinger-year-quarter)')

    schrodinger_version_year = str(schrodinger_dir.split('schrodinger')[-1].split('-')[0])
    schrodinger_version_quarter = str(schrodinger_dir.split('schrodinger')[-1].split('-')[-1])

    cmd = [f'{str(schrodinger_dir)}/jsc', 'local-server-status']

    # schrodinger jsc is implemented from 2025-4 onwards
    if int(schrodinger_version_year) < 2025:
        return None
    if schrodinger_version_year == '2025' and schrodinger_version_quarter != '4':
        return None

    schrodinger_status = subprocess.run(cmd, capture_output=True, text=True, check=True).stdout
    if 'STOPPED' in schrodinger_status:
        subprocess.run([f'{str(schrodinger_dir)}/jsc', 'local-server-start'],
                       capture_output=True, text=True, check=True)

    return None


def prepare_genbench3d_input(input_data: config_input.InputConfig,
                             genbench3d_data: config_genbench3d.GenBench3DConfig) -> tuple[str, str]:
    """move the input folder to a new subfolder inside the genbench3d folder

    Args:
        input_data (config_input.InputConfig): InputConfig
        genbench3d_data (config_genbench3d.GenBench3DConfig): GenBench3DConfig

    Returns:
        tuple[str, str]: new input protein file and ligand file
    """

    GENBENCH3D_INPUT_DIR = 'structural_input'
    genbench3d_input_dirpath = Path(genbench3d_data.output_value) / GENBENCH3D_INPUT_DIR
    genbench3d_input_dirpath.mkdir(exist_ok=True, parents=True)

    input_pdb = str(genbench3d_input_dirpath / Path(input_data.pdb_value).name)
    input_sdf = str(genbench3d_input_dirpath / Path(input_data.sdf_value).name)
    subprocess.run(['cp', input_data.pdb_value, genbench3d_input_dirpath], text=True, check=True)
    subprocess.run(['cp', input_data.sdf_value, genbench3d_input_dirpath], text=True, check=True)

    return input_pdb, input_sdf


def initialise_genbench3d_checkpoint(checkpoint_fname: str,
                                     models: list[str],
                                     minimisations: list[str]) -> CheckpointManager:
    """initialize the checkpoint manager for genbench3d

    Args:
        checkpoint_fname (str): checkpoint filename to be saved
        models (list[str]): name of models used to do genbench3d
        minimisations (list[str]): list of 'unminimised' or 'minimised' 

    Returns:
        CheckpointManager: class to track genbench3d progress for each model's parameter
    """
    checkpoint_manager = CheckpointManager()
    checkpoint_filepath = Path(checkpoint_fname)

    if not checkpoint_filepath.is_file():
        checkpoint = defaultdict(dict)
        for model in models:
            for minimisation in minimisations:
                checkpoint[model][minimisation] = CheckpointStatus.PENDING.value

        checkpoint_manager = checkpoint_manager.from_dict(checkpoint)
        checkpoint_manager.checkpoint_fname = checkpoint_filepath
        checkpoint_manager.save_state()
    else:
        checkpoint_manager = checkpoint_manager.from_json(checkpoint_filepath)

    return checkpoint_manager


def prepare_genbench3d_shared_cmd(input_data: config_input.InputConfig,
                                  genbench3d_data: config_genbench3d.GenBench3DConfig) -> list[str]:
    """prepare genbench3d command line that will be shared to all genbench3d parameter

    Args:
        input_data (config_input.InputConfig): InputConfig
        genbench3d_data (config_genbench3d.GenBench3DConfig): GenBench3DConfig

    Returns:
        list[str]: command line to run genbench3d
    """

    if input_data:
        input_pdb, input_sdf = prepare_genbench3d_input(input_data=input_data,
                                                        genbench3d_data=genbench3d_data)

        cmd  = ["conda", "run", "-n", genbench3d_data.conda_env_value,
                'python', str(Path(genbench3d_data.genbench3d_dir_value) / 'sb_benchmark_mols.py'),
                '-c', genbench3d_data.genbench3d_config_value,
                '-p', input_pdb,
                '-n', input_sdf,
                '--do_conf_analysis']

        if genbench3d_data.do_cancel_protonation_by_genbench3d_value:
            cmd += ['--cancel_protonation']

    else:
        cmd = ["conda", "run", "-n", genbench3d_data.conda_env_value,
                'python', str(Path(genbench3d_data.genbench3d_dir_value) / 'benchmark_mols.py'),
                '-c', genbench3d_data.genbench3d_config_value]

    return cmd


def prepare_genbench3d_model_cmd(model: str,
                                 fname: str,
                                 checkpoint_manager: CheckpointManager,
                                 genbench3d_data: config_genbench3d.GenBench3DConfig,
                                 cmd: list[str]) -> list[dict[str, Any]]:
    """prepare the genbench3d model command lines

    Args:
        model (str): model name
        fname (str): input sdf file
        checkpoint_manager (CheckpointManager): CheckpointManager
        genbench3d_data (config_genbench3d.GenBench3DConfig): GenBench3DConfig
        cmd (list[str]): shared or global cmd created in prepare_genbench3d_shared_cmd()

    Returns:
        list[dict[str, Any]]: list containing information as dictionary to run genbench3d.
    """

    do_sbdd_analysis = False
    if '-p' in cmd and '-n' in cmd:
        do_sbdd_analysis = True

    json_output_dirpath = Path(genbench3d_data.output_value) / 'json_output'
    json_output_dirpath.mkdir(parents=True, exist_ok=True)

    cmd_list = []
    for minimisation in checkpoint_manager.state[model].keys():

        if checkpoint_manager.state[model][minimisation] == CheckpointStatus.COMPLETED.value:
            logging.info(f'{fname} {model} has been completed with {minimisation} config. \
                            Skipping...')
            continue

        model_cmd = deepcopy(cmd)
        fname_stem = Path(fname).stem

        model_cmd += ['-i', fname]
        model_cmd += ['-o', f'{json_output_dirpath}/{fname_stem}_{minimisation}.json']
        model_cmd += ['--log_output',
                      f'{json_output_dirpath}/{fname_stem}_{minimisation}{"_sb" if do_sbdd_analysis else ""}_benchmark.log']

        if do_sbdd_analysis:
            if minimisation == 'minimised':
                model_cmd += ['-m']

            if genbench3d_data.schrodinger_dir_value:
                check_schrodinger_status(genbench3d_data.schrodinger_dir_value)
                glide_output =  Path(genbench3d_data.output_value) / 'Glide' / model / minimisation
                glide_output.mkdir(parents=True, exist_ok=True)

                model_cmd += ['--glide', '--output_glide_dir', str(glide_output)]

            vina_output = Path(genbench3d_data.output_value) / 'Vina' / model / minimisation
            vina_output.mkdir(parents=True, exist_ok=True)
            model_cmd += ['--output_vina_dir', str(vina_output)]

        cmd_data = {}
        cmd_data['minimisation'] = minimisation
        cmd_data['model'] = model
        cmd_data['cmd'] = model_cmd
        cmd_list.append(cmd_data)

    return cmd_list


def prepare_genbench3d_cmd(config_data: dict[str, dict]) -> tuple[list[list[str]], str]:
    """create a list of command lines to run genbench3d along with the checkpoint file

    Args:
        config_data (dict[str, dict]): user's prepared YAML config as dictionary

    Returns:
        tuple[list[list[str]], str]: list of command lines to run genbench3d along with the checkpoint file
    """
    genbench3d_data = config_genbench3d.GenBench3DConfig(config_data[config_constant.ANALYSIS_KEY])
    genbench3d_data.validate_config()

    parameter_data = config_parameter.ParameterConfig(config_data=config_data)
    valid_molecule_file_dict = fetch_model_file_from_model_dir(dir_path=genbench3d_data.input_value,
                                                               model_list=parameter_data.model_list_value)

    complex_minimisation = ['unminimised']

    input_data = False
    if config_constant.INPUT_KEY in config_data:
        input_data = config_input.InputConfig(config_data)
        if genbench3d_data.do_complex_forcefield_minimisation_value:
            complex_minimisation.append('minimised')

    cmd = prepare_genbench3d_shared_cmd(input_data=input_data,
                                        genbench3d_data=genbench3d_data)

    checkpoint_filepath = Path(genbench3d_data.output_value) / 'genbench3d_checkpoint.json'
    checkpoint_manager = initialise_genbench3d_checkpoint(checkpoint_fname=checkpoint_filepath,
                                                        models=valid_molecule_file_dict.keys(),
                                                        minimisations=complex_minimisation)

    cmd_list = []
    for model, fname in valid_molecule_file_dict.items():
        cmd_list += prepare_genbench3d_model_cmd(checkpoint_manager=checkpoint_manager,
                                                 model=model,
                                                 genbench3d_data=genbench3d_data,
                                                 cmd=cmd,
                                                 fname=fname)

    return cmd_list, checkpoint_filepath


def run_single_genbench3d(cmd_data: dict[str, Any]) -> tuple[str, str, str]:
    """to execute genbench3d on a single parameter. checkpoint is not directly saved here, because in the
    case where you want to execute multiprocessing, you cannot share the checkpoint file together. Hence,
    the checkpoint will be updated in another function.

    Args:
        cmd_data (dict[str, Any]): a small dictionary containing the model name, minimisation parameter,
                                   and the command line to run

    Returns:
        tuple[str, str, str]: job status, model name, minimisation parameter
    """

    try:
        model = cmd_data['model']
        minimisation = cmd_data['minimisation']
        cmd = cmd_data['cmd']

        logging.info(f"Running Genbench. cmd: {cmd}")
        subprocess.run(cmd, capture_output=True, text=True, check=True)

        output_fname = cmd[cmd.index('-o') + 1]

        if Path(output_fname).is_file():
            return CheckpointStatus.COMPLETED.value, model, minimisation
        return CheckpointStatus.FAILED.value, model, minimisation

    except (subprocess.CalledProcessError, PermissionError, KeyboardInterrupt) as e:
        logging.exception(e)
        return CheckpointStatus.FAILED.value, model, minimisation


def run_genbench3d(config_data: dict[str, dict]):
    """main function of p4_analyse_genbench3d.py. Initial attempts were made to execute multiprocessing,
    but tracking Glide SP licenses are difficult here. So, I decided to just use single CPU.

    Args:
        config_data (dict[str, dict]): user's prepared YAML config as dictionary
    """
    cmd_list, checkpoint_fname = prepare_genbench3d_cmd(config_data=config_data)
    checkpoint_manager = CheckpointManager().from_json(checkpoint_fname)

    # MAX_CPU_USED = 5
    # CPU_BUFFER = 5
    # num_cpu_available = len(os.sched_getaffinity(0)) - CPU_BUFFER
    # used_cpu = max(1, min(MAX_CPU_USED, num_cpu_available))
    # with Pool(used_cpu) as pool:
    #     for genbench3d_status, model, minimisation in pool.imap_unordered(run_single_genbench3d, cmd_list):

    for cmd in cmd_list:
        genbench3d_status, model, minimisation = run_single_genbench3d(cmd)
        checkpoint_manager.state[model][minimisation] = genbench3d_status
        checkpoint_manager.save_state()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Perform genbench3d analysis after molecule preparation")

    parser.add_argument("-i", "--input", help="config yaml input file or txt file containing yaml filepath",
                        required=True, type=str)
    args = parser.parse_args()

    logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='%(asctime)s - %(module)s: - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

    logging.info('Running p4_analyse_genbench3d.py')
    logging.info('Reading the config filename :)')


    yaml_file_list = read_input(args.input)

    for yaml_file in yaml_file_list:
        checkpoint_json = Path(yaml_file).parent
        with open(yaml_file, 'r', encoding='utf-8') as yaml_f:
            data_input = yaml.safe_load(yaml_f)
        run_genbench3d(config_data=data_input)
