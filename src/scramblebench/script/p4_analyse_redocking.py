"""perform protonation, docking, and PLIF of the de novo ligand.
Most of the functions use temporary directory to prevent user from seeing intermediates files and
because some of the program (Schrodinger) needs files to be in one folder, and the current working
directory is preferred to be in that specific folder as well. Hence, please be patient if the speed
is suboptimal."""
import argparse
import logging
import sys
import subprocess
import os
import tempfile
import gzip
import shutil

from typing import Optional
from pathlib import Path
from copy import deepcopy
from collections import defaultdict

import rdkit
import yaml
import pandas as pd

from rdkit import Chem

from scramblebench.script.config_preparation import config_constant, config_input, config_redocking, config_parameter

from scramblebench.script.utils.process_data import read_input, fetch_model_file_from_model_dir
from scramblebench.script.utils.docking_utils.prepare_protein import VinaProtein, GlideProtein
from scramblebench.script.utils.docking_utils.prepare_docking import calculate_easydock_cpu, fetch_glide_sp_cpu,\
                                                        fetch_ligprep_cpu, prepare_easydock_grid,\
                                                        prepare_easydock_vina_parameter, prepare_ligprep_inp_file,\
                                                        prepare_glide_inp_file, prepare_easydock_ligand_input,\
                                                        process_easydock_docking_output
from scramblebench.script.utils.job_manager import CheckpointManager, CheckpointStatus

logger = logging.getLogger(__name__)


def run_protonation_easydock(protonation_data: config_redocking.RedockingProtonationConfig,
                             valid_molecule_file_dict: dict[str, str]) -> None:
    """run protonation for easydock. No docking will be done. Here, I am running using a tempdir
    as I think that the .db file is something that people may not want to see.

    Args:
        protonation_data (config_redocking.RedockingProtonationConfig): RedockingProtonationConfig
        valid_molecule_file_dict (ItemsView): dictionary of model: ligand filename
    """

    cmd = ['conda', 'run', '-n', protonation_data.environment_value,
            'easydock', '-c', str(calculate_easydock_cpu())]

    for model, fname in valid_molecule_file_dict.items():

        protonation_output_dir = Path(Path(protonation_data.output_value) / model)
        protonation_output_dir.mkdir(parents=True, exist_ok=True)
        output_easydock_sdf = str(f'{protonation_output_dir / Path(fname).stem}_protonated.sdf')
        if Path(output_easydock_sdf).is_file():
            logging.info(f'Protonation with {output_easydock_sdf} has been done. Skipping...')
            continue

        with tempfile.TemporaryDirectory() as temp_protonation_dir:
            output_easydock_db = str(Path(temp_protonation_dir) / f'{Path(fname).stem}_protonated.db')
            model_cmd = deepcopy(cmd)
            model_cmd += ['-i', fname, '-o', output_easydock_db, '--protonation', protonation_data.method_value]

            logging.info(f'Running easydock protonation with cmd: {model_cmd=}')
            subprocess.run(model_cmd, text=True, check=True)

            Path(Path(protonation_data.output_value) / model).mkdir(parents=True, exist_ok=True)
            FETCH_SDF_SCRIPT = Path(__file__).parent / 'utils/docking_utils' / 'easydock_fetch_sdf_from_db.py'
            protonated_cmd = ['conda', 'run', '-n', protonation_data.environment_value,
                            'python', FETCH_SDF_SCRIPT, '-i', output_easydock_db, '-o', output_easydock_sdf]

            subprocess.run(protonated_cmd, text=True, check=True)


def prepare_docking_fnames_easydock(docking_data: config_redocking.EasyDockConfig,
                                    input_data: config_input.InputConfig,
                                    output_dir: str) -> str:
    """prepare the necessary filenames required for easydock, namely config files, grid files.
    Here, I have only set up vina program for docking. The function assumes that the user have
    used the correct config.yml for easydock, and this function will only update the protein and
    grid file for each docking program. In the case where this function needs to support more docking
    program, a class may be needed.

    Args:
        docking_data (config_redocking.EasyDockConfig): EasyDockConfig
        input_data (config_input.InputConfig): InputConfig
        output_dir (str): directory for the new config.yml for easydock

    Returns:
        str: filename for the new config.yml file
    """
    config_fname = docking_data.config_value
    protein_pdb = input_data.pdb_value
    grid_fname = str(Path(protein_pdb).parent / f'{input_data.protein_value}_grid.txt')

    prepare_easydock_grid(pocket_coordinates=input_data.pocket_coord_value,
                          write_fname=grid_fname)


    with open(config_fname, 'r', encoding='utf-8') as config_f:
        config_data = yaml.load(config_f, Loader=yaml.SafeLoader)

        prepared_protein = VinaProtein(pdb_filepath=protein_pdb,
                                        prepare_receptor_bin_path=docking_data.protein_pdbqt_executable_value,
                                        protonation_method=docking_data.protein_preparation_value,
                                        preparation_method=docking_data.protein_pdbqt_preparation_value)

        if docking_data.docking_value != 'server':

            config_data['protein'] = prepared_protein.pdbqt_filepath
            config_data['protein_setup'] = grid_fname

        elif config_data['init_server'].get('program') == 'vina-gpu':
            config_data['init_server']['protein'] = prepared_protein.pdbqt_filepath
            config_data['init_server']['protein_setup'] = grid_fname

        else:
            config_data['init_server']['protein'] = prepared_protein.protein_clean_filepath

    new_config_fname = Path(output_dir) / f'easydock_{Path(input_data.protein_value).stem}.yml'
    if str(docking_data.docking_value).lower().strip() == 'vina':
        config_data = prepare_easydock_vina_parameter(config_data=config_data)
    else:
        logging.warning(f'You are using a different docking program other than vina. \
                        Ensure that the config at {new_config_fname} has the correct format.\
                        Please check easydock github for reference!')


    with open(new_config_fname, 'w', encoding='utf-8') as config_f:
        yaml.dump(config_data, config_f)

    return new_config_fname


def run_plif_easydock(environment: str, db_fname: str, protein_fname: str,
                      reference_fname: str, reference_plif: Optional[list[str]]=None) -> str:
    """run easydock_plif after docking. This function will only works if docking has been done.
    Since the input files are based on the .db files, not the .sdf files.

    Args:
        environment (str): conda environment to run easydock
        db_fname (str): the intermediates .db files storing all information regarding easydock docking
        protein_fname (str): PDB file for protein
        reference_fname (str): SDF file for the reference ligand (not the de novo ligand)
        reference_plif (Optional[list[str]], optional): reference interaction that needs to be matched
        with the de novo generated ligands. Must be in a list with format 
        [residue_name][residue_number].[chain].[interaction_type]. Defaults to None.

    Returns:
        str: CSV file containing 'id', 'stereo_id', 'pose', 'plif_sim' columns
    """

    if not isinstance(reference_plif, list):

        reference_output = Path(db_fname).parent / 'ref_similarity.csv'
        subprocess.run(['conda', 'run', '-n', environment,
                        'easydock_plif', '-i', reference_fname, '-p', protein_fname,
                        '-o', reference_output, '-c', str(4)], check=True)

        reference_plif = [col for col in \
                          pd.read_csv(reference_output, delimiter='\t').columns.tolist() if col != 'Name']

    plif_result = Path(db_fname).parent / 'similarity.csv'
    cmd = ['conda', 'run', '-n', environment,
           'easydock_plif', '-i', db_fname, '-p', protein_fname, '-o',
            plif_result, '--ref_plif']

    cmd += reference_plif

    logging.info(f'Running easydock PLIF with cmd: {cmd}')
    subprocess.run(cmd, check=True)

    return plif_result


def analyse_plif_easydock(plif_result: str, db_fname: str,
                          similarity: float | int, environment: str) -> Optional[str]:
    """analyse the plif output by only saving ligands that have passed the plif_similarity
    filter based on the reference ligand.

    Args:
        plif_result (str): CSV file containing 'id', 'stereo_id', 'pose', 'plif_sim' columns
        db_fname (str): the intermediates .db files storing all information regarding easydock docking
        similarity (float | int): PLIF similarity threshold ranging from 0 to 1
        environment (str): conda environment to run easydock

    Returns:
        Optional[str]: sdf output containing ligands that pass the plif_similarity filter.
                       each ligand has plif similarity and docking score result
    """

    plif_df = pd.read_csv(plif_result, delimiter='\t').set_index('id')

    plif_df = plif_df[plif_df['plif_sim'] >= similarity]

    total_mol = []
    for pose, df in plif_df.groupby('pose'):
        output_pose_sdf = Path(db_fname).parent / f'filtered_sdf_{pose}.sdf'

        subprocess.run(['conda', 'run', '-n', environment,
                        'get_sdf_from_easydock', '-i', str(db_fname),
                        '-d'] + df.index.tolist() + ['--poses', str(pose), # assuming only 1 stereoisomer
                        '-o', str(output_pose_sdf), '--fields', 'docking_score'], check=True)

        total_mol += list(Chem.SDMolSupplier(output_pose_sdf, removeHs=False))


    if len(total_mol) > 0:
        temp_output = Path(db_fname).parent / 'filtered_sdf.sdf'
        with Chem.SDWriter(temp_output) as output_f:
            for mol in total_mol:
                similarity_score = plif_df.loc[mol.GetProp('_Name').rsplit(':', 1)[0], 'plif_sim']
                mol.SetProp('plif_similarity', str(similarity_score))
                output_f.write(mol)

        processed_temp_output = process_easydock_docking_output(temp_output)
        return processed_temp_output

    logging.warning('No molecule matching with PLIF')
    return None


def generate_shared_cmd_docking_easydock(docking_data: config_redocking.EasyDockConfig,
                                         input_data: config_input.InputConfig) -> list[str]:
    """create a shared command line that will be used to run easydock

    Args:
        docking_data (config_redocking.EasyDockConfig): EasyDockConfig
        input_data (config_input.InputConfig): InputConfig

    Returns:
        list[str]: command line to run easydock as a list of string
    """

    cmd = ['conda', 'run', '-n', docking_data.environment_value,
            'easydock',  '--sdf']

    Path(docking_data.output_value).mkdir(parents=True, exist_ok=True)

    config_fname = prepare_docking_fnames_easydock(docking_data, input_data, output_dir=docking_data.output_value)
    cmd += ['--config', str(config_fname)]


    if docking_data.docking_value:
        cmd += ['--program', docking_data.docking_value]
    else:
        logging.warning('You did not specify which docking program to use for easydock')

    if docking_data.cpu_value:
        if hasattr(os, 'sched_getaffinity'):
            num_cpus = len(os.sched_getaffinity(0))
        else:
            num_cpus = os.cpu_count() or 1
        cmd += ['-c', str(min(int(docking_data.cpu_value), num_cpus - 5))]

    else:
        cmd += ['-c', str(calculate_easydock_cpu())]

    return cmd


def run_model_docking_easydock(fname: str, model: str,
                               docking_data: config_redocking.EasyDockConfig,
                               temp_docking_dir: str, cmd: list[str]) -> str:
    """run a single easydock redocking job for a single model. Temporary directory is used to hide
    the .db database

    Args:
        fname (str): sdf input containing valid de novo generated ligand
        model (str): model name used to generate the ligand
        docking_data (config_redocking.EasyDockConfig): EasyDockConfig
        temp_docking_dir (str): temporary directory 
        cmd (list[str]): shared command line that is created by generate_shared_cmd_docking_easydock().

    Returns:
        str: the intermediates .db files storing all information regarding easydock docking
    """

    temp_output_easydock_db = str(Path(temp_docking_dir) / f'{Path(fname).stem}_redocked.db')
    model_cmd = deepcopy(cmd)

    preprocessed_fname = deepcopy(fname)
    if docking_data.protonation_value:
        if 'schrodinger' not in docking_data.protonation_value:
            model_cmd += ['--protonation', docking_data.protonation_value]
        else:
            preprocessed_fname = run_protonation_ligprep(schrodinger_dir=docking_data.protonation_value,
                                    output_dir=temp_docking_dir,
                                    valid_molecule_file_dict={model: fname})

    easydock_input_fname = prepare_easydock_ligand_input(input_fname=preprocessed_fname,
                                                            output_dir=temp_docking_dir)

    model_cmd += ['-i', str(easydock_input_fname), '-o', temp_output_easydock_db]

    logging.info(f'Easydock docking with cmd: {model_cmd=}')
    subprocess.run(model_cmd, text=True, check=True)

    return temp_output_easydock_db


def execute_docking_easydock(docking_data: config_redocking.EasyDockConfig,
                             parameter_data: config_parameter.ParameterConfig,
                             input_data: config_input.InputConfig) -> None:
    """run easydock redocking. Temporary directory is used to hide the .db information. I am not sure
    if people like to see the .db file, as it requires user to have a DB viewer/ SQL knowledge. Since the .sdf
    is the most important thing, I have only decided to show the .sdf file output to user.

    Args:
        docking_data (config_redocking.EasyDockConfig): EasyDockConfig
        parameter_data (config_parameter.ParameterConfig): ParameterConfig
        input_data (config_input.InputConfig): InputConfig
    """


    checkpoint_manager = initialise_redocking_checkpoint(docking_data=docking_data,
                                                         models=parameter_data.model_list_value)

    for model, fname in fetch_model_file_from_model_dir(dir_path=docking_data.input_value,
                                                        model_list=parameter_data.model_list_value).items():

        output_easydock_sdf = str(Path(docking_data.output_value) / model / f'{Path(fname).stem}_easydock_redocked.sdf')
        Path(output_easydock_sdf).parent.mkdir(parents=True, exist_ok=True)

        if docking_data.plif_value:
            output_easydock_plif = str(Path(docking_data.output_value) / model / 'plif' / f'{Path(fname).stem}_easydock_plif_processed.sdf')
            Path(output_easydock_plif).parent.mkdir(parents=True, exist_ok=True)
            if checkpoint_manager.state[model][docking_data.name]['plif'] == CheckpointStatus.COMPLETED.value:
                logging.info(f'Easydock Docking and PLIF with {output_easydock_sdf} has been done. Skipping...')
                continue

        if Path(output_easydock_sdf).is_file() and not docking_data.plif_value:
            logging.info(f'Easydock Docking with {output_easydock_sdf} has been done. Skipping...')
            continue

        try:
            with tempfile.TemporaryDirectory() as temp_docking_dir:

                temp_output_easydock_sdf = str(Path(temp_docking_dir) / f'{Path(fname).stem}_redocked.sdf')
                checkpoint_manager.state[model][docking_data.name]['docking'] = CheckpointStatus.RUNNING.value
                checkpoint_manager.save_state()
                temp_output_easydock_db = run_model_docking_easydock(fname=fname,
                                                                    model=model,
                                                                    docking_data=docking_data,
                                                                    temp_docking_dir=temp_docking_dir,
                                                                    cmd=generate_shared_cmd_docking_easydock(\
                                                                        docking_data=docking_data,
                                                                        input_data=input_data))

                post_processed_docking_fname = process_easydock_docking_output(temp_output_easydock_sdf)
                subprocess.run(['cp', post_processed_docking_fname, output_easydock_sdf], text=True, check=True)
                if Path(output_easydock_sdf).is_file():
                    checkpoint_manager.state[model][docking_data.name]['docking'] = CheckpointStatus.COMPLETED.value
                else:
                    checkpoint_manager.state[model][docking_data.name]['docking'] = CheckpointStatus.FAILED.value
                checkpoint_manager.save_state()


                if docking_data.plif_value:
                    checkpoint_manager.state[model][docking_data.name]['plif'] = CheckpointStatus.RUNNING.value
                    checkpoint_manager.save_state()
                    plif_result = run_plif_easydock(environment=docking_data.environment_value,
                                    db_fname=temp_output_easydock_db,
                                    protein_fname=input_data.pdb_value,
                                    reference_fname=input_data.sdf_value,
                                    reference_plif=docking_data.plif_value)

                    processed_temp_output = analyse_plif_easydock(plif_result=plif_result,
                                                                db_fname=temp_output_easydock_db,
                                                                similarity=docking_data.plif_similarity_value,
                                                                environment=docking_data.environment_value)
                    if processed_temp_output:
                        subprocess.run(['cp', processed_temp_output, output_easydock_plif],
                                    text=True, check=True)
                    checkpoint_manager.state[model][docking_data.name]['plif'] = CheckpointStatus.COMPLETED.value
                    checkpoint_manager.save_state()
        except (subprocess.CalledProcessError, PermissionError, KeyboardInterrupt) as exc:
            if checkpoint_manager.state[model][docking_data.name]['plif'] == CheckpointStatus.RUNNING.value:
                checkpoint_manager.state[model][docking_data.name]['plif'] = CheckpointStatus.FAILED.value
            elif checkpoint_manager.state[model][docking_data.name]['docking'] == CheckpointStatus.RUNNING.value:
                checkpoint_manager.state[model][docking_data.name]['docking'] = CheckpointStatus.FAILED.value
            checkpoint_manager.save_state()

            raise exc


def run_protonation_ligprep(schrodinger_dir: str, output_dir: str,
                            valid_molecule_file_dict: dict[str, str]) -> str:
    """run ligand protonation through schrodinger's ligprep. Temporary files are used because
    schrodinger wants to have all files into one folder. So, tempdir is used for smooth execution,
    although performance time are sacrificed, as these files are copied to another tempfiles by 
    schrodinger again.

    Args:
        schrodinger_dir (str): schrodigner root directory
        output_dir (str): directory for new protonated ligands
        valid_molecule_file_dict (dict[str, str]): dictionary containing the model name and filename
        for valid de novo generated ligand

    Returns:
        str: protonated ligands filename
    """

    cmd = [f'{schrodinger_dir}/ligprep']

    current_dir = Path.cwd()
    with tempfile.TemporaryDirectory() as temp_protonation_dir:
        os.chdir(temp_protonation_dir)
        for model, fname in valid_molecule_file_dict.items():
            model_cmd = deepcopy(cmd)

            ligprep_inp = str(Path(temp_protonation_dir) / f'{Path(fname).stem}_ligprep.inp')
            ligprep_output_dir = Path(Path(output_dir) / model)
            ligprep_output_dir.mkdir(parents=True, exist_ok=True)
            ligprep_output_fname = ligprep_output_dir / f'{Path(fname).stem}_protonated.sdf'

            if Path(ligprep_output_fname).is_file():
                logging.info(f'Ligprep protonation with {fname} has been done. Skipping...')
                continue

            subprocess.run(['cp', fname, temp_protonation_dir], text=True, check=True)
            prepare_ligprep_inp_file(ligand_fname=fname,
                                     output_fname=ligprep_output_fname,
                                     write_fname=ligprep_inp)

            protonation_cpu = str(fetch_ligprep_cpu(schrodinger_dir=schrodinger_dir))
            model_cmd += ['-inp', ligprep_inp,
                          '-NJOBS', protonation_cpu, '-JOBNAME', 
                          f'ligprep_{Path(temp_protonation_dir).name}', 
                          '-HOST', f'localhost:{protonation_cpu}', '-WAIT']

            logging.info(f'Ligprep is running with input: {fname}, output: {ligprep_output_fname}, cmd: {model_cmd}')
            subprocess.run(model_cmd, text=True, check=True)
            subprocess.run(['cp', Path(ligprep_output_fname).name, ligprep_output_fname], text=True, check=True)

    os.chdir(current_dir)

    return ligprep_output_fname


def prepare_hypothesis_fname_phase(schrodinger_env: str, rec_file: str,
                                   lig_file: str, center_coord:str) -> str:
    """automatically prepare hypothesis of protein-ligand. Note that it is very difficult to have
    a parameter to adjust number of matched sites, as it seems that schrodinger only allow user to
    give an appropriate matched sites based on the total pharmacophore sites. Hence, I will be using
    the default option to match all sites.

    Similarly, tempfiles are used because of schrodinger's requirement for everything inside one folder.

    Args:
        schrodinger_env (str): schrodinger root directory
        rec_file (str): pdb file of the protein
        lig_file (str): sdf file of the reference ligand, not de novo ligand
        center_coord (str): center of reference ligand separated by comma (e.g., "2.3,1,-4")

    Returns:
        str: phypo file of the hypothesis
    """

    current_tmpdir = Path.cwd()
    with tempfile.TemporaryDirectory() as tempdir:
        os.chdir(tempdir)
        rec_tempfile = str(Path(tempdir) / 'rec_file.pdb')
        lig_tempfile = str(Path(tempdir) / 'lig_file.sdf')

        subprocess.run(['cp', rec_file, rec_tempfile], check=True)
        subprocess.run(['cp', lig_file, lig_tempfile], check=True)

        jobname = f'{Path(rec_file).stem}_phypo'

        cmd = [f'{schrodinger_env}/utilities/epharmacophores', '-rec_file',
            rec_tempfile, '-lig_file', lig_tempfile,
            f'-site_center={center_coord}', '-in_place', '-j', jobname, 
            '-f', '7', '-site_dist', '2.0', '-pair_dist', '4.0', '-xvol', '-scale', 
            '0.5', '-buff', '2.0', '-limit', '5.0', '-HOST', 'localhost', '-WAIT']

        subprocess.run(cmd, check=True)

        phypo_file = f'{jobname}.phypo'

        output_fname = Path(rec_file).parent / phypo_file
        subprocess.run(['cp', phypo_file, output_fname], check=True)

    os.chdir(current_tmpdir)

    return str(output_fname)


def run_plif_phase(schrodinger_dir: str,
                input_data: config_input.InputConfig,
                library_fname: str,
                hypo_input: Optional[str]=None) -> str:
    """run phase screening acting as PLIF after glide. Docked ligands are required,
    because this function will not find other conformation (i.e., use inplace conformation).

    Similarly, tempfiles are used because of schrodinger's requirement for everything inside one folder.

    Args:
        schrodinger_dir (str): schrodinger root directory
        input_data (config_input.InputConfig): InputConfig
        library_fname (str): sdf filename for docked de novo ligands
        hypo_input (Optional[str], optional): hypothesis files for phase. Defaults to None.

    Raises:
        ValueError: if hypothesis file is not .zip or .phypo

    Returns:
        str: sdf filenames containing ligands that pass the pharmacophore filter
    """

    if not hypo_input or isinstance(hypo_input, bool):
        hypo_input = prepare_hypothesis_fname_phase(schrodinger_env=schrodinger_dir,
                                        rec_file=input_data.pdb_value,
                                        lig_file=input_data.sdf_value,
                                        center_coord=input_data.pocket_coord_value)
    elif not isinstance(hypo_input, str):
        raise ValueError(f'hypothesis for Phase screening expecting bool (generate own file), \
                         user phypo zip input, not {type(hypo_input)}')

    jobname = 'phase_screen_1'
    temp_output = f'{Path(jobname).stem}-hits.sdfgz'


    subprocess.run(['cp', hypo_input, 'hypo.phypo'], check=True)

    cmd = [f'{schrodinger_dir}/phase_screen', library_fname, 'hypo.phypo',
           jobname, '-inplace', '-keep', '999999999', '-report', '1',
           '-ex', '-HOST', 'localhost:10', '-TMPLAUNCHDIR', '-osd', '-WAIT']

    subprocess.run(cmd, check=True)

    return temp_output


def run_model_docking_glide(docking_data: config_redocking.GlideConfig,
                            tempdir: str, model: str, fname: str, grid_filepath: str) -> str:
    """run a single Glide SP redocking job for a single model. Temporary directory is used because
    schrodinger required all files to be in a folder, and current directory to be inside that folder
    as well.

    Args:
        docking_data (config_redocking.GlideConfig): GlideConfig
        tempdir (str): temporary directory
        model (str): model name used to generate the ligand
        fname (str): sdf filename of the valid generated ligand
        grid_filepath (str): grid.zip filepath prepared through schrodinger grid

    Returns:
        str: sdf file containing docked generated ligands.
    """

    if docking_data.protonation_value:
        run_protonation_ligprep(schrodinger_dir=docking_data.dir_value,
                                output_dir=tempdir,
                                valid_molecule_file_dict={model: fname})
        protonated_fname = Path(tempdir) / model / f'{Path(fname).stem}_protonated.sdf'
    else:
        protonated_fname = fname

    os.chdir(Path(protonated_fname).parent)
    jobname = f'glide_{Path(fname).stem}'
    inp_fname = Path(tempdir) / f'{jobname}.inp'
    prepare_glide_inp_file(grid_fname=grid_filepath,
                            ligand_fname=protonated_fname,
                            intra_hbonds=docking_data.reward_intra_hbonds_value,
                            write_fname=inp_fname)

    cmd = [f"{docking_data.dir_value}/glide", inp_fname, '-OVERWRITE',
        '-adjust', '-HOST', f'localhost:{fetch_glide_sp_cpu(docking_data.dir_value)}', '-TMPLAUNCHDIR', '-WAIT']

    temp_glide_sdf_output = f'{jobname}_lib.sdfgz'
    logging.info(f'Glide docking with cmd: {cmd=}')
    subprocess.run(cmd, text=True, check=True)

    return temp_glide_sdf_output


def unzip_sdf_gz(zipped_file: str, output_file: str):
    """unzip sdf.gz files within python. This is because there is some issue when
    trying to use gunzip, so I tried to use this way instead

    Args:
        zipped_file (str): zipped sdf.gz file
        output_file (str): unzipped .sdf output file
    """

    with gzip.open(zipped_file, 'rb') as f_in:
        with open(output_file, "wb") as f_out:
            shutil.copyfileobj(f_in, f_out)

def execute_docking_glide(docking_data: config_redocking.GlideConfig,
                          parameter_data: config_parameter.ParameterConfig,
                          input_data: config_input.InputConfig) -> None:
    """run Glide SP redocking and phase screening

    Args:
        docking_data (config_redocking.GlideConfig): GlideConfig
        parameter_data (config_parameter.ParameterConfig): ParameterConfig
        input_data (config_input.InputConfig): InputConfig
    """
    checkpoint_manager = initialise_redocking_checkpoint(docking_data=docking_data,
                                                         models=parameter_data.model_list_value)

    grid_filepath = GlideProtein(pdb_filepath=input_data.pdb_value,
                                native_ligand=list(Chem.SDMolSupplier(input_data.sdf_value, removeHs=False))[0],
                                grid_output_dirpath=str(Path(input_data.pdb_value).parent),
                                schrodinger_dirpath=docking_data.dir_value,
                                protein_preparation=docking_data.protein_preparation_value).grid_filepath

    current_dir = Path.cwd()
    for model, fname in fetch_model_file_from_model_dir(dir_path=docking_data.input_value,
                                                        model_list=parameter_data.model_list_value).items():

        output_glide_fname = Path(docking_data.output_value) / model / f'{Path(fname).stem}_glide.sdf'
        if docking_data.protonation_value:
            output_glide_fname = Path(docking_data.output_value) / model / f'{Path(fname).stem}_protonated_glide.sdf'

        Path(output_glide_fname).parent.mkdir(parents=True, exist_ok=True)

        if docking_data.plif_value:
            output_phase_plif = str(Path(docking_data.output_value) / model / 'plif' / f'{Path(fname).stem}_phase_plif.sdf')
            Path(output_phase_plif).parent.mkdir(parents=True, exist_ok=True)

            if checkpoint_manager.state[model][docking_data.name]['plif'] == CheckpointStatus.COMPLETED.value:
                logging.info(f'Phase Pharmacophore  with {output_phase_plif} has been done. Skipping...')
                continue

        if output_glide_fname.is_file() and not docking_data.plif_value:
            logging.info(f'Glide docking with {Path(fname).name} has been done. Skipping...')
            continue

        try:
            with tempfile.TemporaryDirectory() as tempdir:
                os.chdir(tempdir)

                if checkpoint_manager.state[model][docking_data.name]['docking'] != CheckpointStatus.COMPLETED.value:
                    checkpoint_manager.state[model][docking_data.name]['docking'] = CheckpointStatus.RUNNING.value
                    checkpoint_manager.save_state()
                    temp_glide_sdf_output = run_model_docking_glide(docking_data=docking_data,
                                                                tempdir=tempdir,
                                                                model=model,
                                                                fname=fname,
                                                                grid_filepath=grid_filepath)

                    unzip_sdf_gz(temp_glide_sdf_output, output_glide_fname)

                    if Path(output_glide_fname).is_file():
                        checkpoint_manager.state[model][docking_data.name]['docking'] = CheckpointStatus.COMPLETED.value
                    else:
                        checkpoint_manager.state[model][docking_data.name]['docking'] = CheckpointStatus.FAILED.value
                    checkpoint_manager.save_state()

                if docking_data.plif_value:
                    checkpoint_manager.state[model][docking_data.name]['plif'] = CheckpointStatus.RUNNING.value
                    checkpoint_manager.save_state()
                    temp_output = run_plif_phase(schrodinger_dir=docking_data.dir_value,
                                                input_data=input_data,
                                                library_fname=output_glide_fname,
                                                hypo_input=docking_data.plif_value)

                    unzip_sdf_gz(temp_output, output_phase_plif)

                    checkpoint_manager.state[model][docking_data.name]['plif'] = CheckpointStatus.COMPLETED.value
                    checkpoint_manager.save_state()

        except (subprocess.CalledProcessError, PermissionError, KeyboardInterrupt) as exc:
            if checkpoint_manager.state[model][docking_data.name]['plif'] == CheckpointStatus.RUNNING.value:
                checkpoint_manager.state[model][docking_data.name]['plif'] = CheckpointStatus.FAILED.value
            elif checkpoint_manager.state[model][docking_data.name]['docking'] == CheckpointStatus.RUNNING.value:
                checkpoint_manager.state[model][docking_data.name]['docking'] = CheckpointStatus.FAILED.value
            checkpoint_manager.save_state()

            raise exc

    os.chdir(current_dir)


def initialise_redocking_checkpoint(docking_data: config_redocking.GlideConfig | config_redocking.EasyDockConfig,
                                    models: list[str]) -> CheckpointManager:
    """initialize the checkpoint manager for genbench3d

    Args:
        docking_data (config_redocking.GlideConfig | config_redocking.EasyDockConfig): GlideConfig or EasyDockConfig
        models (list[str]): name of models used to do genbench3d

    Returns:
        CheckpointManager: class to track genbench3d progress for each model's parameter
    """

    checkpoint_manager = CheckpointManager()

    Path(docking_data.output_value).mkdir(parents=True, exist_ok=True)
    checkpoint_filepath = Path(docking_data.output_value) / 'redocking_checkpoint.json'

    if not checkpoint_filepath.is_file():
        checkpoint = defaultdict(dict)
        for model in models:
            checkpoint[model][docking_data.name] = {}
            checkpoint[model][docking_data.name]['docking'] = CheckpointStatus.PENDING.value
            if docking_data.plif_value:
                checkpoint[model][docking_data.name]['plif'] = CheckpointStatus.PENDING.value

        checkpoint_manager = checkpoint_manager.from_dict(checkpoint)
        checkpoint_manager.checkpoint_fname = checkpoint_filepath
        checkpoint_manager.save_state()
    else:
        checkpoint_manager = checkpoint_manager.from_json(checkpoint_filepath)

    return checkpoint_manager


def run_protonation(config_data: dict[str, dict]):
    """run protonation (main function of the p4_analyse_redocking.py)

    Args:
        config_data (dict[str, dict]): user's prepared YAML config as dictionary
    """
    redocking_data = config_redocking.RedockingConfig(config_data[config_constant.ANALYSIS_KEY])
    parameter_data = config_parameter.ParameterConfig(config_data=config_data)
    redocking_data.validate_config()
    protonation_data = redocking_data.protonation_value
    valid_molecule_file_dict = fetch_model_file_from_model_dir(dir_path=protonation_data.input_value,
                                                               model_list=parameter_data.model_list_value)

    if 'schrodinger' in str(Path(protonation_data.environment_value).name).lower():
        run_protonation_ligprep(schrodinger_dir=protonation_data.environment_value,
                                output_dir=protonation_data.output_value,
                                valid_molecule_file_dict=valid_molecule_file_dict)
    else:
        run_protonation_easydock(protonation_data=protonation_data,
                                 valid_molecule_file_dict=valid_molecule_file_dict)


def run_docking(config_data: dict[str, dict]):
    """run redocking (main function of the p4_analyse_redocking.py)

    Args:
        config_data (dict[str, dict]): user's prepared YAML config as dictionary
    """
    redocking_data = config_redocking.RedockingConfig(config_data[config_constant.ANALYSIS_KEY])
    parameter_data = config_parameter.ParameterConfig(config_data=config_data)
    input_data = config_input.InputConfig(config_data)
    docking_data_list = redocking_data.docking_value.valid_key_list


    for docking_data in docking_data_list:
        if isinstance(docking_data, config_redocking.EasyDockConfig):
            execute_docking_easydock(docking_data,
                                     parameter_data=parameter_data,
                                     input_data=input_data)
        elif isinstance(docking_data, config_redocking.GlideConfig):
            execute_docking_glide(docking_data, parameter_data=parameter_data, input_data=input_data)


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Redocking analysis after molecule preparation")

    parser.add_argument("-i", "--input", help="config yaml input file or txt file containing yaml filepath",
                        required=True, type=str)
    args = parser.parse_args()

    logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='%(asctime)s - %(module)s: - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')

    logging.info('Running p4_analyse_redocking.py')
    logging.info('Reading the config filename :)')


    yaml_file_list = read_input(args.input)

    for yaml_file in yaml_file_list:
        checkpoint_json = Path(yaml_file).parent
        with open(yaml_file, 'r', encoding='utf-8') as yaml_f:
            data_input = yaml.safe_load(yaml_f)

        analysis_dict = data_input[config_constant.ANALYSIS_KEY][config_constant.ANALYSIS_REDOCKING_KEY]
        if config_constant.ANALYSIS_REDOCKING_PROTONATION_KEY in analysis_dict:
            run_protonation(config_data=data_input)
        if config_constant.ANALYSIS_REDOCKING_DOCKING_KEY in analysis_dict:
            run_docking(config_data=data_input)
