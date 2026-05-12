from typing import Any, Callable
from pathlib import Path

import yaml
import argparse
import logging
import sys
import subprocess
import os

from scramblebench.script.config_preparation import config_constant, config_genbench3d, config_input, config_redocking, config_parameter
from scramblebench.script.utils.error_handler import DirNotFoundError
from scramblebench.script.utils.process_data import read_input, fetch_model_folder_name, find_file_name_through_regex

from copy import deepcopy
import json
from collections import defaultdict
from enum import Enum
from multiprocessing import Pool
import tempfile

from openmm.app import PDBFile
from pdbfixer import PDBFixer
from openbabel import pybel


def fetch_valid_prepared_molecule_file(config_class, parameter_class: config_parameter.ParameterConfig) -> dict[str, str]:
    model_for_generation_list = parameter_class.model_list_value

    generation_folder_dirpath = Path(config_class.input_value)
    if not generation_folder_dirpath.is_dir():
        logging.exception('You have prepared your molecule yet!')
        raise DirNotFoundError(f'{generation_folder_dirpath} is not found! Please make sure you run p3_prepare_molecule.py')
    
    valid_molecule_file_dict = {}
    for model in model_for_generation_list:

        matched_fname_list = find_file_name_through_regex(character=model, file_format='.sdf', dirname=Path(generation_folder_dirpath) / model)
        if len(matched_fname_list) > 1:
            logging.exception(f'We found more than 1 matching file for {model} model: {matched_fname_list}. Please ensure only 1 is detected')
            raise ValueError(f'We found more than 1 matching file for {model} model: {matched_fname_list}. Please ensure only 1 is detected')
        elif len(matched_fname_list) == 0:
            logging.warning(f'There are no matched file for {model} model in {generation_folder_dirpath}. Make sure this is intended')

        valid_molecule_file_dict[model] = str(matched_fname_list[0])
    
    return valid_molecule_file_dict


def run_ligprep_protonation(protonation_data, valid_molecule_file_dict):

    cmd = [f'{protonation_data.environment_value}/ligprep']
    if not Path(protonation_data.output_value).is_dir():
        Path(protonation_data.output_value).mkdir(parents=True, exist_ok=True)

    for fname in valid_molecule_file_dict.values():

        with tempfile.TemporaryDirectory() as tempfile_dir:      

            model_cmd = deepcopy(cmd)

            ligprep_inp = str(Path(tempfile_dir) / f'{Path(fname).stem}_ligprep.inp')
            ligprep_output_fname = str(f'{Path(protonation_data.output_value) / Path(fname).stem}_protonated.sdf')
            ligprep_inp_data = f'''
            INPUT_FILE_NAME   {fname}
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
            OUT_SD  {ligprep_output_fname}

            '''

            with open(ligprep_inp, 'w') as ligprep_fname:
                ligprep_fname.write(ligprep_inp_data)
            
            model_cmd += ['-inp', ligprep_inp, '-NJOBS', '1', '-JOBNAME', f'ligprep_{Path(tempfile_dir).name}', '-HOST', 'localhost:1', '-WAIT']
            subprocess.run(model_cmd, text=True)

def run_easydock_protonation(protonation_data, valid_molecule_file_dict):

    cmd = ['conda', 'run', '-n', protonation_data.environment_value,
            'easydock']

    if not Path(protonation_data.output_value).is_dir():
        Path(protonation_data.output_value).mkdir(parents=True, exist_ok=True)

    for fname in valid_molecule_file_dict.values():

        with tempfile.TemporaryDirectory() as tempfile_dir:      
            output_easydock = str(Path(tempfile_dir) / f'{Path(fname).stem}_protonated.db')
            model_cmd = deepcopy(cmd)                               
            model_cmd += ['-i', fname, '-o', output_easydock, '--protonation', protonation_data.method_value]

            print(f'{model_cmd=}')
            subprocess.run(model_cmd, text=True)

            FETCH_SDF_SCRIPT = Path(__file__).parent / 'utils' / 'easydock_fetch_sdf_from_db.py'
            protonated_cmd = ['conda', 'run', '-n', protonation_data.environment_value,
                            'python', FETCH_SDF_SCRIPT, '-i', output_easydock, '-o', str(f'{Path(protonation_data.output_value) / Path(fname).stem}_protonated.sdf')]

            subprocess.run(protonated_cmd, text=True)

def run_protonation(config_data):
    redocking_data = config_redocking.RedockingConfig(config_data[config_constant.ANALYSIS_KEY])
    parameter_data = config_parameter.ParameterConfig(config_data=config_data)
    redocking_data.validate_config()
    protonation_data = redocking_data.protonation_value
    valid_molecule_file_dict = fetch_valid_prepared_molecule_file(config_class=protonation_data,
                                                                  parameter_class=parameter_data)

    if 'schrodinger' in str(protonation_data.environment_value).lower():
        run_ligprep_protonation(protonation_data=protonation_data,
                                valid_molecule_file_dict=valid_molecule_file_dict)
    else:
        run_easydock_protonation(protonation_data=protonation_data,
                                 valid_molecule_file_dict=valid_molecule_file_dict)


class VinaProtein():
    
    def __init__(self, 
                 pdb_filepath: str,
                 prepare_receptor_bin_path: str) -> None:
        self.protein_filepath = pdb_filepath
        self.prepare_receptor_bin_path = prepare_receptor_bin_path
        self._pdbqt_filepath = pdb_filepath.replace('.pdb', 
                                                   '.pdbqt')

    @property
    def pdbqt_filepath(self):
        if not os.path.exists(self._pdbqt_filepath):
            self.vina_prepare_receptor(output_pdbqt_filepath=self._pdbqt_filepath) # Using default configuration
        assert os.path.exists(self._pdbqt_filepath), \
            'Something went wrong during Vina receptor preparation'
        return self._pdbqt_filepath
        
    def clean_protein(self,
                      input_pdb_filepath: str,
                      output_pdb_filepath: str,
                      pH: float = 7.4,
                      ) -> None:
        """
        Adds hydrogens to given pH for a pdb protein (input_filepath)
        """
        logging.info(f'Cleaning protein from {input_pdb_filepath} to {output_pdb_filepath}')
        
        fixer = PDBFixer(filename=input_pdb_filepath)
        # We only use PDBFixer to reinitialize the chains and segment id combinations
        # (e.g. chain A and B, having both segment A and B, will be renamed chain A to D)
        # fixer.findMissingResidues() # it cannot find them with CrossDocked because it does not contain sequence info
        # fixer.findNonstandardResidues()
        # fixer.replaceNonstandardResidues()
        # fixer.findMissingAtoms()
        # fixer.addMissingAtoms()
        # fixer.addMissingHydrogens(pH=pH)
        # fixer.removeHeterogens(keepWater=False)
        intermediate_filepath = output_pdb_filepath.replace('.pdb', '_no_h.pdb')
        PDBFile.writeFile(fixer.topology, fixer.positions, open(intermediate_filepath, 'w'))
            
        # pdbfixer might miss some hydrogens because of missing atoms
        molecule = list(pybel.readfile("pdb", str(intermediate_filepath)))[0]
        # # import pdb;pdb.set_trace()
        # molecule.OBMol.CorrectForPH(pH) # the correct functions fails for some groups: -NH3+ becomes -NH4+
        # molecule.removeh() # some pdb files have hydrogens, these might mess up the next step
        molecule.addh() 
        molecule.write("pdb", str(output_pdb_filepath), overwrite=True)  

    def vina_prepare_receptor(self,
                              output_pdbqt_filepath: str,
                              ligand_name: str = None,
                                chain: str = None,
                                preparation_method: str = 'adfr',
                                pH: float = 7.4
                                ) -> None:
        """
        inspired from teachopencadd talktorial 15 on protein_ligand_docking
        """
        
        if preparation_method == 'adfr':
            self.adfr_receptor_preparation(input_pdb_filepath=self.protein_filepath,
                                           output_pdbqt_filepath=output_pdbqt_filepath)
        else:
            self.pdb_to_pdbqt()

           
    def adfr_receptor_preparation(self,
                                  input_pdb_filepath: str,
                                  output_pdbqt_filepath: str,
                                  ) -> None:
        """
        input_pdb_filepath must be a pbd file that only contains the protein with
        hydrogens
        """
        logging.info(f'Preparing protein from {input_pdb_filepath} to {output_pdbqt_filepath}')
        arg_list = [self.prepare_receptor_bin_path,
                    f'-r {input_pdb_filepath}',
                    f'-o {output_pdbqt_filepath}']
        cmd = ' '.join(arg_list)
        os.system(cmd)
        
        
    def pdb_to_pdbqt(self,
                     ph: float = 7.4,
                     ) -> None:
        molecule = list(pybel.readfile("pdb", str(self._protein_filepath)))[0]
        self.ob_mol_to_pdbqt(molecule, ph)


    def ob_mol_to_pdbqt(self,
                        molecule,
                        ph: float = 7.4,
                        ) -> None:
        
        # add hydrogens at given pH
        molecule.OBMol.CorrectForPH(ph)
        molecule.addh()
        # add partial charges to each atom
        for atom in molecule.atoms:
            atom.OBAtom.GetPartialCharge()

        molecule.write("pdbqt", str(self._pdbqt_filepath), overwrite=True)
        
        # Only keep ATOM and TER lines in pdbqt file
        with open(self.pdbqt_filepath, 'r') as f:
            lines = [line.strip() for line in f.readlines()]
        new_lines = [line 
                     for line in lines 
                     if line.startswith('ATOM') or line.startswith('TER')]
        with open(self.pdbqt_filepath, 'w') as f:
            for line in new_lines:
                f.write(line)
                f.write('\n')


def prepare_easydock_config(config_fname, input_data: config_input.InputStructure, output_dir):
    protein_pdb = input_data.pdb_value
    pocket_coordinate = input_data.pocket_coord_value.split(',')

    with open(config_fname, 'r') as config_f:
        config_data = yaml.load(config_f, Loader=yaml.SafeLoader)
    

    grid_txt = f'''
    center_x = {pocket_coordinate[0]}
    center_y = {pocket_coordinate[1]}
    center_z = {pocket_coordinate[2]}
    size_x = 25.0
    size_y = 25.0
    size_z = 25.0
    '''
    grid_fname = str(Path(protein_pdb).parent / f'{input_data.protein_name}_grid.txt')
    with open(grid_fname, 'w') as grid_f:
        grid_f.write(grid_txt)
    
    ADFR_PATH = Path(__file__).parent / 'utils/ADFRsuite_x86_64Linux1.0/clean_install/bin/prepare_receptor'

    config_data['protein'] = VinaProtein(pdb_filepath=protein_pdb, prepare_receptor_bin_path=ADFR_PATH).pdbqt_filepath
    config_data['protein_setup'] = grid_fname
    if config_data.get('exhaustiveness') is None:
        config_data['exhaustiveness'] = 32
    if config_data.get('n_poses') is None:
        config_data['n_poses'] = 5
    if config_data.get('ncpu') is None:
        config_data['ncpu'] = max(1, min(len(os.sched_getaffinity(0)) - 5, 5))
    if config_data.get('seed') is None:
        config_data['seed'] = 42

    new_config_fname = Path(output_dir) / config_fname
    with open(new_config_fname, 'w') as config_f:
        yaml.dump(config_data, config_f)

    return new_config_fname


def run_easydock_docking(docking_data: config_redocking.EasyDockConfig, parameter_data, input_data):
    
    cmd = ['conda', 'run', '-n', docking_data.environment_value,
            'easydock',  '--sdf']

    if not Path(docking_data.output_value).is_dir():
        Path(docking_data.output_value).mkdir(parents=True, exist_ok=True)

    config_fname = prepare_easydock_config(docking_data.config_value, input_data, output_dir=docking_data.output_value)
    cmd += ['--config', config_fname]

    valid_molecule_file_dict = fetch_valid_prepared_molecule_file(config_class=docking_data,
                                                                  parameter_class=parameter_data)
    if docking_data.protonation_value:
        cmd += ['--protonation', docking_data.protonation_value]
    
    if docking_data.docking_value:
        cmd += ['--program', docking_data.docking_value]

    for _, fname in valid_molecule_file_dict.items():

        with tempfile.TemporaryDirectory() as tempfile_dir:      
            temp_output_easydock_db = str(Path(tempfile_dir) / f'{Path(fname).stem}_redocked.db')
            temp_output_easydock_sdf = str(Path(tempfile_dir) / f'{Path(fname).stem}_redocked.sdf')
            output_easydock_sdf =  str(Path(docking_data.output_value) / f'{Path(fname).stem}_easydock_redocked.sdf')
            model_cmd = deepcopy(cmd)                               
            model_cmd += ['-i', fname, '-o', temp_output_easydock_db, ]


            print(f'{model_cmd=}')
            subprocess.run(model_cmd, text=True)
            subprocess.run(['cp', temp_output_easydock_sdf, output_easydock_sdf], text=True)


def run_docking(config_data):
    redocking_data = config_redocking.RedockingConfig(config_data[config_constant.ANALYSIS_KEY])
    parameter_data = config_parameter.ParameterConfig(config_data=config_data)
    input_data = config_input.InputStructure(config_data[config_constant.INPUT_KEY])
    docking_data_list = redocking_data.docking_value.valid_key_list

    for docking_data in docking_data_list:
        if isinstance(docking_data, config_redocking.EasyDockConfig):
            run_easydock_docking(docking_data, parameter_data=parameter_data, input_data=input_data)
        elif isinstance(docking_data, config_redocking.GlideConfig):
            print('booo')


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Prepare generated molecules for downstream analysis")

    parser.add_argument("-i", "--input", help="config yaml input file or txt file containing yaml filepath", required=True, type=str)
    args = parser.parse_args()

    logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='%(asctime)s - %(module)s: - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
    
    logging.info('Running p4_analyse_redocking.py')
    logging.info('Reading the config filename :\)')


    yaml_file_list = read_input(args.input)

    for yaml_file in yaml_file_list:
        checkpoint_json = Path(yaml_file).parent
        data_input = yaml.safe_load(open(yaml_file, 'r'))

        if config_constant.ANALYSIS_REDOCKING_PROTONATION_KEY in data_input[config_constant.ANALYSIS_KEY][config_constant.ANALYSIS_REDOCKING_KEY]:
            run_protonation(config_data=data_input)
        if config_constant.ANALYSIS_REDOCKING_DOCKING_KEY in data_input[config_constant.ANALYSIS_KEY][config_constant.ANALYSIS_REDOCKING_KEY]:
            run_docking(config_data=data_input)