from typing import Any, Optional
from pathlib import Path

from scramblebench.script.utils.error_handler import FileDataError, FileTypeError, DirNotFoundError
from scramblebench.script.config_preparation import config_constant

import Bio
import rdkit
from Bio.PDB import PDBParser
from rdkit import Chem
import os
from oddt.toolkits.extras.rdkit import fixer
import numpy as np
import logging
import re
import subprocess

logger = logging.getLogger(__name__)


def check_file_start_with_number(fname: str):
    try:
        float(Path(fname).name[0])
        raise TypeError(f'please do not use number as starting filename. change {fname}')
    except ValueError:
        pass

def check_conda_env(conda_env: str) -> None:

    # Run 'conda env list' command and capture the output
    output = subprocess.run(
        "conda env list | awk '{ print $1}'",
        shell=True,
        capture_output=True,
        text=True,
        check=True
    ).stdout.split('\n')

    if not conda_env in output:
        if conda_env == 'not_applicable' or conda_env is None:
            logging.warning('Please note that some of your conda is explicitly left out')
        else:
            raise ValueError(f'Conda environment {conda_env} does not exist')
        

class RedockingConfig:

    def __init__(self, config_data: dict[str, Any]):

        self.name = config_constant.ANALYSIS_REDOCKING_KEY
        redocking_data = config_data[self.name]

        self.protonation_value = None
        if config_constant.ANALYSIS_REDOCKING_PROTONATION_KEY in redocking_data.keys():
            self.protonation_value = RedockingProtonationConfig(redocking_data)

        self.docking_value = DockingConfig(redocking_data)

    def update(self, key, value):
        pass

    def validate_config(self):

        if self.protonation_value:
            self.protonation_value.validate_config()
        self.docking_value.validate_config()

    def write(self, prefix_dir=None):

        data = {}
        if self.protonation_value:
            data = data | self.protonation_value.write(prefix_dir)
        data = data | self.docking_value.write(prefix_dir)

        return {self.name   : data}


class RedockingProtonationConfig:

    def __init__(self, config_data: dict[str, Any]):

        self.name = config_constant.ANALYSIS_REDOCKING_PROTONATION_KEY
        protonation_data = config_data[self.name]

        self.input_name = 'input'
        self.output_name = 'output'
        self.environment_name = 'env'
        self.method_name = 'method'

        self.input_value = protonation_data[self.input_name]
        self.output_value = protonation_data[self.output_name]
        self.environment_value = protonation_data[self.environment_name]
        self.method_value = protonation_data[self.method_name]

    def update(self, key, value):
        pass

    def validate_config(self):

        EASYDOCK_PROTONATION_PROGRAM = ['molgpka', 'unipka', 'chemaxon']
        SCHRODINGER_PROTONATION_PROGRAM = ['ligprep']

        check_file_start_with_number(self.input_value)
        check_file_start_with_number(self.output_value)

        if not isinstance(self.input_value, str):
            raise TypeError(f'Please put input folder as string, not {type(self.input_value)}')
        if not isinstance(self.output_value, str):
            raise TypeError(f'Please put output folder as string, not {type(self.output_value)}')

        if not isinstance(self.method_value, str):
            raise TypeError(f'Please put protonation method as string, not {type(self.input_value)}')
        
        if self.method_value.lower() in EASYDOCK_PROTONATION_PROGRAM:
            check_conda_env(self.environment_value)
        elif self.method_value.lower() in SCHRODINGER_PROTONATION_PROGRAM:
            schrodinger_dirpath = Path(self.environment_value).resolve()

            if not schrodinger_dirpath.is_dir():
                raise DirNotFoundError(f'schrodinger directory is not found in {str(schrodinger_dirpath)}')
            if 'schrodinger' not in str(schrodinger_dirpath):
                raise ValueError(f'We did not found the schrodinger directory in {str(schrodinger_dirpath)}')

            self.environment_value = Path(self.environment_value).resolve()
        else:
            raise ValueError(f'Your protonation programme {self.method_value} is not supported by our program. Currently, we only support {",".join(EASYDOCK_PROTONATION_PROGRAM + SCHRODINGER_PROTONATION_PROGRAM)}')
    
    def update_input_output(self, prefix_dir):
        if '/' in self.input_value[0]:
            logging.warning('prefix dir for generation is provide. Ignoring absolute path provided by input key and only take the last dir')
            self.input_value = Path(prefix_dir) / Path(self.input_value).name
        else:
            self.input_value = Path(prefix_dir) / self.input_value
        
        if '/' in self.output_value[0]:
            logging.warning('prefix dir for generation is provide. Ignoring absolute path provided by output key and only take the last dir')
            self.output_value = Path(prefix_dir) / Path(self.output_value).name
        else:
            self.output_value = Path(prefix_dir) / self.output_value  

    def write(self, prefix_dir=None):
        if prefix_dir:
            self.update_input_output(prefix_dir)

        return {self.name   : {self.input_name: str(self.input_value),
                               self.output_name: str(self.output_value),
                               self.environment_name: str(self.environment_value),
                               self.method_name: str(self.method_value)}}

class DockingConfig:

    def __init__(self, config_data: dict[str, Any]):

        self.name = config_constant.ANALYSIS_REDOCKING_DOCKING_KEY

        if self.name not in config_data.keys():
            raise ValueError(f'{self.name} must be in the {config_constant.ANALYSIS_REDOCKING_KEY} key, because you want to redock')
        
        docking_data = config_data[self.name]

        self.valid_key_list = []
        if config_constant.ANALYSIS_DOCKING_EASYDOCK_KEY in docking_data.keys():
            self.valid_key_list.append(EasyDockConfig(docking_data))
        if config_constant.ANALYSIS_DOCKING_SCHRODINGER_KEY in docking_data.keys():
            self.valid_key_list.append(GlideConfig(docking_data))

    def update(self, key, value):
        pass

    def validate_config(self):
        for config in self.valid_key_list:
            config.validate_config()

    def write(self, prefix_dir=None):
        data = {}
        for config in self.valid_key_list:
            data = data | config.write(prefix_dir)
            
        return {self.name : data}


class EasyDockConfig:

    def __init__(self, config_data: dict[str, Any]):

        self.name = config_constant.ANALYSIS_DOCKING_EASYDOCK_KEY
        easydock_data = config_data[self.name]

        self.input_name = 'input'
        self.output_name = 'output'
        self.environment_name = 'conda_env'
        self.docking_name = 'docking_program'
        self.protonation_name = 'protonation'
        self.config_name = 'config_fname'
        self.cpu_name = 'ncpu'
        self.protein_pdbqt_preparation_name = 'protein_pdbqt_preparation'
        self.protein_preparation_name = 'protein_preparation'
        self.protein_pdbqt_executable_name = 'protein_pdbqt_executable'
        self.plif_name = 'ref_plif'
        self.plif_similarity_name = 'plif_similarity_threshold'

        self.input_value = easydock_data[self.input_name]
        self.output_value = easydock_data[self.output_name]
        self.environment_value = easydock_data[self.environment_name]
        self.docking_value = easydock_data.get(self.docking_name) or 'vina'
        self.protonation_value = easydock_data.get(self.protonation_name) or None
        self.config_value = easydock_data.get(self.config_name) 
        self.cpu_value = easydock_data.get(self.cpu_name)
        self.protein_preparation_value = easydock_data.get(self.protein_preparation_name) or 'pdbfixer'
        self.protein_pdbqt_preparation_value = easydock_data.get(self.protein_pdbqt_preparation_name) or 'obabel'
        self.protein_pdbqt_executable_value = easydock_data.get(self.protein_pdbqt_executable_name) or 'obabel'
        self.plif_value = easydock_data.get(self.plif_name)
        self.plif_similarity_value = easydock_data.get(self.plif_similarity_name)

        if isinstance(self.plif_value, str):
            self.plif_value = self.plif_value.split(',')
        if isinstance(self.plif_similarity_value, str):
            self.plif_similarity_value = float(self.plif_similarity_value)

        if self.config_value is None:
            if self.docking_value == 'vina':
                self.config_value = str(Path(__file__).parent.parent / 'utils/docking_utils/easydock_config.yml') 
            else:
                raise ValueError(f'You are using {self.docking_value} program to dock. Please fill in the config fname')
    
    def update(self, key, value):
        pass

    def validate_config(self):

        EASYDOCK_PROTONATION_PROGRAM = ['molgpka', 'unipka', 'chemaxon', None]
        EASYDOCK_DOCKING_PROGRAM = ['vina', 'gnina', 'smina', 'vina-gpu', 'qvina', 'server']

        SUPPORTED_PDBQT_PREPARATION_PROGRAM = ['adfr', 'obabel']
        SUPPORTED_PROTEIN_PREPARATION_PROGRAM = ['pdbfixer', 'obabel']

        check_file_start_with_number(self.input_value)
        check_file_start_with_number(self.output_value)

        if not isinstance(self.input_value, str):
            raise TypeError(f'Please put input folder as string, not {type(self.input_value)}')
        if not isinstance(self.output_value, str):
            raise TypeError(f'Please put output folder as string, not {type(self.output_value)}')
        if not isinstance(self.config_value, str):
            raise TypeError(f'Please put config fname as string, not {type(self.config_value)}')

        config_filepath = Path(self.config_value)
        if not config_filepath.is_file():
            raise FileNotFoundError(f'{str(config_filepath)} is not found')
        if config_filepath.suffix not in ['.yaml', '.yml']:
            raise ValueError(f'{str(config_filepath)} should have a .yaml or .yml suffix, not {config_filepath.suffix}')

        if self.protein_pdbqt_preparation_value:
            if self.protein_pdbqt_preparation_value not in SUPPORTED_PDBQT_PREPARATION_PROGRAM:
                raise ValueError(f'{self.protein_pdbqt_preparation_value} is not supported by scramblebench')
        
        if self.protein_pdbqt_executable_value:
            if not Path(self.protein_pdbqt_executable_value).is_file() and self.protein_pdbqt_executable_value != 'obabel':
                raise FileNotFoundError(f'{self.protein_pdbqt_executable_value} is not found')

        if self.protein_preparation_value:
            if isinstance(self.protein_preparation_value, bool):
                self.protein_preparation_value = 'pdbfixer'

            elif self.protein_preparation_value not in SUPPORTED_PROTEIN_PREPARATION_PROGRAM:
                if not (Path(self.protein_preparation_value).is_dir() and 'schrodinger' in Path(self.protein_preparation_value).name):
                    raise ValueError(f'{self.protein_preparation_value} is not supported by scramblebench')

        if self.protonation_value:
            if not isinstance(self.protonation_value, str):
                raise TypeError(f'Please put protonation method as string, not {type(self.input_value)}')
            if self.protonation_value.lower() not in EASYDOCK_PROTONATION_PROGRAM and 'schrodinger' not in self.protonation_value:
                raise ValueError(f'{self.protonation_value} is not supported by easydock')

        if self.plif_value:
            if isinstance(self.plif_value, bool):
                pass
            elif isinstance(self.plif_value, list):
                for plif in self.plif_value:
                    if len(plif.split('.')) != 3:
                        raise ValueError(f'plif value should follow the residue.chain.interaction format, not {plif}')
            else:
                raise ValueError(f'unsupported type {type(self.plif_value)}')


        if self.plif_similarity_value and not isinstance(self.plif_similarity_value, (int, float)):
            raise ValueError(f'please write plif similarity value as float')

        if not isinstance(self.docking_value, str):
            raise TypeError(f'Please put docking method as string, not {type(self.input_value)}')
        if self.docking_value.lower() not in EASYDOCK_DOCKING_PROGRAM:
            raise ValueError(f'{self.docking_value} is not supported by easydock')
        
        if self.cpu_value:
            self.cpu_value = int(self.cpu_value)

        check_conda_env(self.environment_value)

    def update_input_output(self, prefix_dir):
        if '/' in self.input_value[0]:
            logging.warning('prefix dir for generation is provide. Ignoring absolute path provided by input key and only take the last dir')
            self.input_value = Path(prefix_dir) / Path(self.input_value).name
        else:
            self.input_value = Path(prefix_dir) / self.input_value
        
        if '/' in self.output_value[0]:
            logging.warning('prefix dir for generation is provide. Ignoring absolute path provided by output key and only take the last dir')
            self.output_value = Path(prefix_dir) / Path(self.output_value).name
        else:
            self.output_value = Path(prefix_dir) / self.output_value  

    def write(self, prefix_dir=None):

        DEFAULT_PLIF_SIMILARITY_THRESHOLD = 0.8
        if prefix_dir:
            self.update_input_output(prefix_dir)
        
        data =  {self.input_name: str(self.input_value),
                               self.output_name: str(self.output_value),
                               self.environment_name: str(self.environment_value),
                               self.docking_name: str(self.docking_value),
                               self.config_name: str(Path(self.config_value).resolve())}

        if self.protonation_value:
            if 'schrodinger' in self.protonation_value:
                self.check_schrodinger_ligprep()
            
        data[self.protonation_name] = self.protonation_value
            
        if not self.protein_pdbqt_executable_value:
            if not self.protein_pdbqt_preparation_value:
                self.protein_pdbqt_preparation_value = 'adfr'
                self.protein_pdbqt_executable_value = str(Path(__file__).parent.parent / 'utils/prepare_receptor')
            else:
                self.protein_pdbqt_executable_value = str(Path(__file__).parent.parent / 'utils/prepare_receptor')
        
        if self.protein_pdbqt_preparation_value == 'obabel' or self.protein_preparation_value == 'obabel':
            self.protein_pdbqt_preparation_value = 'obabel'
            self.protein_preparation_value = 'obabel'

        if self.cpu_value:
            data[self.cpu_name] = self.cpu_value
        
        if self.plif_value:
            if isinstance(self.plif_value, list):
                data[self.plif_name] = ','.join(self.plif_value)
            elif isinstance(self.plif_value, bool):
                data[self.plif_name] = True

            if not self.plif_similarity_value:
                self.plif_similarity_value = DEFAULT_PLIF_SIMILARITY_THRESHOLD
            
            data[self.plif_similarity_name] = self.plif_similarity_value

        data[self.protein_pdbqt_preparation_name] = self.protein_pdbqt_preparation_value
        data[self.protein_preparation_name] = self.protein_preparation_value
        data[self.protein_pdbqt_executable_name] = self.protein_pdbqt_executable_value

        return {self.name : data}
    

    def check_schrodinger_ligprep(self):
        if 'ligprep' in self.protonation_value:
            assert 'schrodinger' in Path(self.protonation_value).parent.name
            self.protonation_value = str(Path(self.protonation_value).parent.resolve())
            print(f'{self.protonation_value}')
        else:
            assert 'schrodinger' in Path(self.protonation_value).name
            self.protonation_value = str(Path(self.protonation_value).resolve())

class GlideConfig:

    def __init__(self, config_data: dict[str, Any]):

        self.name = config_constant.ANALYSIS_DOCKING_SCHRODINGER_KEY
        glide_data = config_data[self.name]
    
        self.input_name = 'input'
        self.output_name = 'output'
        self.dir_name = 'schrodinger_dir'
        self.reward_intra_hbonds_name = 'reward_intra_hbonds'
        self.protonation_name = 'protonation'
        self.protein_preparation_name = 'protein_preparation'
        self.plif_name = 'plif_input'

        self.input_value = glide_data[self.input_name]
        self.output_value = glide_data[self.output_name]
        self.dir_value = glide_data[self.dir_name]
        self.reward_intra_hbonds_value = glide_data.get(self.reward_intra_hbonds_name) or False
        self.protonation_value = glide_data.get(self.protonation_name) or None
        self.protein_preparation_value = glide_data.get(self.protein_preparation_name) or None
        self.plif_value = glide_data.get(self.plif_name)

    def update(self, key, value):
        pass

    def validate_config(self):

        SCHRODINGER_PROTONATION_PROGRAM = ['ligprep', True, False, None]
        SCHRODINGER_PROTEIN_PREPARATION_PROGRAM = ['protwizard', True, False, None]

        check_file_start_with_number(self.input_value)
        check_file_start_with_number(self.output_value)

        if not isinstance(self.input_value, str):
            raise TypeError(f'Please put input folder as string, not {type(self.input_value)}')
        if not isinstance(self.output_value, str):
            raise TypeError(f'Please put output folder as string, not {type(self.output_value)}')

        schrodinger_dirpath = Path(self.dir_value)
        if not schrodinger_dirpath.is_dir():
            raise DirNotFoundError(f'schrodinger directory is not found in {str(schrodinger_dirpath)}')
        if 'schrodinger' not in str(schrodinger_dirpath):
            raise ValueError(f'We did not found the schrodinger directory in {str(schrodinger_dirpath)}')
        
        if self.protonation_value:
            if not isinstance(self.protonation_value, (bool, str)):
                raise TypeError(f'Please put protonation method as string or boolean, not {type(self.protonation_value)}')
            if self.protonation_value.lower() not in SCHRODINGER_PROTONATION_PROGRAM:
                raise ValueError(f'{self.protonation_value} is not supported by easydock')

        if self.protein_preparation_value:
            if not isinstance(self.protein_preparation_value, (bool, str)):
                raise TypeError(f'Please put protonation method as string or boolean, not {type(self.protein_preparation_value)}')
            if self.protein_preparation_value.lower() not in SCHRODINGER_PROTEIN_PREPARATION_PROGRAM:
                raise ValueError(f'{self.protein_preparation_value} is not supported by easydock')
            
        if self.reward_intra_hbonds_value is not None:
            assert isinstance(self.reward_intra_hbonds_value, bool)

        if self.plif_value:
            if isinstance(self.plif_value, str):
                assert Path(self.plif_value).is_file()
            elif not isinstance(self.plif_value, bool):
                raise ValueError(f'plif must be either string, NoneType, or bool, not {type(self.plif_value)}')
            

    def update_input_output(self, prefix_dir):
        if '/' in self.input_value[0]:
            logging.warning('prefix dir for generation is provide. Ignoring absolute path provided by input key and only take the last dir')
            self.input_value = Path(prefix_dir) / Path(self.input_value).name
        else:
            self.input_value = Path(prefix_dir) / self.input_value
        
        if '/' in self.output_value[0]:
            logging.warning('prefix dir for generation is provide. Ignoring absolute path provided by output key and only take the last dir')
            self.output_value = Path(prefix_dir) / Path(self.output_value).name
        else:
            self.output_value = Path(prefix_dir) / self.output_value  

    def write(self, prefix_dir=None):

        if prefix_dir:
            self.update_input_output(prefix_dir)

        data = {self.input_name: str(self.input_value),
                self.output_name: str(self.output_value),
                self.dir_name: str(Path(self.dir_value).resolve()),
                self.protonation_name: self.protonation_value,
                self.protein_preparation_name: self.protein_preparation_value}
        
        if self.reward_intra_hbonds_value:
            data[self.reward_intra_hbonds_name] = self.reward_intra_hbonds_value

        if self.plif_value:
            data[self.plif_name] = self.plif_value
            
        return {self.name   : data}

