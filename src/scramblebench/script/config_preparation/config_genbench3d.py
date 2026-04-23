from scramblebench.script.config_preparation import config_constant

from pathlib import Path
from typing import Any
import subprocess

class GenBench3DConfig:
    def __init__(self, config_data: dict[str, Any]):
        config_data = config_data[config_constant.GENBENCH3D_KEY]
        self.input_name = 'input'
        self.output_name = 'output'
        self.genbench_dir_name = 'genbench_dir'
        self.conda_env_name = 'conda_env'
        self.schrodinger_dir_name = 'schrodinger_dir'
        self.genbench_config_name = 'genbench_config'
        self.do_complex_forcefield_minimisation_name = 'do_complex_forcefield_minimisation'
        self.do_docking_forcefield_minimisation_name = 'do_docking_forcefield_minimisation'
        self.do_cancel_protonation_by_genbench_name = 'do_cancel_protonation_by_obabel_or_adfr'

        self.input_value = config_data[self.input_name]
        self.output_value = config_data[self.output_name]
        self.genbench_dir_value = config_data[self.genbench_dir_name]
        self.conda_env_value = config_data[self.conda_env_name]
        self.schrodinger_dir_value = config_data[self.schrodinger_dir_name]
        self.genbench_config_value = config_data[self.genbench_config_name]
        self.do_complex_forcefield_minimisation_value = config_data[self.do_complex_forcefield_minimisation_name]
        self.do_docking_forcefield_minimisation_value = config_data[self.do_docking_forcefield_minimisation_name]
        self.do_cancel_protonation_by_genbench_value = config_data[self.do_cancel_protonation_by_genbench_name]

    def update(self, key: str, value: str):
        if key == self.input_name:
            self.input_value = value
        elif key == self.output_name:
            self.output_value = value
        elif key == self.genbench_dir_name:
            self.genbench_dir_name = value
        elif key == self.conda_env_name:
            self.conda_env_name = value
        elif key == self.schrodinger_dir_name:
            self.schrodinger_dir_name = value
        elif key == self.genbench_config_name:
            self.genbench_config_name = value
        elif key == self.do_complex_forcefield_minimisation_name:
            self.do_complex_forcefield_minimisation_name = value
        elif key == self.do_docking_forcefield_minimisation_name:
            self.do_docking_forcefield_minimisation_name = value
        elif key == self.do_cancel_protonation_by_genbench_name:
            self.do_cancel_protonation_by_genbench_name = value        
        else:
            raise TypeError(f'no key called {key}')

        return self

    def validate_config(self):
        check_file_start_with_number(self.input_value)
        check_file_start_with_number(self.output_value)
        check_conda_env(self.conda_env_value)

        if not isinstance(self.input_value, str):
            raise TypeError(f'Please put input folder as string, not {type(self.input_value)}')
        if not isinstance(self.output_value, str):
            raise TypeError(f'Please put output folder as string, not {type(self.output_value)}')
        if not Path(self.genbench_dir_value).is_dir():
            raise ValueError(f'{self.genbench_dir_value} is not a directory. Please recheck your genbench input')
        if not Path(self.schrodinger_dir_value).is_dir():
            raise ValueError(f'{self.schrodinger_dir_value} is not a directory. Please recheck your schrodinger input')
        if not Path(self.genbench_config_value).is_file():
            raise ValueError(f'{self.genbench_config_value} is not a file. Please recheck your config input')
        if not isinstance(self.do_complex_forcefield_minimisation_value, bool):
            raise ValueError(f'The parameter do_complex_forcefield_minimisation is not a boolean. Please type True or False')
        if not isinstance(self.do_docking_forcefield_minimisation_value, bool) :
            raise ValueError(f'The parameter do_docking_forcefield_minimisation is not a boolean. Please type True or False')
        if not isinstance(self.do_cancel_protonation_by_genbench_value, bool):
            raise ValueError(f'The parameter do_cancel_protonation_by_genbench is not a boolean. Please type True or False')
        
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
            print('Please note that your conda is explicitly left out')
        else:
            raise ValueError(f'Conda environment {conda_env} does not exist')
        
def check_file_start_with_number(fname: str):
    try:
        float(Path(fname).name[0])
        raise TypeError(f'please do not use number as starting filename. change {fname}')
    except ValueError:
        pass