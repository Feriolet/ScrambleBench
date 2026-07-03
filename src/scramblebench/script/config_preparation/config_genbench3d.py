"""protonate and dock de novo ligand."""

import logging

from pathlib import Path
from typing import Any, Optional
from typing_extensions import Self

from scramblebench.script.config_preparation import config_constant
from scramblebench.script.config_preparation.config_utils import check_conda_env, check_file_start_with_number


logger = logging.getLogger(__name__)

class GenBench3DConfig:
    """class containing information regarding genbench3d
    """
    def __init__(self, config_data: dict[str, Any]):
        """initialize class. There are a lot of attributes here, but as genbench3d has a lot of optimization,
        I have decided to just dump everything here. Maybe it is better if I can read the genbench3d.yaml file, but
        to reduce the attributes. May or may not be better.

        Args:
            config_data (dict[str, Any]): user's prepared YAML config with analysis genbench3d
                                          key as dictionary
        """

        config_data = config_data[config_constant.ANALYSIS_GENBENCH3D_KEY]
        self.input_name = 'input'
        self.output_name = 'output'
        self.genbench3d_dir_name = 'genbench3d_dir'
        self.conda_env_name = 'conda_env'
        self.schrodinger_dir_name = 'schrodinger_dir'
        self.genbench3d_config_name = 'genbench3d_config'
        self.do_complex_forcefield_minimisation_name = 'do_complex_forcefield_minimisation'
        self.do_docking_forcefield_minimisation_name = 'do_docking_forcefield_minimisation'
        self.do_cancel_protonation_by_genbench3d_name = 'skip_genbench3d_protonation'

        self.input_value = config_data[self.input_name]
        self.output_value = config_data[self.output_name]
        self.genbench3d_dir_value = config_data[self.genbench3d_dir_name]
        self.conda_env_value = config_data[self.conda_env_name]
        self.genbench3d_config_value = config_data[self.genbench3d_config_name]
        self.schrodinger_dir_value = config_data.get(self.schrodinger_dir_name)
        self.do_complex_forcefield_minimisation_value = config_data.get(self.do_complex_forcefield_minimisation_name) or False
        self.do_docking_forcefield_minimisation_value = config_data.get(self.do_docking_forcefield_minimisation_name) or False
        self.do_cancel_protonation_by_genbench3d_value = config_data.get(self.do_cancel_protonation_by_genbench3d_name) or False


    def update(self, key: str, value: Any) -> Self:
        """update value of a key

        Args:
            key (str): valid GenBench3DConfig key
            value (Any): new value that will be updated

        Raises:
            TypeError: invalid GenBench3DConfig key

        Returns:
            Self: GenBench3DConfig
        """
        if key == self.input_name:
            self.input_value = value
        elif key == self.output_name:
            self.output_value = value
        elif key == self.genbench3d_dir_name:
            self.genbench3d_dir_name = value
        elif key == self.conda_env_name:
            self.conda_env_name = value
        elif key == self.schrodinger_dir_name:
            self.schrodinger_dir_name = value
        elif key == self.genbench3d_config_name:
            self.genbench3d_config_name = value
        elif key == self.do_complex_forcefield_minimisation_name:
            self.do_complex_forcefield_minimisation_name = value
        elif key == self.do_docking_forcefield_minimisation_name:
            self.do_docking_forcefield_minimisation_name = value
        elif key == self.do_cancel_protonation_by_genbench3d_name:
            self.do_cancel_protonation_by_genbench3d_name = value
        else:
            raise TypeError(f'no key called {key}')

        return self


    def validate_config(self):
        """check if GenBench3DConfig is valid

        Raises:
            TypeError: if input folder is not a string
            TypeError: if output folder is not a string
            ValueError: if genbench3d folder does not exist
            ValueError: if schrodinger root directory does not exist
            ValueError: if genbench3d config fname does not exist
            ValueError: if MMFF94 minimisation is not a boolean
            ValueError: if docking mininplace is not a boolean
            ValueError: if cancel genbench3d protonation is not a boolean
        """

        check_file_start_with_number(self.input_value)
        check_file_start_with_number(self.output_value)
        check_conda_env(self.conda_env_value)

        if not isinstance(self.input_value, str):
            raise TypeError(f'Please put input folder as string, not {type(self.input_value)}')
        if not isinstance(self.output_value, str):
            raise TypeError(f'Please put output folder as string, not {type(self.output_value)}')

        if not Path(self.genbench3d_dir_value).is_dir():
            raise ValueError(f'{self.genbench3d_dir_value} is not a directory. \
                             Please recheck your genbench input')

        if self.schrodinger_dir_value:
            if not (Path(self.schrodinger_dir_value).is_dir() and \
                    'schrodinger' in self.schrodinger_dir_value):
                raise ValueError(f'{self.schrodinger_dir_value} is not a directory. \
                                 Please recheck your schrodinger input')

        if not Path(self.genbench3d_config_value).is_file():
            raise ValueError(f'{self.genbench3d_config_value} is not a file. \
                             Please recheck your config input')

        if not isinstance(self.do_complex_forcefield_minimisation_value, bool):
            raise ValueError('The parameter do_complex_forcefield_minimisation is not a boolean. \
                             Please type True or False')
        if not isinstance(self.do_docking_forcefield_minimisation_value, bool) :\
            raise ValueError('The parameter do_docking_forcefield_minimisation is not a boolean. \
                             Please type True or False')
        if not isinstance(self.do_cancel_protonation_by_genbench3d_value, bool):
            raise ValueError('The parameter do_cancel_protonation_by_genbench is not a boolean. \
                             Please type True or False')


    def update_input_output(self, prefix_dir: str) -> None:
        """update the input and output folder for post generation. This is utilized in p1_generate_config.py,
        because user may provide multiple batch parameters. Hence, the output folder directory path depends on these
        parameters.

        Args:
            prefix_dir (str): the new directory that will replace/append the previous directory
        """

        if '/' in self.input_value[0]:
            logging.warning('prefix dir for generation is provided. \
                            Ignoring absolute path provided by input key and only take the last dir')
            self.input_value = Path(prefix_dir) / Path(self.input_value).name
        else:
            self.input_value = Path(prefix_dir) / self.input_value

        if '/' in self.output_value[0]:
            logging.warning('prefix dir for generation is provided. \
                            Ignoring absolute path provided by output key and only take the last dir')
            self.output_value = Path(prefix_dir) / Path(self.output_value).name
        else:
            self.output_value = Path(prefix_dir) / self.output_value


    def write(self, prefix_dir: Optional[str]=None) -> dict[str, dict]:
        """"Write the standardised config format for GenBench3DConfig.

        Args:
            prefix_dir (Optional[str], optional): the new directory that will replace/append the previous directory. 
            Defaults to None.

        Returns:
            dict[str, dict]: GenBench3DConfig as a dictionary
        """
        if prefix_dir:
            self.update_input_output(prefix_dir=prefix_dir)

        data =  {self.input_name: str(self.input_value),
                self.output_name: str(self.output_value),
                self.conda_env_name: self.conda_env_value,
                self.genbench3d_dir_name: str(Path(self.genbench3d_dir_value).resolve()),
                self.genbench3d_config_name: str(Path(self.genbench3d_config_value).resolve()),
                self.do_complex_forcefield_minimisation_name: self.do_complex_forcefield_minimisation_value,
                self.do_docking_forcefield_minimisation_name: self.do_docking_forcefield_minimisation_value,
                self.do_cancel_protonation_by_genbench3d_name: self.do_cancel_protonation_by_genbench3d_value}

        if self.schrodinger_dir_value:
            data[self.schrodinger_dir_name] = self.schrodinger_dir_value

        return {config_constant.ANALYSIS_GENBENCH3D_KEY: data}
