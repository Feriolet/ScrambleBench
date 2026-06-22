"""This file handles the generation key of user's YAML config. 
Mainly used for p1_generate_config.py and other downstream path to
access the subkey field easily.
"""

import logging

from pathlib import Path
from typing import Any
from typing_extensions import Self

from scramblebench.script.config_preparation import config_constant

logger = logging.getLogger(__name__)


class GenerationParameterConfig:
    """class for generation parameter used for downstream analysis. This class will save
    the parameter used to generate the de novo ligand such as box size, num sample, and 
    job title (deprecated). Not to be confused with ParameterConfig class
    """

    def __init__(self, parameter_data: dict[str, Any]):
        """initialize the class

        Args:
            parameter_data (dict[str, Any]): user's prepared YAML config with 
                                             generation parameter key as dictionary

        Raises:
            ValueError: if box size is not a list, int, float, or string of digits separated by comma
            ValueError: if num sample is not a list, int, float, or string of digits separated by comma
        """

        if config_constant.GENERATION_PARAMETER_KEY in parameter_data:
            parameter_data = parameter_data[config_constant.GENERATION_PARAMETER_KEY]
        else:
            parameter_data = {}

        self.box_size_name = 'box_size'
        self.num_sample_name = 'num_sample'
        self.title_name = 'name'

        if parameter_data.get(self.box_size_name) is None:
            logging.warning('Unspecified parameter box size. Setting default to 16 A')
            self.box_size_value = 16
        elif isinstance(parameter_data[self.box_size_name], (int, float)):
            self.box_size_value = parameter_data[self.box_size_name]
        elif isinstance(parameter_data[self.box_size_name], str):
            self.box_size_value = [int(num.strip()) for num in parameter_data[self.box_size_name].split(',')]
        elif isinstance(parameter_data[self.box_size_name], list):
            self.box_size_value = parameter_data[self.box_size_name]
        else:
            raise ValueError(f'unsupported type of {parameter_data[self.box_size_name]}: \
                             {type(parameter_data[self.box_size_name])}')

        if parameter_data.get(self.num_sample_name) is None:
            logging.warning('Unspecified parameter num sample. Setting default to 100')
            self.num_sample_value = 100
        elif isinstance(parameter_data[self.num_sample_name], int):
            self.num_sample_value = parameter_data[self.num_sample_name]
        elif isinstance(parameter_data[self.num_sample_name], str):
            self.num_sample_value = [int(num.strip()) for num in parameter_data[self.num_sample_name].split(',')]
        elif isinstance(parameter_data[self.num_sample_name], list):
            self.num_sample_value = parameter_data[self.num_sample_name]
        else:
            raise ValueError(f'unsupported type of {parameter_data[self.num_sample_name]}: \
                             {type(parameter_data[self.num_sample_name])}')

        self.title_value = parameter_data.get(self.title_name) or 'generic_title'


    def update(self, key: str, value: str) -> Self:
        """update value of key of the class

        Args:
            key (str): valid GenerationParameterConfig key
            value (str): new value that will be updated

        Raises:
            ValueError: if box size is not int, float or iterables that cannot be converted to numbers
            ValueError: if num sample is not int, float or iterables that cannot be converted to numbers
            TypeError: invalid GenerationParameterConfig key

        Returns:
            Self: GenerationParameterConfig
        """
        if key == self.box_size_name:
            try:
                self.box_size_value = [int(num.strip()) for num in str(value).split(',')]
            except ValueError as exc:
                raise ValueError(f'num sample should be integer or list of integer, not {value}') from exc

        elif key == self.num_sample_name:
            try:
                self.num_sample_value = [int(num.strip()) for num in str(value).split(',')]
            except ValueError as exc:
                raise ValueError(f'num sample should be integer or list of integer, not {value}') from exc
        elif key in self.title_name:
            self.title_value[key] = value

        else:
            raise TypeError(f'no key called {key}')

        return self


    def validate_config(self):
        """validate whether the user GenerationParameterConfig is valid.

        Raises:
            ValueError: if box size is not int, float or iterables that cannot be converted to numbers
            ValueError: if box size is not int, float or iterables that cannot be converted to numbers
            ValueError: if num sample  is not int, float or iterables that cannot be converted to numbers
            ValueError: if num sample is not int, float or iterables that cannot be converted to numbers
        """
        assert isinstance(self.box_size_value, (int, float))

        if isinstance(self.box_size_value, list):
            if not all(isinstance(num, (int, float)) for num in self.box_size_value):
                raise ValueError(f'{self.box_size_value} needs to be an integer or list of integers')
        elif isinstance(self.box_size_value, (int, float)):
            pass
        else:
            raise ValueError(f'{self.box_size_value} needs to be an integer or list of integers')


        if isinstance(self.num_sample_value, list):
            if not all(isinstance(num, int) for num in self.num_sample_value):
                raise ValueError(f'{self.num_sample_value} needs to be an integer or list of integers')
        elif isinstance(self.num_sample_value, int):
            pass
        else:
            raise ValueError(f'{self.num_sample_value} needs to be an integer or list of integers')

    def write(self) -> dict[str, dict]:
        """Write the standardised config format for GenerationParameterConfig

        Returns:
            dict[str, dict]: GenerationParameterConfig as dictionary
        """
        return {config_constant.GENERATION_PARAMETER_KEY: {self.box_size_name : self.box_size_value,
                                                           self.num_sample_name : self.num_sample_value,
                                                           self.title_name : self.title_value }}


class GenerationConfig:
    """"class containing generation information and GenerationParameterConfig
    """

    def __init__(self, config_data: dict[str, Any]):
        """initialize class

        Args:
            config_data (dict[str, Any]): user's YAML config with generation key as dictionary
        """
        config_data = config_data[config_constant.GENERATION_KEY]

        self.name = config_constant.GENERATION_KEY
        self.input_name = 'input'
        self.output_name = 'output'
        self.script_name = 'script_pathfile'

        self.input_value = config_data[self.input_name]
        self.output_value = config_data[self.output_name]
        self.parameter_value = GenerationParameterConfig(config_data)

        #to do: create tests for default or incorrect generation script
        self.script_value = config_data.get(self.script_name) or config_constant.GENERATION_TEMPLATE_SCRIPT_PATH

    def update(self, key: str, value: str) -> Self:
        """update value of key of the class

        Args:
            key (str): valid GenerationConfig key
            value (str): new value that will be updated

        Returns:
            Self: GenerationConfig
        """
        if key == self.input_name:
            self.input_value = value
        elif key == self.output_name:
            self.output_value = value
        else:
            self.parameter_value = self.parameter_value.update(key, value)

        return self


    def validate_config(self):
        """check whether the user GenerationConfig is valid.

        Raises:
            FileNotFoundError: if the bash script for generation does not exist
            ValueError: if generation script not in .sh file format
        """
        check_file_start_with_number(self.input_value)
        check_file_start_with_number(self.output_value)

        assert isinstance(self.input_value, str)
        assert isinstance(self.output_value, str)

        self.parameter_value.validate_config()

        if not Path(self.script_value).is_file():
            raise FileNotFoundError(f'{self.script_value} was not found. Please check your directory again')
        if Path(self.script_value).suffix != '.sh':
            raise ValueError(f'Script generation file must be in .sh format. Instead, we detect {self.script_value}')


    def update_output(self, prefix_dir: str):
        """this function will update the output directory by appending a prefix directory.
        This prefix directory is typically created by the p1_generate_config.py, where user
        may submit multiple parameters and so will generate a subfolder for each parameter.

        Args:
            prefix_dir (str): prefix directory to append to the self.output_value
        """

        if '/' in self.output_value[0]:
            logging.warning('prefix dir for generation is provided. \
                            Ignoring absolute path provided by output key and only take the last dir')
            self.output_value = Path(prefix_dir) / Path(self.output_value).name
        else:
            self.output_value = Path(prefix_dir) / self.output_value

    def write(self, prefix_dir: str=None) -> dict[str, dict]:
        """Write the standardised config format for GenerationConfig

        Args:
            prefix_dir (str, optional): prefix directory to append to self.output_value. Defaults to None.

        Returns:
            dict[str, dict]: GenerationConfig as dictionary
        """
        if prefix_dir:
            self.update_output(prefix_dir=prefix_dir)

        return {config_constant.GENERATION_KEY : {self.input_name: str(Path(self.input_value).resolve()),
                                                  self.output_name: str(Path(self.output_value).resolve()),
                                                  self.script_name: str(Path(self.script_value).resolve())} |
                                                  self.parameter_value.write()}

def check_file_start_with_number(fname: str):
    """check if filename starts with a number. Having numbers at start may cause an issue in processing

    Args:
        fname (str): filename

    Raises:
        TypeError: if filename starts with a number
    """
    try:
        float(Path(fname).name[0])
        raise TypeError(f'please do not use number as starting filename. change {fname}')
    except ValueError:
        pass
