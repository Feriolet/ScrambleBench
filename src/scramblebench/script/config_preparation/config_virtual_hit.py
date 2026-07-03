"""This file handles the virtual hit key of user's YAML config. 
Mainly used for p1_generate_config.py and other downstream path to
access the subkey field easily.
"""

import logging

from pathlib import Path
from typing import Any, Optional
from typing_extensions import Self

from scramblebench.script.config_preparation import config_constant
from scramblebench.script.config_preparation.config_utils import check_file_start_with_number

logger = logging.getLogger(__name__)

class VirtualHitConfig:
    """class containing information for virtual hit rate. This is utilized in p5_collect_data_analysis.py
    """
    def __init__(self, config_data: dict[str, Any]):
        """initialize class

        Args:
            config_data (dict[str, Any]): user's analysis dictionary with analysis virtual hit
                                          key as dictionary.
        """
        config_data = config_data[config_constant.ANALYSIS_VIRTUAL_HIT_KEY]
        self.input_name = 'input'
        self.output_name = 'output'
        self.query_name = 'query'
        self.filter_name = 'filter'

        self.input_value = config_data[self.input_name]
        self.output_value = config_data[self.output_name]
        self.query_value = config_data.get(self.query_name)
        self.filter_value = config_data.get(self.filter_name)

        if self.query_value:
            self.clean_query_line()


    def update(self, key: str, value: str) -> Self:
        """update new value of a key

        Args:
            key (str): valid VirtualHitConfig key
            value (str): new value that will be updated

        Returns:
            Self: VirtualHitConfig
        """
        return self


    def validate_config(self) -> None:
        """check if VirtualHitConfig is valid. No checking is done for self.query_value because
        it is impossible to check every single string and sanitize the input. Sanitizing the input
        will requiring interpreting the variable, the operator, and the value, which will demand a 
        lot of dependencies.

        Raises:
            TypeError: if input folder is not a string
            TypeError: if output folder is not a string
            ValueError: if chemical filter is not supported
        """
        check_file_start_with_number(self.input_value)
        check_file_start_with_number(self.output_value)

        if not isinstance(self.input_value, str):
            raise TypeError(f'Please put input folder as string, not {type(self.input_value)}')
        if not isinstance(self.output_value, str):
            raise TypeError(f'Please put output folder as string, not {type(self.output_value)}')

        SUPPORTED_CHEMICAL_FILTER = ['PAINS']

        if self.filter_value and self.filter_value not in SUPPORTED_CHEMICAL_FILTER:
            raise ValueError(f'chemical filter {self.filter_value} is not supported')


    def clean_query_line(self):
        """replace the self.query_value for AND and OR, because df.query() only accept the lowercase.
        """
        if 'AND' in self.query_value:
            self.query_value = self.query_value.replace('AND', 'and')
        if 'OR' in self.query_value:
            self.query_value = self.query_value.replace('OR', 'or')


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
        """"Write the standardised config format for GlideConfig.

        Args:
            prefix_dir (Optional[str], optional): the new directory that will replace/append the previous directory. 
            Defaults to None.

        Returns:
            dict[str, dict]: GlideConfig as a dictionary
        """
        if prefix_dir:
            self.update_input_output(prefix_dir=prefix_dir)

        return {config_constant.ANALYSIS_VIRTUAL_HIT_KEY: {self.input_name: str(self.input_value),
                                                          self.output_name: str(self.output_value),
                                                          self.query_name: self.query_value,
                                                          self.filter_name: self.filter_value}}
