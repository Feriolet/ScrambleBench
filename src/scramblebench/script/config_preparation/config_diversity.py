"""This file handles the diversity analysis key of user's YAML config. 
Mainly used for p1_generate_config.py and other downstream path to
access the subkey field easily.
"""

import logging

from pathlib import Path
from typing import Any, Optional
from typing_extensions import Self

from scramblebench.script.config_preparation import config_constant
from scramblebench.script.config_preparation.config_utils import check_conda_env, check_file_start_with_number

logger = logging.getLogger(__name__)


class DiversityConfig:
    """class containing information regarding diversity analysis
    """


    def __init__(self, config_data: dict[str, Any]):
        """initialize classes

        Args:
            config_data (dict[str, Any]): user's prepared YAML config with analysis diversity
                                          key as dictionary
        """
        config_data = config_data[config_constant.ANALYSIS_DIVERSITY_KEY]
        self.input_name = 'input'
        self.output_name = 'output'
        self.environment_name = 'conda_env'

        self.method_name = 'method'
        self.method_distance_name = 'distance'
        self.method_diversity_name = 'diversity'

        self.input_value = config_data[self.input_name]
        self.output_value = config_data[self.output_name]
        self.method_value = config_data.get(self.method_name)
        self.environment_value = config_data[self.environment_name]

        if config_data.get(self.method_name) is None:
            self.method_value = [{self.method_distance_name: 'ecfp',
                                  self.method_diversity_name: 'hamdiv'}]


    def update(self, key: str, value: str) -> Self:
        """update new value of a key. I have not made this yet, since I am too lazy zzzzzz

        Args:
            key (str): _description_
            value (str): _description_

        Returns:
            Self: _description_
        """
        return self


    def validate_config(self):
        """check whether DiversityConfig is valid

        Raises:
            TypeError: if input folder starts with a number 
            TypeError: if output folder starts with a number
            ValueError: if diversity method is not a list
            ValueError: if the list of diversity method is not a dict
        """

        check_file_start_with_number(self.input_value)
        check_file_start_with_number(self.output_value)
        check_conda_env(self.environment_value)

        if not isinstance(self.input_value, str):
            raise TypeError(f'Please put input folder as string, not {type(self.input_value)}')
        if not isinstance(self.output_value, str):
            raise TypeError(f'Please put output folder as string, not {type(self.output_value)}')


        if not isinstance(self.method_value, list):
            raise ValueError(f'Your value in the diversity.method section is not a list of dictionary, \
                             but rather {type(self.method_value)}')
        for method in self.method_value:
            if not isinstance(method, dict):
                raise ValueError(f'Your value in the diversity.method section is not a list of dictionary,\
                                but rather list of {type(method)}')

            self.validate_diversity_distance_metric_pair(diversity=method[self.method_diversity_name],
                                                         distance=method[self.method_distance_name])


    def validate_diversity_distance_metric_pair(self, diversity: str, distance: str) -> None:
        """validate whether ScrambleBench supports pairs of diversity and distance. Mostly inspired from the
        HamDiv paper

        Args:
            diversity (str): diversity metric - relating to how to manipulate distance 
            distance (str): distance metric - relating to how different is a pair of molecule

        Raises:
            ValueError: diversity metric is not a string
            ValueError: distance metric is not a string
            ValueError: wrong diversity and distance metric
        """

        SUPPORTED_DIVERSITY_DISTANCE_PAIR = [('hamdiv', 'ecfp'), ('hamdiv', 'mces'), ('average', 'ecfp'),
                                             ('richness',), ('rs',), ('fg',),('bm',), ('generic_bm',), ('intdiv',),
                                             ('sumdiv',), ('diam',), ('sumdiam',), ('sumbot',), ('bot',), ('dpp',) ]
        if not isinstance(diversity, str):
            raise ValueError(f'You diversity value must be a string, not {type(diversity)}')

        if distance is None:
            pair = (diversity.lower(),)
        else:
            if not isinstance(distance, str):
                raise ValueError(f'You distance value must be a string, not {type(diversity)}')
            pair = (diversity.lower(), distance.lower())

        if pair not in SUPPORTED_DIVERSITY_DISTANCE_PAIR:
            raise ValueError(f'The pair {diversity} and {distance} is not supported')


    def update_input_output(self, prefix_dir: str):
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
        """"Write the standardised config format for DiversityConfig.

        Args:
            prefix_dir (Optional[str], optional): the new directory that will replace/append the previous directory. 
            Defaults to None.

        Returns:
            dict[str, dict]: DiversityConfig as a dictionary
        """
        if prefix_dir:
            self.update_input_output(prefix_dir=prefix_dir)

        return {config_constant.ANALYSIS_DIVERSITY_KEY: {self.input_name: str(self.input_value),
                                                          self.output_name: str(self.output_value),
                                                          self.environment_name: self.environment_value,
                                                          self.method_name: self.method_value}}
