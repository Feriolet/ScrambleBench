"""This file handles the analysis key of user's YAML config. 
Mainly used for p1_generate_config.py and other downstream path to
access the subkey field easily.
"""

import logging

from typing import Any, Optional

from scramblebench.script.config_preparation import config_constant
from scramblebench.script.config_preparation.config_genbench3d import GenBench3DConfig
from scramblebench.script.config_preparation.config_redocking import RedockingConfig
from scramblebench.script.config_preparation.config_diversity import DiversityConfig
from scramblebench.script.config_preparation.config_virtual_hit import VirtualHitConfig


logger = logging.getLogger(__name__)

class AnalysisConfig:
    """class containing subclass of analysis. Subclasses are created as the type of analyses
    depends on the user. So, I just use some dependency injection.
    """


    def __init__(self, config_data: dict[str, Any]):
        """initialize class

        Args:
            config_data (dict[str, Any]): user's prepared YAML config with analysis key as dictionary
        """

        self.name = config_constant.ANALYSIS_KEY
        analysis_data = config_data[self.name]

        self.valid_key_list = []
        if config_constant.ANALYSIS_GENBENCH3D_KEY in analysis_data.keys():
            self.valid_key_list.append(GenBench3DConfig(analysis_data))
        if config_constant.ANALYSIS_REDOCKING_KEY in analysis_data.keys():
            self.valid_key_list.append(RedockingConfig(analysis_data))
        if config_constant.ANALYSIS_DIVERSITY_KEY in analysis_data.keys():
            self.valid_key_list.append(DiversityConfig(analysis_data))
        if config_constant.ANALYSIS_VIRTUAL_HIT_KEY in analysis_data.keys():
            self.valid_key_list.append(VirtualHitConfig(analysis_data))


    def update(self, key: str, value: Any):
        """update new value of a key. I have not made this yet, since I am too lazy zzzzzz

        Args:
            key (str): _description_
            value (Any): _description_
        """
        pass


    def validate_config(self):
        """check if AnalysisConfig is valid by looping each analysis subclass
        """
        for config in self.valid_key_list:
            config.validate_config()


    def write(self, prefix_dir: Optional[str]=None) -> dict[str, dict]:
        """Write the standardised config format for AnalysisConfig by looping analysis subclass
        Args:
            prefix_dir (str, optional): the new directory that will replace/append the previous directory. 
            Defaults to None.

        Returns:
            dict[str, dict]:: AnalysisConfig as dictionary
        """

        data = {}
        for config in self.valid_key_list:
            data = data | config.write(prefix_dir)

        return {self.name  : data}
