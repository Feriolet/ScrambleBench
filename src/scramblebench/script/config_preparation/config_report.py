"""This file handles the report key of user's YAML config. 
Mainly used for p1_generate_config.py and other downstream path to
access the subkey field easily.
"""

import logging

from typing import Any
from scramblebench.script.config_preparation import config_constant


logger = logging.getLogger(__name__)

class ReportConfig:
    """class containing information regarding what to put in report
    """

    def __init__(self, config_data: dict[str, Any]=None):
        """initialize class

        Args:
            config_data (dict[str, Any], optional): user's prepared YAML config with report
                                                    key as dictionary. Defaults to None.

        Raises:
            KeyError: if report key is not in the config dictionary
        """
        self.name = config_constant.REPORT_KEY

        if self.name not in config_data:
            raise KeyError(f'You forgot to put "{self.name}" key in your config file!')

        report_data = config_data[self.name]

        self.rmsd_name = 'rmsd'
        self.docking_score_name = 'docking_score'
        self.plot_name = 'plot'
        self.qed_name = 'qed'
        self.validity3d_name = 'validity3d'
        self.virtual_hit_name = 'virtual_hit'

        self.rmsd_value = report_data.get(self.rmsd_name) or False
        self.docking_score_value = report_data.get(self.docking_score_name) or False
        self.plot_value = report_data.get(self.plot_name) or 'violin'
        self.qed_value = report_data.get(self.qed_name) or False
        self.validity3d_value = report_data.get(self.validity3d_name) or False
        self.virtual_hit_value = report_data.get(self.virtual_hit_name) or False

    def validate_config(self) -> None:
        """check if ReportConfig is valid

        Raises:
            TypeError: if rmsd subfield is not a bool or NoneType
            TypeError: if docking score subfield is not a bool or NoneType
            TypeError: if qed value subfield is not a bool or NoneType
            TypeError: if validity3d subfield is not a bool or NoneType
            TypeError: if virtual hit rate subfield is not a bool or NoneType
            ValueError: if plot type subfield is not supported
        """

        SUPPORTED_PLOT_TYPE = ['violin', 'box', 'raincloud']

        if not isinstance(self.rmsd_value, bool):
            raise TypeError(f'Please put rmsd value as bool, not {type(self.rmsd_value)}')
        if not isinstance(self.docking_score_value, bool):
            raise TypeError(f'Please put docking score as bool, not {type(self.docking_score_value)}')
        if not isinstance(self.qed_value, bool):
            raise TypeError(f'Please put qed score as bool, not {type(self.qed_value)}')
        if not isinstance(self.validity3d_value, bool):
            raise TypeError(f'Please put validity3d metric as bool, not {type(self.validity3d_value)}')
        if not isinstance(self.virtual_hit_value, bool):
            raise TypeError(f'Please put virtual hit metric as bool, not {type(self.validity3d_value)}')

        if self.plot_value.lower() not in SUPPORTED_PLOT_TYPE:
            raise ValueError(f'ScrambleBench do not support {self.plot_value}')

    def write(self) -> dict[str, dict]:
        """write a standardized config format for ReportConfig

        Returns:
            dict[str, dict]: ReportConfig as a dictionary
        """

        data = {self.plot_name: self.plot_value}

        if self.rmsd_value:
            data[self.rmsd_name] = self.rmsd_value
        if self.docking_score_value:
            data[self.docking_score_name] = self.docking_score_value
        if self.qed_value:
            data[self.qed_name] = self.qed_value
        if self.validity3d_value:
            data[self.validity3d_name] = self.validity3d_value
        if self.virtual_hit_value:
            data[self.virtual_hit_name] = self.virtual_hit_value

        return {self.name: data}
