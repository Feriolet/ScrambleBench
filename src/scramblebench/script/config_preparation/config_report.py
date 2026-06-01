from scramblebench.script.config_preparation import config_constant, config_input, config_model, config_generation
from pathlib import Path
from typing import Any
import logging
from copy import deepcopy


logger = logging.getLogger(__name__)

class ReportConfig:

    def __init__(self, config_data: dict[str, Any]=None):
        self.name = config_constant.REPORT_KEY

        if self.name not in config_data:
            raise KeyError(f'You forgot to put "{self.name}" key in your config file!')
        
        report_data = config_data[self.name]

        self.rmsd_name = 'rmsd'
        self.docking_score_name = 'docking_score'
        self.plot_name = 'plot'
        self.qed_name = 'qed'

        self.rmsd_value = report_data.get(self.rmsd_name) or False
        self.docking_score_value = report_data.get(self.docking_score_name) or False
        self.plot_value = report_data.get(self.plot_name) or 'violin'
        self.qed_value = report_data.get(self.qed_name) or False
    
    
    def validate_config(self):

        SUPPORTED_PLOT_TYPE = ['violin', 'box', 'raincloud']

        if not isinstance(self.rmsd_value, bool):
            raise TypeError(f'Please put rmsd value as bool, not {type(self.input_value)}')
        if not isinstance(self.docking_score_value, bool):
            raise TypeError(f'Please put docking score as bool, not {type(self.input_value)}')
        if not isinstance(self.qed_value, bool):
            raise TypeError(f'Please put qed score as bool, not {type(self.input_value)}')
        
        if self.plot_value.lower() not in SUPPORTED_PLOT_TYPE:
            raise ValueError(f'ScrambleBench do not support {self.plot_value}')

    def write(self):

        data = {self.plot_name: self.plot_value}

        if self.rmsd_value:
            data[self.rmsd_name] = self.rmsd_value
        if self.docking_score_value:
            data[self.docking_score_name] = self.docking_score_value
        if self.qed_value:
            data[self.qed_name] = self.qed_value

        return {self.name: data}
