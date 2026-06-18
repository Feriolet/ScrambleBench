from scramblebench.script.config_preparation import config_constant, config_input, config_model, config_generation
from pathlib import Path
from typing import Any
import logging
from copy import deepcopy


logger = logging.getLogger(__name__)

class VirtualHitConfig:
    def __init__(self, config_data: dict[str, Any]):
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

    def update(self, key: str, value: str):
        return self

    def validate_config(self):
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
        if 'AND' in self.query_value:
            self.query_value = self.query_value.replace('AND', 'and')
        if 'OR' in self.query_value:
            self.query_value = self.query_value.replace('OR', 'or')


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
            self.update_input_output(prefix_dir=prefix_dir)
            
        return {config_constant.ANALYSIS_VIRTUAL_HIT_KEY: {self.input_name: str(self.input_value),
                                                          self.output_name: str(self.output_value),
                                                          self.query_name: self.query_value,
                                                          self.filter_name: self.filter_value}}

def check_file_start_with_number(fname: str):
    try:
        float(Path(fname).name[0])
        raise TypeError(f'please do not use number as starting filename. change {fname}')
    except ValueError:
        pass
