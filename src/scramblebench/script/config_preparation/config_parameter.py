from scramblebench.script.config_preparation import config_constant, config_input, config_model, config_generation
from pathlib import Path
from typing import Any
import logging
from copy import deepcopy


logger = logging.getLogger(__name__)

class ParameterConfig:

    def __init__(self, config_data: dict[str, Any]=None):
        self.name = config_constant.PARAMETER_KEY


        self.protein_name = 'protein_name'
        self.model_list_name = 'model_list'
        self.num_sample_name = 'num_sample'
        self.batch_parameter_name = 'batch_parameter'

        self.protein_value = None
        self.model_list_value = []
        self.batch_parameter_dict = {}

        if config_data:

            if self.name not in config_data:
                raise KeyError(f'You forgot to put "{self.name}" key in your config file!')
            
            parameter_data = config_data[self.name]
            self.protein_value = parameter_data[self.protein_name]
            self.model_list_value = parameter_data[self.model_list_name]
            self.batch_parameter_dict = parameter_data[self.batch_parameter_name]
            self.num_sample_value = parameter_data[self.num_sample_name]

    def create(self, config_data: dict[str, Any], batch_parameter):

        self.name = config_constant.PARAMETER_KEY
        input_data = config_input.InputStructure(input_dict=config_data[config_constant.INPUT_KEY])
        self.protein_name = 'protein_name'
        self.protein_value = input_data.protein_value

        model_data = config_model.ModelConfig(config_data=config_data)
        self.model_list_name = 'model_list'
        self.model_list_value = model_data.get_model_list()

        self.num_sample_value = config_generation.GenerationConfig(config_data=config_data)\
                                                            .parameter_value.num_sample_value
        
        self.batch_parameter_name = 'batch_parameter'
        self.parameter_dict = batch_parameter
    
        return self
    
    
    def write(self):
        return {self.name: {self.protein_name: self.protein_value,
                            self.model_list_name: self.model_list_value,
                            self.num_sample_name: self.num_sample_value,
                            self.batch_parameter_name: self.parameter_dict}}
    
    def deep_get(self, dictionary: dict, nested_key: list):
        copied_dict = deepcopy(dictionary)
        for key in nested_key:
            copied_dict = copied_dict.get(key)

            if copied_dict is None:
                logging.exception(f'The dictionary {dictionary} has no value for key {key}')
                raise KeyError(f'Key {key} not found for dictionary {dictionary}')
            if not isinstance(copied_dict, dict) and key != nested_key[-1]:
                logging.warning(f'Your dictionary {dictionary} overshoot the nested key. The value for {key} is not a dictionary. Ignoring subsequent keys')
                return copied_dict
        
        return copied_dict
    