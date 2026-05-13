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
        self.repeat_parameter_name = 'repeat_parameter'

        self.protein_value = None
        self.model_list_value = []
        self.repeat_parameter_dict = {}

        if config_data:
            parameter_data = config_data[self.name]
            self.protein_value = parameter_data[self.protein_name]
            self.model_list_value = parameter_data[self.model_list_name]
            self.repeat_parameter_dict = parameter_data[self.repeat_parameter_name]


    def create(self, config_data: dict[str, Any]):

        self.name = config_constant.PARAMETER_KEY
        input_data = config_input.InputConfig(config_data=config_data)
        self.protein_name = 'protein_name'
        self.protein_value = list(input_data.inputstructure_dict.keys())[0]

        model_data = config_model.ModelConfig(config_data=config_data)
        self.model_list_name = 'model_list'
        self.model_list_value = model_data.get_model_list()

        parameter_data = config_generation.GenerationConfig(config_data=config_data).parameter_value
        parameter_dict = parameter_data.repeat_parameter_value

        self.repeat_parameter_name = parameter_data.repeat_parameter_name
        self.parameter_dict = {}

        for key, value in parameter_dict.items():
            value = value.split(',')
            self.parameter_dict[key] = self.deep_get(config_data, value)
    
        return self
    
    
    def write(self):
        return {self.name: {self.protein_name: self.protein_value,
                            self.model_list_name: self.model_list_value,
                            self.repeat_parameter_name: self.parameter_dict}}
    
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