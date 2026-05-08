
from scramblebench.script.config_preparation import config_constant
from pathlib import Path
from typing import Any
import logging

logger = logging.getLogger(__name__)


class GenerationParameterConfig:

    def __init__(self, parameter_data: dict[str, Any]):
        parameter_data = parameter_data[config_constant.GENERATION_PARAMETER_KEY]

        self.box_size_name = 'box_size'
        self.num_sample_name = 'num_sample'
        self.title_name = 'name'
        self.repeat_parameter_name = 'repeat_parameter'

        if parameter_data.get(self.box_size_name) is None:
            logging.warning('Unspecified parameter box size. Setting default to 16 A')
            self.box_size_value = 10
        elif isinstance(parameter_data[self.box_size_name], (int, float)):
            self.box_size_value = parameter_data[self.box_size_name]
        elif isinstance(parameter_data[self.box_size_name], str):
            self.box_size_value = [int(num.strip()) for num in parameter_data[self.box_size_name].split(',')]
        elif isinstance(parameter_data[self.box_size_name], list):
            self.box_size_value = parameter_data[self.box_size_name]
        else:
            raise ValueError(f'unsupported type of {parameter_data[self.box_size_name]}: {type(parameter_data[self.box_size_name])}')
        
        if parameter_data.get(self.num_sample_name) is None:
            logging.warning('Unspecified parameter num sample. Setting default to 100')
            self.num_sample_value = 100
        if isinstance(parameter_data[self.num_sample_name], int):
            self.num_sample_value = parameter_data[self.num_sample_name]
        elif isinstance(parameter_data[self.num_sample_name], str):
            self.num_sample_value = [int(num.strip()) for num in parameter_data[self.num_sample_name].split(',')]
        elif isinstance(parameter_data[self.num_sample_name], list):
            self.num_sample_value = parameter_data[self.num_sample_name]
        else:
            raise ValueError(f'unsupported type of {parameter_data[self.num_sample_name]}: {type(parameter_data[self.num_sample_name])}')
        
        self.title_value = parameter_data.get(self.title_name) or 'generic_title'
        self.repeat_parameter_value = parameter_data.get(self.repeat_parameter_name) or {}

    def update(self, key: str, value: str):
        if key == self.box_size_name:
            self.box_size_value = value
        elif key == self.num_sample_name:
            try:
                self.num_sample_value = [int(num.strip()) for num in str(value).split(',')]
            except ValueError:
                raise ValueError(f'num sample should be integer or list of integer, not {value}')
        elif key in self.title_name:
            self.title_value[key] = value
        
        elif key in self.repeat_parameter_name:
            self.repeat_parameter_value = value
        else:
            raise TypeError(f'no key called {key}')

        return self

    def validate_config(self):
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

    def write(self):
        return {config_constant.GENERATION_PARAMETER_KEY: {self.box_size_name : self.box_size_value,
                                                           self.num_sample_name : self.num_sample_value,
                                                           self.title_name : self.title_value,
                                                           self.repeat_parameter_name: self.repeat_parameter_value}}
    
class GenerationConfig:

    def __init__(self, config_data: dict[str, Any]):
        config_data = config_data[config_constant.GENERATION_KEY]

        self.name = config_constant.GENERATION_KEY
        self.input_name = 'input'
        self.output_name = 'output'
        self.script_name = 'script_pathfile'

        self.input_value = config_data[self.input_name]
        self.output_value = config_data[self.output_name]
        self.parameter_value = GenerationParameterConfig(config_data)

        #todo: create tests for default or incorrect generation script
        self.script_value = config_data.get(self.script_name) or config_constant.GENERATION_TEMPLATE_SCRIPT_PATH

    def update(self, key: str, value: str):
        if key == self.input_name:
            self.input_value = value
        elif key == self.output_name:
            self.output_value = value
        else:
            self.parameter_value = self.parameter_value.update(key, value)

        return self

    def validate_config(self):
        check_file_start_with_number(self.input_value)
        check_file_start_with_number(self.output_value)
        
        assert isinstance(self.input_value, str)
        assert isinstance(self.output_value, str)

        self.parameter_value.validate_config()
        
        if not Path(self.script_value).is_file():
            raise FileNotFoundError(f'{self.script_value} was not found. Please check your directory again')
        if not Path(self.script_value).suffix == '.sh':
            raise ValueError(f'Script generation file must be in .sh format. Instead, we detect {self.script_value}')

    def write(self):

        return {config_constant.GENERATION_KEY : {self.input_name: str(Path(self.input_value).resolve()),
                                                  self.output_name: str(Path(self.output_value).resolve()),
                                                  self.script_name: str(Path(self.script_value).resolve())} | 
                                                  self.parameter_value.write()}

def check_file_start_with_number(fname: str):
    try:
        float(Path(fname).name[0])
        raise TypeError(f'please do not use number as starting filename. change {fname}')
    except ValueError:
        pass
    

