
from scramblebench.script.config_preparation import config_constant
from pathlib import Path
from typing import Any

class GenerationParameterConfig:

    def __init__(self, parameter_data: dict[str, Any]):
        parameter_data = parameter_data[config_constant.GENERATION_PARAMETER_KEY]

        self.box_size_name = 'box_size'
        self.num_sample_name = 'num_sample'
        self.title_name = 'name'

        self.box_size_value = parameter_data[self.box_size_name]
        self.num_sample_value = [int(num.strip()) for num in parameter_data[self.num_sample_name].split(',')]
        self.title_value = parameter_data[self.title_name]

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
        else:
            raise TypeError(f'no key called {key}')

        return self

    def validate_config(self):
        assert isinstance(self.parameter_value.box_size_value, (int, float))

        if isinstance(self.parameter_value.num_sample_value, list):
            if not all(isinstance(num, int) for num in self.parameter_value.num_sample_value):
                raise ValueError(f'{self.parameter_value.num_sample_value} needs to be an integer or list of integers')
        elif isinstance(self.parameter_value.num_sample_value, int):
            pass
        else:
            raise ValueError(f'{self.parameter_value.num_sample_value} needs to be an integer or list of integers')
        
class GenerationConfig:

    def __init__(self, config_data: dict[str, Any]):
        config_data = config_data[config_constant.GENERATION_KEY]

        self.name = config_constant.GENERATION_KEY
        self.input_name = 'input'
        self.output_name = 'output'
    
        self.input_value = config_data[self.input_name]
        self.output_value = config_data[self.output_name]
        self.parameter_value = GenerationParameterConfig(config_data)

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

        assert isinstance(self.parameter_value.box_size_value, (int, float))

        if isinstance(self.parameter_value.num_sample_value, list):
            if not all(isinstance(num, int) for num in self.parameter_value.num_sample_value):
                raise ValueError(f'{self.parameter_value.num_sample_value} needs to be an integer or list of integers')
        elif isinstance(self.parameter_value.num_sample_value, int):
            pass
        else:
            raise ValueError(f'{self.parameter_value.num_sample_value} needs to be an integer or list of integers')
        

def check_file_start_with_number(fname: str):
    try:
        float(Path(fname).name[0])
        raise TypeError(f'please do not use number as starting filename. change {fname}')
    except ValueError:
        pass
    

