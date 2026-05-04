
from scramblebench.script.config_preparation import config_constant
from scramblebench.script.config_preparation.config_model import ModelConfig

from pathlib import Path
from typing import Any
import logging

logger = logging.getLogger(__name__)

class PostGenerationConfig:

    def __init__(self, config_data: dict[str, Any]):
        post_generation_data = config_data[config_constant.POST_GENERATION_KEY]

        self.name = config_constant.POST_GENERATION_KEY
        self.input_name = 'input'
        self.output_name = 'output'
        self.pick_random_name = 'pick_random'
        self.pick_last_name = 'pick_last'
    
        self.pick_random_value = None
        self.pick_last_value = None
        self.input_value = post_generation_data[self.input_name]

        self.output_value = post_generation_data[self.output_name]

        if post_generation_data.get(self.pick_random_name):
            self.pick_random_value = [model.strip() for model in post_generation_data.get(self.pick_random_name).split(',')]
        if post_generation_data.get(self.pick_last_value):
            self.pick_last_value = [model.strip() for model in post_generation_data.get(self.pick_last_name).split(',')]

        self.reference_model = ModelConfig(config_data)

    def update(self, key: str, value: str):
        if key == self.input_name:
            self.input_value = value
        elif key == self.output_name:
            self.output_value = value
        elif key == self.pick_random_name:
            self.pick_random_value = value
        elif key == self.pick_last_name:
            self.pick_last_value = value
        else:
            raise TypeError(f'no key called {key}')

        return self

    def validate_model_pick(self):
        validated_model = []

        model_reference_list = self.reference_model.get_model_list()
        logging.info(f'Model name detected in the config file: {model_reference_list}')
        if self.pick_last_value:
            for model in self.pick_last_value:
                if model not in model_reference_list:
                    raise ValueError(f'The model {model} is not found on the {config_constant.MODEL_KEY} name parameter')
                if model in validated_model:
                    raise ValueError(f'The model {model} has at least been listed twice in pick random and pick last')
                validated_model.append(model)

        if self.pick_random_value:
            for model in self.pick_random_value:
                if model not in model_reference_list:
                    raise ValueError(f'The model {model} is not found on the {config_constant.MODEL_KEY} name parameter')
                if model in validated_model:
                    raise ValueError(f'The model {model} has at least been listed twice in pick random and pick last')
                validated_model.append(model)

    def validate_config(self):
        check_file_start_with_number(self.input_value)
        check_file_start_with_number(self.output_value)
        
        assert isinstance(self.input_value, str)
        assert isinstance(self.output_value, str)

        self.validate_model_pick()
        
def check_file_start_with_number(fname: str):
    try:
        float(Path(fname).name[0])
        raise TypeError(f'please do not use number as starting filename. change {fname}')
    except ValueError:
        pass
    
