
from typing import Any, Callable
from pathlib import Path
from dataclasses import dataclass
import yaml

from scramblebench.script.config_preparation.config_input import InputConfig
from scramblebench.script.config_preparation.config_model import ModelConfig
from scramblebench.script.config_preparation.config_generation import GenerationConfig
from scramblebench.script.config_preparation.config_post_generation import PostGenerationConfig
from scramblebench.script.config_preparation.config_genbench3d import GenBench3DConfig
from scramblebench.script.config_preparation import config_constant


def load_config(config_fname: str) -> dict[str, Any]:
    with open(config_fname, 'r') as config_fn:
        return yaml.safe_load(config_fn)


def check_correct_input_output_folder(prestep: Callable, poststep: Callable) -> None:

    if not str(prestep.output_value) == str(poststep.input_value):
        raise TypeError(f'The output folder of {prestep.name} and input folder of {poststep.name} is not the same')


def validate_config(config_data: dict[str, Any]) -> None:

    compulsory_config_key_list = [
        config_constant.INPUT_KEY,
        config_constant.MODEL_KEY,
        config_constant.GENERATION_KEY,
        config_constant.POST_GENERATION_KEY
    ]

    for config_key in compulsory_config_key_list:
        if config_key not in config_data.keys():
            raise ValueError(f'{config_key} not in your config file. Please add it with case sensitive letters.')
        
    InputConfig(config_data).validate_config()
    ModelConfig(config_data).validate_config()
    GenerationConfig(config_data).validate_config()
    PostGenerationConfig(config_data).validate_config()

    check_correct_input_output_folder(prestep=GenerationConfig(config_data),
                                      poststep=PostGenerationConfig(config_data))
    
    if config_constant.GENBENCH3D_KEY in config_data.keys():
        GenBench3DConfig(config_data).validate_config()
        check_correct_input_output_folder(prestep=PostGenerationConfig(config_data),
                                        poststep=GenBench3DConfig(config_data))
        
if __name__ == '__main__':
    configf = '/opt/veincent/GenAI_manuscript/ScrambleBench/src/scramblebench/test/config.yml'
    data_input = yaml.safe_load(open(configf, 'r'))

    print(validate_config(data_input))