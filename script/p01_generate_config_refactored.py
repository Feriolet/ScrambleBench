
from typing import Any
from pathlib import Path

import yaml
from ScrambleBench.script.error_handler import FileDataError, FileTypeError, DirNotFound
from ScrambleBench.script.config_preparation.config_input import InputConfig, check_input
from ScrambleBench.script.config_preparation.config_model import ModelConfig, check_model

def load_config(config_fname: str) -> dict[str, Any]:
    with open(config_fname, 'r') as config_fn:
        return yaml.safe_load(config_fn)
    
def read_config(config_data: dict[str, Any]):
    config_key_class_dict = {
        'input': [InputConfig, check_input],
        'model': [ModelConfig, check_model]
    }

    inject_config_class_list = []
    for config_key in config_data.keys():
        if config_key in config_key_class_dict:
            inject_config_class_list.append(config_key_class_dict[config_key])
    
    return inject_config_class_list

if __name__ == '__main__':
    configf = '/opt/veincent/GenAI_manuscript/ScrambleBench/test/config.yml'
    data_input = yaml.safe_load(open(configf, 'r'))

    print(read_config(data_input))