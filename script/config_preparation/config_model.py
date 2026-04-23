
from typing import Any
from pathlib import Path

from ScrambleBench.script.error_handler import FileDataError, FileTypeError, DirNotFound
import subprocess

class ModelStructure:

    def __init__(self, model_dict: dict[str, str]):

        self.model_name = 'name'
        self.dir_name = 'dir'
        self.conda_env_name = 'conda_env'
    
        self.model_name_value = model_dict[self.model_name]
        self.dir_value = model_dict[self.dir_name]
        self.conda_env_value = model_dict[self.conda_env_name]
    
    def update(self, key: str, value: str):
        if key == self.model_name:
            self.model_name_value = value
        elif key == self.dir_name:
            self.dir_value = value
        elif key == self.conda_env_name:
            self.conda_env_value = value
        else:
            raise TypeError(f'no key called {key}')

        return self 
    
class ModelConfig:

    def __init__(self, config_data: dict[str, Any]):
        config_data = config_data['model']

        self.modelstructure_dict = {model_keyname: ModelStructure(model_data) for model_keyname, model_data in config_data.items()}

    def __iter__(self):
        """Returns the iterator object itself."""
        return iter(list(self.modelstructure_dict.items()))
    
    def update(self, model_keyname, model_key, value):
        try:
            self.modelstructure_dict[model_keyname] = self.modelstructure_dict[model_keyname].update(model_key, value)
        except KeyError:
            raise KeyError(f'no key called {model_keyname}')

        return self
    

class ModelConfig:

    def __init__(self, config_data: dict[str, Any]):
        config_data = config_data['model']

        self.modelstructure_dict = {model_keyname: ModelStructure(model_data) for model_keyname, model_data in config_data.items()}

    def __iter__(self):
        """Returns the iterator object itself."""
        return iter(list(self.modelstructure_dict.items()))
    
    def update(self, model_keyname, model_key, value):
        try:
            self.modelstructure_dict[model_keyname] = self.modelstructure_dict[model_keyname].update(model_key, value)
        except KeyError:
            raise KeyError(f'no key called {model_keyname}')

        return self
    

def check_conda_env(conda_env: str) -> None:

    # Run 'conda env list' command and capture the output
    output = subprocess.run(
        "conda env list | awk '{ print $1}'",
        shell=True,
        capture_output=True,
        text=True,
        check=True
    ).stdout.split('\n')

    if not conda_env in output:
        if conda_env == 'not_applicable' or conda_env is None:
            print('Please note that your conda is explicitly left out')
        else:
            raise ValueError(f'Conda environment {conda_env} does not exist')
    

def check_modelstruct(modelstruct_data: ModelStructure) -> None:
    if not isinstance(modelstruct_data.model_name_value, str):
        raise TypeError(f'{modelstruct_data.model_name_value} is preferably set to a string')

    if not Path(modelstruct_data.dir_value).is_dir():
        if modelstruct_data.dir_value is None or modelstruct_data.dir_value == 'not_applicable':
            print(f'Please note that the directory of {modelstruct_data.model_name_value} is left out for generation')
        else:   
            raise DirNotFound(f'Directory {modelstruct_data.dir_value=} does not exist')

    check_conda_env(modelstruct_data.conda_env_value)
    

def check_model(model_data: ModelConfig) -> None:
    for model in model_data:
        check_modelstruct(model[-1]) 
