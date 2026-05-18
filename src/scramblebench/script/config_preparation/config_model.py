
from typing import Any
from pathlib import Path

from scramblebench.script.utils.error_handler import DirNotFoundError
from scramblebench.script.config_preparation import config_constant
import subprocess
import logging

logger = logging.getLogger(__name__)

class ModelStructure:

    def __init__(self, model_dict: dict[str, str]):

        self.model_name = 'name'
        self.dir_name = 'dir'
        self.conda_env_name = 'conda_env'
    
        self.model_name_value = model_dict[self.model_name]
        self.dir_value = model_dict.get(self.dir_name)
        self.conda_env_value = model_dict.get(self.conda_env_name)
    
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
    
    def validate_config(self) -> None:
        if not isinstance(self.model_name_value, str):
            raise TypeError(f'{self.model_name_value} is preferably set to a string')

        if not Path(self.dir_value).is_dir():
            if self.dir_value is None or self.dir_value == 'not_applicable':
                logging.warning(f'Please note that the directory of {self.model_name_value} is left out for generation')
            else:   
                raise DirNotFoundError(f'Directory {self.dir_value=} does not exist')

        check_conda_env(self.conda_env_value)

    def write(self):
        if self.dir_value is None or self.dir_value == 'not_applicable':
            dirname = 'not_applicable'
        else:
            dirname = str(Path(self.dir_value).resolve())

        return { self.model_name: self.model_name_value,
                 self.dir_name: dirname,
                 self.conda_env_name: self.conda_env_value}

class ModelConfig:

    def __init__(self, config_data: dict[str, Any]):

        self.name = config_constant.MODEL_KEY
        model_data = config_data[config_constant.MODEL_KEY]

        self.modelstructure_dict = {model_keyname: ModelStructure(structmodel_data) for model_keyname, structmodel_data in model_data.items()}

    def __iter__(self):
        """Returns the iterator object itself."""
        return iter(list(self.modelstructure_dict.items()))
    
    def update(self, model_keyname, model_key, value):
        try:
            self.modelstructure_dict[model_keyname] = self.modelstructure_dict[model_keyname].update(model_key, value)
        except KeyError:
            raise KeyError(f'no key called {model_keyname}')

        return self
    
    # see __iter__ method for how this is possible.
    def validate_config(self) -> None:
        for model in self:
            (model[-1].validate_config()) 

    def get_model_list(self):
        return [modelstruct.model_name_value for modelstruct in self.modelstructure_dict.values()]

    def write(self) -> dict[str, Any]:
        model_data = {}
        for key, modelstruct in self.modelstructure_dict.items():
            model_data[key] = modelstruct.write()
            
        return {config_constant.MODEL_KEY: model_data}
    
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
            logging.warning('Please note that your conda is explicitly left out')
        else:
            raise ValueError(f'Conda environment {conda_env} does not exist')
    


    


