"""This file handles the model key of user's YAML config. 
Mainly used for p1_generate_config.py and other downstream path to
access the subkey field easily.
"""

import logging

from pathlib import Path
from collections.abc import Iterable, ItemsView
from typing import Any
from typing_extensions import Self

from scramblebench.script.utils.error_handler import DirNotFoundError
from scramblebench.script.config_preparation import config_constant
from scramblebench.script.config_preparation.config_utils import check_conda_env

logger = logging.getLogger(__name__)

class ModelStructure:
    """class containing each model information wrapped inside ModelConfig class
    """

    def __init__(self, model_dict: dict[str, str]):
        """initialize the class to save model name, directory name (typically Github repo),
        and conda env

        Args:
            model_dict (dict[str, str]): user's prepared YAML config with model key as dictionary
        """

        self.model_name = 'name'
        self.dir_name = 'dir'
        self.conda_env_name = 'conda_env'

        self.model_name_value = model_dict[self.model_name]
        self.dir_value = model_dict.get(self.dir_name)
        self.conda_env_value = model_dict.get(self.conda_env_name)


    def update(self, key: str, value: str) -> Self:
        """update value of key of the class

        Args:
            key (str): valid ModelStructure key
            value (str): new value that will be updated

        Raises:
            TypeError: invalid key of ModelStructure

        Returns:
            Self: ModelStructure
        """
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
        """validate whether the user ModelStructure is valid.

        Raises:
            TypeError: if the model name is not a string
            DirNotFoundError: if the model directory does not exist
        """
        if not isinstance(self.model_name_value, str):
            raise TypeError(f'{self.model_name_value} is preferably set to a string')

        if not Path(self.dir_value).is_dir():
            if self.dir_value is None or self.dir_value == 'not_applicable':
                logging.warning(f'Please note that the directory of {self.model_name_value} is left out for generation')
            else:
                raise DirNotFoundError(f'Directory {self.dir_value=} does not exist')

        check_conda_env(self.conda_env_value)


    def write(self) -> dict[str, dict]:
        """Write the standardised config format for ModelStructure

        Returns:
            dict[str, dict]: ModelStructure as dictionary
        """
        if self.dir_value is None or self.dir_value == 'not_applicable':
            dirname = 'not_applicable'
        else:
            dirname = str(Path(self.dir_value).resolve())

        return { self.model_name: self.model_name_value,
                 self.dir_name: dirname,
                 self.conda_env_name: self.conda_env_value}


class ModelConfig:
    """class containing model name and the corresponding ModelStructure
    """

    def __init__(self, config_data: dict[str, Any]):
        """initialize the class. ModelStructure will not know the model name written
        in the subkey of model config, defined as the "model_keyname" variable

        Args:
            config_data (dict[str, Any]): user's YAML config with model key as dictionary

        Raises:
            KeyError: if model key != config_constant.MODEL_KEY
        """

        self.name = config_constant.MODEL_KEY

        if self.name not in config_data:
            raise KeyError(f'You forgot to put "{self.name}" key in your config file!')

        model_data = config_data[config_constant.MODEL_KEY]

        self.modelstructure_dict = {model_keyname: ModelStructure(structmodel_data) \
                                    for model_keyname, structmodel_data in model_data.items()}


    def __iter__(self) -> Iterable[ItemsView[str, ModelStructure]]:
        """custom iteration method to return the list of ModelStructure as dictionary.
        
        Returns:
            Iterable[ItemsView[str, ModelStructure]]: list of user's (model name, ModelStructure)
        """
        return iter(list(self.modelstructure_dict.items()))


    def update(self, model_keyname, model_key, value) -> Self:
        """update value of key of the class

        Args:
            model_keyname (_type_): model name (not to be confused with ModelStructure.model_name)
            model_key (_type_): valid key name of ModelStructure
            value (_type_): new value that will be updated

        Raises:
            KeyError: if model_key is not valid

        Returns:
            Self: ModelConfig
        """
        try:
            self.modelstructure_dict[model_keyname] = self.modelstructure_dict[model_keyname].update(model_key, value)
        except KeyError as exc:
            raise KeyError(f'no key called {model_keyname}') from exc

        return self


    # see __iter__ method for how this is possible.
    def validate_config(self) -> None:
        """check whether ModelConfig is valid by looping all ModelStructure
        """
        for model in self:
            (model[-1].validate_config())


    def get_model_list(self) -> list[str]:
        """get list of model name

        Returns:
            list[str]: list of ModelStruct.model_name_value
        """
        return [modelstruct.model_name_value for modelstruct in self.modelstructure_dict.values()]


    def write(self) -> dict[str, Any]:
        """Return standardised dictionary for ModelConfig from list of ModelStructure

        Returns:
            dict[str, Any]: standardised dictionary for user's model config
        """
        model_data = {}
        for key, modelstruct in self.modelstructure_dict.items():
            model_data[key] = modelstruct.write()

        return {config_constant.MODEL_KEY: model_data}
