"""This file handles the post generation key of user's YAML config. 
Mainly used for p1_generate_config.py and other downstream path to
access the subkey field easily.
"""
import logging

from pathlib import Path
from typing import Any, Optional
from typing_extensions import Self

from scramblebench.script.config_preparation import config_constant
from scramblebench.script.config_preparation.config_model import ModelConfig
from scramblebench.script.config_preparation.config_parameter import ParameterConfig
from scramblebench.script.config_preparation.config_utils import check_file_start_with_number

logger = logging.getLogger(__name__)


class PostGenerationConfig:
    """class containing the information regarding post generation. The class also relies on
    ModelConfig and ParameterConfig to check the validity of the model that is provided to this 
    class.
    """


    def __init__(self, config_data: dict[str, Any]):
        """initialize class

        Args:
            config_data (dict[str, Any]): user's prepared YAML config with post generation key as dictionary

        Raises:
            KeyError: if config dictionary does not contain the post_generation key
        """

        self.name = config_constant.POST_GENERATION_KEY

        if self.name not in config_data:
            raise KeyError(f'You forgot to put "{self.name}" key in your config file!')

        post_generation_data = config_data[self.name]

        self.input_name = 'input'
        self.output_name = 'output'
        self.pick_random_name = 'pick_random'
        self.pick_last_name = 'pick_last'

        self.pick_random_value = None
        self.pick_last_value = None
        self.input_value = post_generation_data[self.input_name]

        self.output_value = post_generation_data[self.output_name]

        if post_generation_data.get(self.pick_random_name):
            self.pick_random_value = [model.strip() for model in \
                                      post_generation_data.get(self.pick_random_name).split(',')]
        if post_generation_data.get(self.pick_last_name):
            self.pick_last_value = [model.strip() for model in \
                                    post_generation_data.get(self.pick_last_name).split(',')]

        self.reference_model = []
        if config_constant.PARAMETER_KEY in config_data:
            self.reference_model = ParameterConfig(config_data).model_list_value
        else:
            self.reference_model = ModelConfig(config_data).get_model_list()


    def update(self, key: str, value: str) -> Self:
        """update value of a key

        Args:
            key (str): valid PostGenerationConfig key
            value (str): new value that will be updated

        Raises:
            TypeError: invalid PostGenerationConfig key

        Returns:
            Self: PostGenerationConfig
        """
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
        """check whether the model in the class chosen is valid based on the provided models in ModelConfig class

        Raises:
            ValueError: if model in pick_last is not found in the ModelConfig
            ValueError: if model is listed in both pick_last and pick_random
            ValueError: if model in pick_random is not found in the ModelConfig
            ValueError: if model is listed in both pick_last and pick_random
        """
        validated_model = []

        logging.info(f'Model name detected in the config file: {self.reference_model}')
        if self.pick_last_value:
            for model in self.pick_last_value:
                if model not in self.reference_model:
                    raise ValueError(f'Model {model} is not found on the {config_constant.MODEL_KEY} parameter')
                if model in validated_model:
                    raise ValueError(f'Model {model} has at least been listed twice in pick random and pick last')
                validated_model.append(model)

        if self.pick_random_value:
            for model in self.pick_random_value:
                if model not in self.reference_model:
                    raise ValueError(f'Model {model} is not found on the {config_constant.MODEL_KEY} parameter')
                if model in validated_model:
                    raise ValueError(f'The model {model} has at least been listed twice in pick random and pick last')
                validated_model.append(model)


    def validate_config(self):
        """check if PostGenerationConfig is valid.
        """
        check_file_start_with_number(self.input_value)
        check_file_start_with_number(self.output_value)

        assert isinstance(self.input_value, str)
        assert isinstance(self.output_value, str)

        self.validate_model_pick()


    def update_input_output(self, prefix_dir: str) -> None:
        """update the input and output folder for post generation. This is utilized in p1_generate_config.py,
        because user may provide multiple batch parameters. Hence, the output folder directory path depends on these
        parameters.

        Args:
            prefix_dir (str): the new directory that will replace/append the previous directory
        """

        if '/' in self.input_value[0]:
            logging.warning('prefix dir for generation is provided. \
                            Ignoring absolute path provided by input key and only take the last dir')
            self.input_value = Path(prefix_dir) / Path(self.input_value).name
        else:
            self.input_value = Path(prefix_dir) / self.input_value

        if '/' in self.output_value[0]:
            logging.warning('prefix dir for generation is provided. \
                            Ignoring absolute path provided by output key and only take the last dir')
            self.output_value = Path(prefix_dir) / Path(self.output_value).name
        else:
            self.output_value = Path(prefix_dir) / self.output_value


    def write(self, prefix_dir: Optional[str]=None) -> dict[str, dict]:
        """Write the standardised config format for PostGenerationConfig

        Args:
            prefix_dir (Optional[str], optional): the new directory that will replace/append the previous directory. 
            Defaults to None.

        Returns:
            dict[str, dict]: PostGenerationConfig as dictionary
        """

        post_generation_dict_data = {}

        if prefix_dir:
            self.update_input_output(prefix_dir=prefix_dir)

        post_generation_dict_data = {self.input_name: str(self.input_value),
                                     self.output_name: str(self.output_value)}

        if self.pick_last_value:
            post_generation_dict_data[self.pick_last_name] = ','.join(self.pick_last_value)
        if self.pick_random_value:
            post_generation_dict_data[self.pick_random_name] = ','.join(self.pick_random_value)

        return {config_constant.POST_GENERATION_KEY : post_generation_dict_data}
