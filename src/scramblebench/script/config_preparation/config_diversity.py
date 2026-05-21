
from scramblebench.script.config_preparation import config_constant
from pathlib import Path
from typing import Any
import logging
import subprocess


logger = logging.getLogger(__name__)


class DiversityConfig:
    def __init__(self, config_data: dict[str, Any]):
        config_data = config_data[config_constant.ANALYSIS_DIVERSITY_KEY]
        self.input_name = 'input'
        self.output_name = 'output'
        self.environment_name = 'conda_env'

        self.method_name = 'method'
        self.method_distance_name = 'distance'
        self.method_diversity_name = 'diversity'

        self.input_value = config_data[self.input_name]
        self.output_value = config_data[self.output_name]
        self.method_value = config_data.get(self.method_name) 
        self.environment_value = config_data[self.environment_name]

        if config_data.get(self.method_name) is None:
            self.method_value = [{self.method_distance_name: 'ecfp',
                                  self.method_diversity_name: 'hamdiv'}]

    def update(self, key: str, value: str):
        return self

    def validate_config(self):
        check_file_start_with_number(self.input_value)
        check_file_start_with_number(self.output_value)
        check_conda_env(self.environment_value)

        if not isinstance(self.input_value, str):
            raise TypeError(f'Please put input folder as string, not {type(self.input_value)}')
        if not isinstance(self.output_value, str):
            raise TypeError(f'Please put output folder as string, not {type(self.output_value)}')
        

        if not isinstance(self.method_value, list):
            raise ValueError(f'Your value in the diversity.method section is not a list of dictionary, but rather {type(self.method_value)}')
        for method in self.method_value:
            if not isinstance(method, dict):
                raise ValueError(f'Your value in the diversity.method section is not a list of dictionary, but rather list of {type(method)}')
            for key, value in method.items():
                if key == self.method_diversity_name:
                    self.validate_diversity_metric(diversity=value)
                elif key == self.method_distance_name:
                    self.validate_distance_metric(distance=value)
                else:
                    raise KeyError(f'unsupported key {key} in the diversity method config')

    def validate_diversity_metric(self, diversity):

        SUPPORTED_DIVERSITY_METRIC = ['average', 'hamdiv', 'richness', 'fg', 'rs', 'bm', 
                                      'generic_bm', 'intdiv', 'sumdiv', 'diam', 'sumdiam', 'bot',
                                      'sumbot', 'dpp']

        if not isinstance(diversity, str):
            raise ValueError(f'You diversity value must be a string, not {type(diversity)}')
        if diversity.lower() not in SUPPORTED_DIVERSITY_METRIC:
            raise ValueError(f'The value of the diversity metric {diversity} is not supported! Please key in either: {",".join(SUPPORTED_DIVERSITY_METRIC)}')
    
    def validate_distance_metric(self, distance):

        SUPPORTED_DISTANCE_METRIC = ['mces', 'ecfp', None] 

        if distance is not None:
            if not isinstance(distance, str):
                raise ValueError(f'You distance value must be a string or nonetype, not {type(distance)}')
            if distance.lower() not in SUPPORTED_DISTANCE_METRIC:
                raise ValueError(f'The value of the distance metric {distance} is not supported! Please key in either: {",".join(SUPPORTED_DISTANCE_METRIC)}')
                
    def update_input_output(self, prefix_dir):
        if '/' in self.input_value[0]:
            logging.warning('prefix dir for generation is provide. Ignoring absolute path provided by input key and only take the last dir')
            self.input_value = Path(prefix_dir) / Path(self.input_value).name
        else:
            self.input_value = Path(prefix_dir) / self.input_value
        
        if '/' in self.output_value[0]:
            logging.warning('prefix dir for generation is provide. Ignoring absolute path provided by output key and only take the last dir')
            self.output_value = Path(prefix_dir) / Path(self.output_value).name
        else:
            self.output_value = Path(prefix_dir) / self.output_value  

    def write(self, prefix_dir=None):
        if prefix_dir:
            self.update_input_output(prefix_dir=prefix_dir)
            
        return {config_constant.ANALYSIS_DIVERSITY_KEY: {self.input_name: str(self.input_value),
                                                          self.output_name: str(self.output_value),
                                                          self.environment_name: self.environment_value,
                                                          self.method_name: self.method_value}}

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
            logging.warning('Please note that some of your conda is explicitly left out')
        else:
            raise ValueError(f'Conda environment {conda_env} does not exist')
        
def check_file_start_with_number(fname: str):
    try:
        float(Path(fname).name[0])
        raise TypeError(f'please do not use number as starting filename. change {fname}')
    except ValueError:
        pass
