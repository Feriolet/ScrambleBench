from typing import Any, Optional
from pathlib import Path

from scramblebench.script.utils.error_handler import FileDataError, FileTypeError
from scramblebench.script.config_preparation import config_constant
from scramblebench.script.config_preparation.config_genbench3d import GenBench3DConfig
from scramblebench.script.config_preparation.config_redocking import RedockingConfig
from scramblebench.script.config_preparation.config_diversity import DiversityConfig

import Bio
import rdkit
from Bio.PDB import PDBParser
from rdkit import Chem
import os
from oddt.toolkits.extras.rdkit import fixer
import numpy as np
import logging
import re


logger = logging.getLogger(__name__)

class AnalysisConfig:

    def __init__(self, config_data: dict[str, Any]):

        self.name = config_constant.ANALYSIS_KEY
        analysis_data = config_data[self.name]

        self.valid_key_list = []
        if config_constant.ANALYSIS_GENBENCH3D_KEY in analysis_data.keys():
            self.valid_key_list.append(GenBench3DConfig(analysis_data))
        if config_constant.ANALYSIS_REDOCKING_KEY in analysis_data.keys():
            self.valid_key_list.append(RedockingConfig(analysis_data))
        if config_constant.ANALYSIS_DIVERSITY_KEY in analysis_data.keys():
            self.valid_key_list.append(DiversityConfig(analysis_data))

    def update(self, key, value):
        pass

    def validate_config(self):
        for config in self.valid_key_list:
            config.validate_config()

    def write(self, prefix_dir=None):
        
        data = {}
        for config in self.valid_key_list:
            data = data | config.write(prefix_dir)

        return {self.name  : data}
