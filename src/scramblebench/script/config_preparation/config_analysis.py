from typing import Any, Optional
from pathlib import Path

from scramblebench.script.utils.error_handler import FileDataError, FileTypeError
from scramblebench.script.config_preparation import config_constant
from scramblebench.script.config_preparation.config_genbench3d import GenBench3DConfig

import Bio
import rdkit
from Bio.PDB import PDBParser
from rdkit import Chem
import os
from scramblebench.script.utils.split_protein_ligand import split_pocket_ligand
from oddt.toolkits.extras.rdkit import fixer
import numpy as np
import logging
import re


logger = logging.getLogger(__name__)

class AnalysisConfig:

    def __init__(self, config_data: dict[str, Any]):

        self.name = config_constant.ANALYSIS_KEY
        analysis_data = config_data[self.name]

        self.genbench_name = config_constant.ANALYSIS_GENBENCH3D_KEY
        self.genbench_value = GenBench3DConfig(analysis_data[self.genbench_name])


    def update(self, key, value):
        if key == self.genbench_name:
            self.genbench_value = self.genbench_value.update(value)

    def validate_config(self):

        self.genbench_value.validate_config()

    def write(self):

        return {config_constant.ANALYSIS_KEY   : {self.genbench_name: self.genbench_value.write()}}
