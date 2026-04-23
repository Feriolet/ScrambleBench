import sys
from os.path import dirname as d
from os.path import abspath, join
root_dir = d(d(abspath(__file__)))
sys.path.append(root_dir)

import pytest
from typing import NamedTuple


class ConfigInputKey(NamedTuple):
    name: str
    complex_name: str
    pdb_name: str
    sdf_name: str
    title_name: str

@pytest.fixture
def input_key() -> ConfigInputKey:
    return ConfigInputKey(
        name = 'input',
        complex_name = 'complex_path',
        pdb_name = 'pdb_path',
        sdf_name = 'sdf_path',
        title_name = 'protein_title',
    )