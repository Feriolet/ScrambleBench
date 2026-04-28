import pytest
import pandas as pd
import pathlib
from pathlib import Path
import subprocess
import yaml
from dataclasses import dataclass

from scramblebench.script.config_preparation.config_input import InputConfig
from scramblebench.script.config_preparation.config_model import ModelConfig
from scramblebench.script.config_preparation.config_generation import GenerationConfig
from scramblebench.script.config_preparation.config_post_generation import PostGenerationConfig
from scramblebench.script.config_preparation.config_genbench3d import GenBench3DConfig

from scramblebench.script.p01_generate_config_refactored import check_correct_input_output_folder
from scramblebench.script.p01_generate_config_refactored import load_config
from scramblebench.script.error_handler import FileTypeError, FileDataError, DirNotFound


pytest.input = 'input'
pytest.output = 'output'
pytest.genbench_dir = 'genbench_dir'
pytest.conda_env = 'conda_env'
pytest.schrodinger_dir = 'schrodinger_dir'
pytest.genbench_config = 'genbench_config'
pytest.do_complex_minimisation = 'do_complex_forcefield_minimisation'
pytest.do_docking_minimisation = 'do_docking_forcefield_minimisation'
pytest.cancel_protonation = 'do_cancel_protonation_by_obabel_or_adfr'

@pytest.fixture(scope='function')
def pytest_load_config():
    file_dir = str(Path(__file__).resolve().parent)
    config_data = f'{file_dir}/config.yml'
    return load_config(config_data)

@dataclass
class ConfigInputFile:
    correct_protein_pdb: str
    correct_complex_pdb: str
    correct_ligand_sdf: str

    complex_without_protein_pdb: str
    incorrect_protein_sdf: str
    incorrect_protein_pdb: str

    incorrect_ligand_sdf: str
    multiple_ligand_sdf: str
    translated_ligand_sdf: str

    empty_pdb: str
    empty_sdf: str

@dataclass
class ConfigInputKey:
    name: str
    complex_name: str
    pdb_name: str
    sdf_name: str
    title_name: str

@dataclass
class ConfigModelKey:
    name: str
    title_example_name: str
    title_name: str
    dir_name: str
    conda_name: str

@dataclass
class ConfigGenerationKey:
    name: str
    input: str
    numsample: str

@dataclass
class ConfigPostGenerationKey:
    pick_last: str
    pick_random: str
    input: str
    name: str

@dataclass
class ConfigGenBench3DKey:
    input: str
    output: str
    genbench_dir: str
    conda_env: str
    schrodinger_dir: str
    genbench_config: str
    do_complex_minimisation: str
    do_docking_minimisation: str
    cancel_protonation: str

@pytest.fixture(scope='session')
def input_file() -> ConfigInputFile:
    file_dir = str(Path(__file__).resolve().parent)
    return ConfigInputFile(
        correct_protein_pdb=f'{file_dir}/input_0/protein.pdb',
        correct_complex_pdb=f'{file_dir}/input_0/complex.pdb',
        correct_ligand_sdf=f'{file_dir}/input_0/ligand.sdf',
        complex_without_protein_pdb = f'{file_dir}/input_0/complex_without_protein.pdb',
        incorrect_protein_sdf = f'{file_dir}/input_0/incorrect_protein.sdf',
        incorrect_protein_pdb = f'{file_dir}/input_0/incorrect_protein.pdb',
        incorrect_ligand_sdf = f'{file_dir}/input_0/incorrect_ligand.sdf',
        multiple_ligand_sdf = f'{file_dir}/input_0/multiple_ligand.sdf',
        translated_ligand_sdf = f'{file_dir}/input_0/translated_ligand.sdf',
        empty_pdb = f'{file_dir}/input_0/empty.pdb',
        empty_sdf = f'{file_dir}/input_0/empty.sdf')

@pytest.fixture(scope='session')
def input_key() -> ConfigInputKey:
    return ConfigInputKey(
        name = 'input',
        complex_name = 'complex_path',
        pdb_name = 'pdb_path',
        sdf_name = 'sdf_path',
        title_name = 'protein_title',
    )

@pytest.fixture(scope='session')
def model_key() -> ConfigModelKey:
    return ConfigModelKey(
        name = 'model',
        title_example_name = 'pocket2mol',
        title_name = 'name',
        dir_name = 'dir',
        conda_name = 'conda_env',
    )

@pytest.fixture(scope='session')
def generation_key() -> ConfigGenerationKey:
    return ConfigGenerationKey(
    name = 'generation',
    input = 'input',
    numsample = 'num_sample'
    )

@pytest.fixture(scope='session')
def post_generation_key() -> ConfigPostGenerationKey:
    return ConfigPostGenerationKey(
    pick_last = 'pick_last',
    pick_random = 'pick_random',
    input = 'input',
    name = 'post_generation'
    )

@pytest.fixture(scope='session')
def genbench3d_key() -> ConfigGenBench3DKey:
    return ConfigGenBench3DKey(
    input = 'input',
    output = 'output',
    genbench_dir = 'genbench_dir',
    conda_env = 'conda_env',
    schrodinger_dir = 'schrodinger_dir',
    genbench_config = 'genbench_config',
    do_complex_minimisation = 'do_complex_forcefield_minimisation',
    do_docking_minimisation = 'do_docking_forcefield_minimisation',
    cancel_protonation = 'do_cancel_protonation_by_obabel_or_adfr'
    )

class TestInputConfig:

    def test_update_config_pass(self, pytest_load_config, input_key: ConfigInputKey):
        InputConfig(pytest_load_config).update(input_key.title_name, 'test').validate_config()

    def test_update_config_incorrect_type_error(self, pytest_load_config, input_key: ConfigInputKey):
        with pytest.raises(TypeError):
            InputConfig(pytest_load_config).update(input_key.complex_name, 20.5).validate_config()

    def test_load_config_incorrect_key_error(self, pytest_load_config, input_key: ConfigInputKey):

        # rename old key to new invalid key
        pytest_load_config[input_key.name]['hiiiiiii'] = pytest_load_config[input_key.name].pop(input_key.complex_name)
        print(pytest_load_config)
        with pytest.raises(KeyError):
            InputConfig(pytest_load_config)

    def test_update_config_incorrect_key_error(self, pytest_load_config):
        with pytest.raises(TypeError):
            InputConfig(pytest_load_config).update('hiiiii', 20.5).validate_config()

    def test_input_file_pass(self, pytest_load_config):
        InputConfig(pytest_load_config).validate_config()

    def test_complex_pdb_not_found(self, pytest_load_config, input_key: ConfigInputKey):
        data = InputConfig(pytest_load_config).update(input_key.complex_name, 'xyz1230987')

        with pytest.raises(FileNotFoundError):
            data.validate_config()

    def test_complex_pdb_wrong_format(self, pytest_load_config, input_key: ConfigInputKey, input_file: ConfigInputFile):
        data = InputConfig(pytest_load_config).update(input_key.complex_name, input_file.incorrect_protein_sdf)

        with pytest.raises(FileTypeError):
            data.validate_config()

    def test_complex_pdb_empty_ligand_error(self, pytest_load_config, input_key: ConfigInputKey, input_file: ConfigInputFile):
        data = InputConfig(pytest_load_config).update(input_key.complex_name, input_file.correct_protein_pdb)

        with pytest.raises(FileDataError):
            data.validate_config()

    def test_complex_pdb_empty_protein_error(self, pytest_load_config, input_key: ConfigInputKey, input_file: ConfigInputFile):
        data = InputConfig(pytest_load_config).update(input_key.complex_name, input_file.complex_without_protein_pdb)

        with pytest.raises(FileDataError):
            data.validate_config()

    def test_complex_pdb_empty_file_error(self, pytest_load_config, input_key: ConfigInputKey, input_file: ConfigInputFile):
        data = InputConfig(pytest_load_config).update(input_key.complex_name, input_file.empty_pdb)

        with pytest.raises(FileDataError):
            data.validate_config()


    def test_protein_pdb_not_found_error(self, pytest_load_config, input_key: ConfigInputKey):
        data = InputConfig(pytest_load_config).update(input_key.pdb_name, 'xyz1230987')

        with pytest.raises(FileNotFoundError):
            data.validate_config()

    def test_protein_pdb_wrong_format_error(self, pytest_load_config, input_key: ConfigInputKey, input_file: ConfigInputFile):
        data = InputConfig(pytest_load_config).update(input_key.pdb_name, input_file.incorrect_protein_sdf)

        with pytest.raises(FileTypeError):
            data.validate_config()

    def test_protein_pdb_contains_ligand_error(self, pytest_load_config, input_key: ConfigInputKey, input_file: ConfigInputFile):
        data = InputConfig(pytest_load_config).update(input_key.pdb_name, input_file.correct_complex_pdb)

        with pytest.raises(FileDataError):
            data.validate_config()

    def test_protein_pdb_mismatch_complex_error(self, pytest_load_config, input_key: ConfigInputKey, input_file: ConfigInputFile):
        data = InputConfig(pytest_load_config).update(input_key.pdb_name, input_file.incorrect_protein_pdb)

        with pytest.raises(FileDataError):
            data.validate_config()

    def test_protein_pdb_empty_file_error(self, pytest_load_config, input_key: ConfigInputKey, input_file: ConfigInputFile):
        data = InputConfig(pytest_load_config).update(input_key.pdb_name, input_file.empty_pdb)

        with pytest.raises(FileDataError):
            data.validate_config()


    def test_ligand_sdf_not_found_error(self, pytest_load_config, input_key: ConfigInputKey):
        data = InputConfig(pytest_load_config).update(input_key.sdf_name, 'xyz1230987')

        with pytest.raises(FileNotFoundError):
            data.validate_config()

    def test_ligand_sdf_wrong_format_error(self, pytest_load_config, input_key: ConfigInputKey, input_file: ConfigInputFile):
        data = InputConfig(pytest_load_config).update(input_key.sdf_name, input_file.empty_pdb)

        with pytest.raises(FileTypeError):
            data.validate_config()

    # def test_ligands_contains_protein_error(global_var):
    #     data = InputConfig(pytest_load_config)).update(pytest.sdf_name, pytest.correct_complex_pdb)

    #     with pytest.raises(FileDataError):
    #         data.validate_config()

    def test_ligand_sdf_empty_file_error(self, pytest_load_config, input_key: ConfigInputKey, input_file: ConfigInputFile):
        data = InputConfig(pytest_load_config).update(input_key.sdf_name, input_file.empty_sdf)

        with pytest.raises(FileDataError):
            data.validate_config()

    def test_ligands_contains_multiple_ligands_error(self, pytest_load_config, input_key: ConfigInputKey, input_file: ConfigInputFile):
        data = InputConfig(pytest_load_config).update(input_key.sdf_name, input_file.multiple_ligand_sdf)

        with pytest.raises(FileDataError):
            data.validate_config()

    def test_ligand_sdf_mismatch_complex_error(self, pytest_load_config, input_key: ConfigInputKey, input_file: ConfigInputFile):
        data = InputConfig(pytest_load_config).update(input_key.sdf_name, input_file.incorrect_ligand_sdf)

        with pytest.raises(FileDataError):
            data.validate_config()

    def test_ligand_sdf_not_in_pocket_error(self, pytest_load_config, input_key: ConfigInputKey, input_file: ConfigInputFile):
        data = InputConfig(pytest_load_config).update(input_key.sdf_name, input_file.translated_ligand_sdf)

        with pytest.raises(FileDataError):
            data.validate_config()

class TestModelConfig:
    def test_model_pass(self, pytest_load_config):
        ModelConfig(pytest_load_config)

    def test_modelcommercial_pass(self, pytest_load_config, model_key: ConfigModelKey):
        ModelConfig(pytest_load_config) \
        .update(model_key.title_example_name, model_key.dir_name, 'not_applicable') \
        .update(model_key.title_example_name, model_key.conda_name, 'not_applicable') \
        .validate_config()

    def test_load_config_incorrect_model_key_error(self, pytest_load_config, model_key: ConfigModelKey):
        pytest_load_config[model_key.name][model_key.title_example_name]['hiiiiiii'] = pytest_load_config[model_key.name][model_key.title_example_name].pop(model_key.title_name)
        with pytest.raises(KeyError):
            ModelConfig(pytest_load_config)

    def test_incorrect_update_model_title_not_as_str_int_or_float(self, pytest_load_config, model_key: ConfigModelKey):
        with pytest.raises(TypeError):
            ModelConfig(pytest_load_config) \
                .update(model_key.title_example_name, model_key.title_name, [1,2,3]) \
                .validate_config()

    def test_incorrect_update_model_dir(self, pytest_load_config, model_key: ConfigModelKey):
        with pytest.raises(DirNotFound):
            ModelConfig(pytest_load_config) \
                .update(model_key.title_example_name, model_key.dir_name, 'xyz12345') \
                .validate_config()

    def test_incorrect_update_conda_dir(self, pytest_load_config, model_key: ConfigModelKey):
        with pytest.raises(ValueError):
            ModelConfig(pytest_load_config) \
                .update(model_key.title_example_name, model_key.conda_name, 'xyz12345') \
                .validate_config()

class TestGenerationConfig:
    def test_correct_generation(self, pytest_load_config):
        GenerationConfig(pytest_load_config).validate_config()

    def test_correct_generation_update_parameter(self, pytest_load_config, generation_key: ConfigGenerationKey):
        GenerationConfig(pytest_load_config).update(generation_key.numsample, 500).validate_config()


    def test_incorrect_generation_numsample(self, pytest_load_config, generation_key: ConfigGenerationKey):
        with pytest.raises(ValueError):
            GenerationConfig(pytest_load_config).update(generation_key.numsample, 'hello').validate_config()
        
    def test_incorrect_generation_file_startwithnum(self, pytest_load_config, generation_key: ConfigGenerationKey):
        with pytest.raises(TypeError):
            GenerationConfig(pytest_load_config).update(generation_key.input, 'hi/5ht2c').validate_config()

class TestPostGenerationConfig:
    def test_correct_post_generation(self, pytest_load_config):
        PostGenerationConfig(pytest_load_config).validate_config()

    def test_incorrect_post_generation_conflict_model(self, pytest_load_config, post_generation_key: ConfigPostGenerationKey):
        with pytest.raises(ValueError):
            PostGenerationConfig(pytest_load_config) \
                .update(post_generation_key.pick_last, 'Pocket2Mol') \
                .update(post_generation_key.pick_random, 'Pocket2Mol') \
                .validate_config()

    def test_incorrect_post_generation_nonexisting_model(self, pytest_load_config, post_generation_key: ConfigPostGenerationKey):
        with pytest.raises(ValueError):
            PostGenerationConfig(pytest_load_config) \
                .update(post_generation_key.pick_last, 'DiffGui') \
                .validate_config()

    def test_correct_generation_and_post_generation_folder(self, pytest_load_config):
        check_correct_input_output_folder(prestep=GenerationConfig(pytest_load_config),
                                          poststep=PostGenerationConfig(pytest_load_config))

    def test_incorrect_generation_and_post_generation_folder(self, pytest_load_config, post_generation_key: ConfigPostGenerationKey):
        with pytest.raises(TypeError):
            check_correct_input_output_folder(prestep=GenerationConfig(pytest_load_config),
                                            poststep=PostGenerationConfig(pytest_load_config).update(post_generation_key.input, 'hahahhaha'))

    def test_incorrect_post_generation_file_startwithnum(self, pytest_load_config, post_generation_key: ConfigPostGenerationKey):
        with pytest.raises(TypeError):
            PostGenerationConfig(pytest_load_config) \
                .update(post_generation_key.input, 'hi/5ht2c') \
                .validate_config()

class TestGenBench3DConfig:

    def test_correct_genbench3d_file(self, pytest_load_config):
        GenBench3DConfig(pytest_load_config).validate_config()

    def test_incorrect_genbench_input_starts_with_num(self, pytest_load_config, genbench3d_key: ConfigGenBench3DKey):
        with pytest.raises(TypeError):
            GenBench3DConfig(pytest_load_config).update(genbench3d_key.input, '555').validate_config()