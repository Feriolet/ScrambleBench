import pytest
import pandas as pd
import pathlib
from pathlib import Path
import subprocess
import yaml

from ScrambleBench.script.config_preparation.config_input import InputConfig, check_input
from ScrambleBench.script.config_preparation.config_model import ModelConfig, check_model
from ScrambleBench.script.p01_generate_config_refactored import load_config
from ScrambleBench.script.error_handler import FileTypeError, FileDataError, DirNotFound

@pytest.fixture
def global_var():

    # INPUT CONFIG
    pytest.correct_protein_pdb = 'input_0/protein.pdb'
    pytest.correct_complex_pdb = 'input_0/complex.pdb'
    pytest.correct_ligand_sdf = 'input_0/ligand.sdf'

    pytest.complex_without_protein_pdb = 'input_0/complex_without_protein.pdb'
    pytest.incorrect_protein_sdf = 'input_0/incorrect_protein.sdf'
    pytest.incorrect_protein_pdb = 'input_0/incorrect_protein.pdb'
    
    pytest.incorrect_ligand_sdf = 'input_0/incorrect_ligand.sdf'
    pytest.multiple_ligand_sdf = 'input_0/multiple_ligand.sdf'
    pytest.translated_ligand_sdf = 'input_0/translated_ligand.sdf'

    pytest.empty_pdb = 'input_0/empty.pdb'
    pytest.empty_sdf = 'input_0/empty.sdf'
    
    pytest.config_fname = 'config.yml'
    pytest.config_input_name = 'input'
    pytest.config_input_complex_name = 'complex_path'
    pytest.config_input_pdb_name = 'pdb_path'
    pytest.config_input_sdf_name = 'sdf_path'
    pytest.config_input_title_name = 'protein_title'

    # MODEL CONFIG
    pytest.config_model_name = 'model'
    pytest.config_model_title_example_name = 'pocket2mol'
    pytest.config_model_title_name = 'name'
    pytest.config_model_dir_name = 'dir'
    pytest.config_model_conda_name = 'conda_env'

    pytest.config_generation_name = 'generation'
    pytest.config_post_generation_name = 'post_generation'

# def test_update_config_pass(global_var):
#     check_input(InputConfig(load_config(pytest.config_fname)).update(pytest.config_input_title_name, 'test'))

# def test_update_config_incorrect_type_error(global_var):
#     with pytest.raises(TypeError):
#         check_input(InputConfig(load_config(pytest.config_fname)).update(pytest.config_input_complex_name, 20.5))

# def test_load_config_incorrect_key_error(global_var):
#     config_data = load_config(pytest.config_fname)
#     config_data[pytest.config_input_name]['hiiiiiii'] = config_data[pytest.config_input_name].pop(pytest.config_input_complex_name)
#     print(config_data)
#     with pytest.raises(KeyError):
#         InputConfig(config_data)

# def test_update_config_incorrect_key_error(global_var):
#     with pytest.raises(TypeError):
#         check_input(InputConfig(load_config(pytest.config_fname)).update('hiiiii', 20.5))

# def test_input_file_pass(global_var):
#     check_input(InputConfig(load_config(pytest.config_fname)))

# def test_complex_pdb_not_found(global_var):
#     data = InputConfig(load_config(pytest.config_fname)).update(pytest.config_input_complex_name, 'xyz1230987')

#     with pytest.raises(FileNotFoundError):
#         check_input(data)

# def test_complex_pdb_wrong_format(global_var):
#     data = InputConfig(load_config(pytest.config_fname)).update(pytest.config_input_complex_name, pytest.incorrect_protein_sdf)

#     with pytest.raises(FileTypeError):
#         check_input(data)

# def test_complex_pdb_empty_ligand_error(global_var):
#     data = InputConfig(load_config(pytest.config_fname)).update(pytest.config_input_complex_name, pytest.correct_protein_pdb)

#     with pytest.raises(FileDataError):
#         check_input(data)

# def test_complex_pdb_empty_protein_error(global_var):
#     data = InputConfig(load_config(pytest.config_fname)).update(pytest.config_input_complex_name, pytest.complex_without_protein_pdb)

#     with pytest.raises(FileDataError):
#         check_input(data)

# def test_complex_pdb_empty_file_error(global_var):
#     data = InputConfig(load_config(pytest.config_fname)).update(pytest.config_input_complex_name, pytest.empty_pdb)

#     with pytest.raises(FileDataError):
#         check_input(data)


# def test_protein_pdb_not_found_error(global_var):
#     data = InputConfig(load_config(pytest.config_fname)).update(pytest.config_input_pdb_name, 'xyz1230987')

#     with pytest.raises(FileNotFoundError):
#         check_input(data)

# def test_protein_pdb_wrong_format_error(global_var):
#     data = InputConfig(load_config(pytest.config_fname)).update(pytest.config_input_pdb_name, pytest.incorrect_protein_sdf)

#     with pytest.raises(FileTypeError):
#         check_input(data)

# def test_protein_pdb_contains_ligand_error(global_var):
#     data = InputConfig(load_config(pytest.config_fname)).update(pytest.config_input_pdb_name, pytest.correct_complex_pdb)

#     with pytest.raises(FileDataError):
#         check_input(data)

# def test_protein_pdb_mismatch_complex_error(global_var):
#     data = InputConfig(load_config(pytest.config_fname)).update(pytest.config_input_pdb_name, pytest.incorrect_protein_pdb)

#     with pytest.raises(FileDataError):
#         check_input(data)

# def test_protein_pdb_empty_file_error(global_var):
#     data = InputConfig(load_config(pytest.config_fname)).update(pytest.config_input_pdb_name, pytest.empty_pdb)

#     with pytest.raises(FileDataError):
#         check_input(data)


# def test_ligand_sdf_not_found_error(global_var):
#     data = InputConfig(load_config(pytest.config_fname)).update(pytest.config_input_sdf_name, 'xyz1230987')

#     with pytest.raises(FileNotFoundError):
#         check_input(data)

# def test_ligand_sdf_wrong_format_error(global_var):
#     data = InputConfig(load_config(pytest.config_fname)).update(pytest.config_input_sdf_name, pytest.empty_pdb)

#     with pytest.raises(FileTypeError):
#         check_input(data)

# # def test_ligands_contains_protein_error(global_var):
# #     data = InputConfig(load_config(pytest.config_fname))).update(pytest.config_input_sdf_name, pytest.correct_complex_pdb)

# #     with pytest.raises(FileDataError):
# #         check_input(data)

# def test_ligand_sdf_empty_file_error(global_var):
#     data = InputConfig(load_config(pytest.config_fname)).update(pytest.config_input_sdf_name, pytest.empty_sdf)

#     with pytest.raises(FileDataError):
#         check_input(data)

# def test_ligands_contains_multiple_ligands_error(global_var):
#     data = InputConfig(load_config(pytest.config_fname)).update(pytest.config_input_sdf_name, pytest.multiple_ligand_sdf)

#     with pytest.raises(FileDataError):
#         check_input(data)

# def test_ligand_sdf_mismatch_complex_error(global_var):
#     data = InputConfig(load_config(pytest.config_fname)).update(pytest.config_input_sdf_name, pytest.incorrect_ligand_sdf)

#     with pytest.raises(FileDataError):
#         check_input(data)

# def test_ligand_sdf_not_in_pocket_error(global_var):
#     data = InputConfig(load_config(pytest.config_fname)).update(pytest.config_input_sdf_name, pytest.translated_ligand_sdf)

#     with pytest.raises(FileDataError):
#         check_input(data)


# def test_model_pass(global_var):
#     ModelConfig(load_config(pytest.config_fname))

# def test_modelcommercial_pass(global_var):
#     check_model(ModelConfig(load_config(pytest.config_fname)).update(pytest.config_model_title_example_name, pytest.config_model_dir_name, 'not_applicable').update(pytest.config_model_title_example_name, pytest.config_model_conda_name, 'not_applicable'))

# def test_load_config_incorrect_model_key_error(global_var):
#     config_data = load_config(pytest.config_fname)
#     config_data[pytest.config_model_name][pytest.config_model_title_example_name]['hiiiiiii'] = config_data[pytest.config_model_name][pytest.config_model_title_example_name].pop(pytest.config_model_title_name)
#     with pytest.raises(KeyError):
#         ModelConfig(config_data)

# def test_incorrect_update_model_title_not_as_str_int_or_float(global_var):
#     with pytest.raises(TypeError):
#         check_model(ModelConfig(load_config(pytest.config_fname)).update(pytest.config_model_title_example_name, pytest.config_model_title_name, [1,2,3]))

# def test_incorrect_update_model_dir(global_var):
#     with pytest.raises(DirNotFound):
#         check_model(ModelConfig(load_config(pytest.config_fname)).update(pytest.config_model_title_example_name, pytest.config_model_dir_name, 'xyz12345'))

# def test_incorrect_update_conda_dir(global_var):
#     with pytest.raises(ValueError):
#         check_model(ModelConfig(load_config(pytest.config_fname)).update(pytest.config_model_title_example_name, pytest.config_model_conda_name, 'xyz12345'))

def test_correct_generation(global_var):
    check_generation(GenerationConfig(load_config(pytest.config_fname)))