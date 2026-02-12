import yaml
import rdkit
from rdkit import Chem
from rdkit.Chem import rdMolTransforms
import os
from split_protein_ligand import split_pocket_ligand
import argparse
import pathlib
import subprocess

def check_complex(complex_pdb):
    # Create pocket file
    if not os.path.isfile(complex_pdb):
        raise TypeError('File does not exist')
    elif complex_pdb[-4:] != '.pdb':
        raise TypeError('protein is not in PDB format')

    return split_pocket_ligand(complex_pdb, cutoff=10)


def check_pdb(protein_pdb):
    # Create pocket file
    if not os.path.isfile(protein_pdb):
        raise TypeError('File does not exist')
    elif protein_pdb[-4:] != '.pdb':
        raise TypeError('protein is not in PDB format')
    

def check_ligand(ligand_sdf):
    # Find ligand center
    if not os.path.isfile(ligand_sdf):
        raise TypeError('File does not exist')
    elif ligand_sdf[-4:] != '.sdf':
        raise TypeError('Ligand is not in SDF format')
    mol = Chem.SDMolSupplier(ligand_sdf, removeHs=False)[0]

    if not mol:
        raise ValueError('error parsing mol')
    if mol.GetNumAtoms() < 2:
        raise ValueError('mol is probably empty')

    try:
        centroid = list(rdMolTransforms.ComputeCentroid(mol.GetConformer(0), ignoreHs=True))
        if not abs(centroid[-1]) > 0.0001:
            raise ValueError('Ligand needs to be in 3D!')

        # As per Pocket2Mol, they mentioned that the first character should be a whitespace
        return ' ' + ','.join([f"{coord:.3f}" for coord in centroid])
    
    except ValueError:
        print('No conformers detected. Ligand needs to be in 3D')


def check_models(model_dict):
    output_dict = {}
    for model, model_dir in model_dict.items():
        if os.path.isdir(model_dir):
            output_dict[model] = os.path.abspath(model_dir)
        else:
            raise ValueError('Directory does not exist')
    
    return output_dict


def check_conda_env(env_dict):
    output_dict = {}

    # Run 'conda env list' command and capture the output
    result = subprocess.run(
        "conda env list | awk '{ print $1}'",
        shell=True,
        capture_output=True,
        text=True,
        check=True
    )
    output = result.stdout.split('\n')

    for model, model_env in env_dict.items():
        if model_env in output:
            output_dict[model] = model_env
        else:
            raise ValueError(f'Conda environment {model_env} does not exist')
    
    return output_dict

#https://github.com/yaml/pyyaml/issues/127
class MyDumper(yaml.SafeDumper):
    # HACK: insert blank lines between top-level objects
    # inspired by https://stackoverflow.com/a/44284819/3786245
    def write_line_break(self, data=None):
        super().write_line_break(data)

        if len(self.indents) == 1:
            super().write_line_break()


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('--config_path', type=str, required=True, help='path to the config file')

    args = parser.parse_args()

    #config_file = '/opt/veincent/GenAI_manuscript/run_config/config_5HT2C_prepared.yml'

    with open(args.config_path, 'r') as config_f:
        config_data = yaml.safe_load(config_f)
    
    pocket_fname = check_complex(config_data['input']['complex_path'])
    pocket_center = check_ligand(config_data['input']['sdf_path'])

    check_pdb(config_data['input']['pdb_path'])

    output_dict: dict = {}
    output_dict['input'] = {}
    output_dict['input']['complex_path'] = os.path.abspath(config_data['input']['complex_path'])
    output_dict['input']['pdb_path'] = os.path.abspath(config_data['input']['pdb_path'])
    output_dict['input']['sdf_path'] = os.path.abspath(config_data['input']['sdf_path'])
    output_dict['input']['pocket_path'] = os.path.abspath(pocket_fname)

    output_dict['parameter'] = config_data['parameter']
    output_dict['parameter']['pocket_coord'] = pocket_center
    
    output_dict['models_dir'] = check_models(config_data['models_dir'])
    output_dict['conda_env'] = check_conda_env(config_data['conda_env'])

    output_dict['output'] = {}
    
    num_sample_l = config_data['parameter']['num_sample'].split(',')

    for num_sample in num_sample_l:
        output_dir = os.path.abspath(f"{config_data['output']['output_dir']}_{num_sample.strip()}")
        pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)

        output_dict['output']['output_dir'] = output_dir
        output_dict['parameter']['num_sample'] = int(num_sample)

        yaml_output_fname = f'{output_dir}/run_generative_ai.yaml'
        with open(yaml_output_fname, 'w') as yaml_f:
            yaml.dump(output_dict, yaml_f, Dumper=MyDumper, sort_keys=False)

