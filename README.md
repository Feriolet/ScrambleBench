# ScrambleBench v0.1.0

## Introduction

![Intro graphic](asset/GenAI_Abstract.png)
Hi! Welcome to ScrambleBench, A Workflow for Comparative Assessment of Structure-based *de novo* Generative Models. This repository contains the code used for our manuscript. Our v0.1.0 version streamlined most of the features used in our manuscript. Please look out for our version 0.2.0 soon!

**Note: v0.1.0 is not a stable version, since it is still updating its code even today (see latest commits). Hence, the codes that are cloned today may be different tomorrow. Please update the repository regularly until the version is stable.**

## Table of Contents
- [ScrambleBench v0.1.0](#scramblebench-v010)
  - [Introduction](#introduction)
  - [Table of Contents](#table-of-contents)
  - [Features](#features)
  - [Future feature for v0.1.0](#future-feature-for-v010)
  - [Installation](#installation)
    - [Method 1: Install through environment.yaml (Recommended)](#method-1-install-through-environmentyaml-recommended)
    - [Method 2: Install manually](#method-2-install-manually)
  - [Model Repository Installation](#model-repository-installation)
    - [PMDM](#pmdm)
    - [DiffSBDD](#diffsbdd)
    - [Pocket2Mol](#pocket2mol)
    - [PocketFlow](#pocketflow)
      - [Automatic (Recommended)](#automatic-recommended)
      - [Manual](#manual)
    - [Lingo3DMol](#lingo3dmol)
      - [Recommended](#recommended)
      - [Manually](#manually)
  - [Other Installation](#other-installation)
    - [Easydock](#easydock)
    - [Genbench3D](#genbench3d)
    - [Diversity](#diversity)
  - [Pre-trained Model Installation](#pre-trained-model-installation)
  - [Usage](#usage)
    - [Examples](#examples)
    - [Step-by step pipeline](#step-by-step-pipeline)
      - [1. Prepare Config File](#1-prepare-config-file)
      - [2. Run Generation](#2-run-generation)
      - [2.5. Parameter Key](#25-parameter-key)
      - [3. Prepare molecule](#3-prepare-molecule)
      - [4. Analysis Key](#4-analysis-key)
      - [4a. GenBench3D analysis](#4a-genbench3d-analysis)
      - [4b. Redocking analysis](#4b-redocking-analysis)
      - [4c. Diversity Analysis](#4c-diversity-analysis)
      - [4d. Pharmacophore-based screening](#4d-pharmacophore-based-screening)
      - [5. Compilation of data analysis](#5-compilation-of-data-analysis)
      - [6. Plotting](#6-plotting)
  - [Data Availability](#data-availability)
  - [Reproducing Figures](#reproducing-figures)
  - [FAQ](#faq)


## Features

- **Streamlined workflow**: from de novo molecule generation, molecule validation, and analysis to a single `.csv` file
- **Integrated software**: This repository integrates multiple docking program supported by `easydock`, conformation analysis by `genbench3d`, and various molecular diversity, such as conventional `tanimoto ecfp`, `hamdiv`, and `generic BM scaffold`.
- **Checkpoints**: Checkpoints are available for most of the worflow to resume process that got interrupted.
- **Batch runs**: Each YAML configs are generated for multiple possible parameters (`protein targets`, `num_sample`) which are suitables for batch runs in job schedulers.
- **Flexible parameters**: `ScrambleBench` can be run with or without protein inputs
- **Multiple program model**: `ScrambleBench` supports open-source program (`easydock`) and commercial program (`Schrödinger`).


## Future feature for v0.1.0

- **Report summary**: Generate a PDF report for each config YAML detailing the results of each analysis of docking, diversity, and conformational analysis.
- **Additional report**: JSON file that reports on `validity3d`, `uniqueness`, and `diversity`
- **Pharmacophore Screening**: implements pharmacophore screening to reflect ideal target binding


## Installation

```sh
git clone https://github.com/Feriolet/ScrambleBench
cd ScrambleBench
```

### Method 1: Install through environment.yaml (Recommended)
```bash
conda env create -f environment.yaml
conda activate scramblebench
pip install -e .
```

### Method 2: Install manually 
```bash
conda create -n scramblebench python=3.10 oddt openbabel
conda activate scramblebench
pip install rdkit numpy matplotlib ptitprince seaborn pandas meeko fastparquet pyarrow uv pytest openmm pdbfixer pyyaml Bio pytest-xdist
pip install -e .
```

## Model Repository Installation

You can write the repository anywhere, but I put mine in the dedicated `models` folder

```bash
mkdir models
cd models
```

### PMDM

We are not installing the qvina environment for our study.
The mol.yml is changed to included --extra-link-url and --find-links for pytorch
In case torch can not be imported, may need to delete the conda/env/pmdm/lib/libcudart.so.11.0 i if libcudart.so.11.8.*.* is also present in the conda/env/pmdm/lib.

```bash
git clone https://github.com/Feriolet/PMDM
cd PMDM
conda env create -f mol.yml

rm /path/to/miniforge/envs/benchmark_pmdm/lib/libcudart.so.11.0
```

### DiffSBDD
```bash
git clone https://github.com/arneschneuing/DiffSBDD
cd DiffSBDD
conda env create -f environment.yaml
mkdir checkpoints
conda activate diffsbdd
pip install "setuptools<11.0"
```

### Pocket2Mol

The config of `sample_for_pdb.yaml` is changed to allow Pocket2Mol to expand its atom
```bash
git clone https://github.com/Feriolet/Pocket2Mol
cd Pocket2Mol
conda env create -f env_cuda113.yml
```

### PocketFlow

#### Automatic (Recommended)
```bash
git clone https://github.com/Feriolet/PocketFlow
cd PocketFlow
conda env create -f environment.yml
```

#### Manual
```bash
conda create -n benchmark_pocketflow python=3.10 pymol-open-source=2.5.0 openbabel -y
conda activate benchmark_pocketflow
pip install torch==1.13.0+cu117 torchvision==0.14.0+cu117 torchaudio==0.13.0 --extra-index-url https://download.pytorch.org/whl/cu117
pip install scipy numpy==1.23.0
pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-1.13.0+cu117.html --no-index
pip install rdkit tensorboard six lmdb easydict
pip install torch_geometric==2.3.1
```

### Lingo3DMol

Note: Lingo3DMol/util/pocket_code_all.py was modified to increase pocket representation length from 500 to 1000, probably because the pocket size is too big

#### Recommended
```bash
git clone https://github.com/Feriolet/Lingo3DMol
cd Lingo3DMol
conda env create -f environment.yml
```

#### Manually
```bash
conda create -n lingo3dmol python=3.8
conda activate lingo3dmol
conda install pytorch==1.10.1 torchvision==0.11.2 torchaudio==0.10.1 cudatoolkit=11.3 -c pytorch -c conda-forge
pip install scipy==1.7.3 pandas==1.5.1 numpy==1.20.3 rdkit==2022.09.1 psutil torch_geometric==2.3.1
#maybe add conda install pyg -c pyg if still does not work
```

## Other Installation

### Easydock

```bash
git clone https://github.com/ci-lab-cz/easydock
cd easydock
conda env create -f env.yml -n benchmark_easydock
```

### Genbench3D
```bash
git clone https://github.com/Feriolet/genbench3d.git
cd genbench3d
conda env create -f environment.yml

#also install adfr if you can
```

Please add this kernel model to the genbench folder

```bash

cp LigBoundConf* [genbench3d root dir]
cd [genbench3d root dir]
tar -xvzf LigBoundConf_geometry_kernel_densities.tar.gz
tar -xvzf LigBoundConf_geometry_values.tar.gz
```

### Diversity
```bash
git clone https://github.com/HXYfighter/HamDiv
cd HamDiv
conda create -n python_tsp numpy rdkit tqdm networkx -y
conda activate python_tsp
pip install python_tsp
```

## Pre-trained Model Installation

As I do not own the model, you can download the corresponding model in the respective owner's Github repository:


Pocket2Mol (to `Pocket2Mol/ckpt/`): [pretrained_Pocket2Mol.pt](https://drive.google.com/drive/folders/1KfdOczjUPITPhIvCuBmnj4xFTV-iI2xB)

PocketFlow (to `PocketFlow/ckpt/`): [ZINC-pretrained-255000.pt](https://github.com/Saoge123/PocketFlow)

Lingo3DMol (to `Lingo3DMol/checkpoint/`): [contact.pkl](https://stonewise-lingo3dmol-public.s3.cn-northwest-1.amazonaws.com.cn/contact.pkl) and [gen_mol.pkl](https://stonewise-lingo3dmol-public.s3.cn-northwest-1.amazonaws.com.cn/gen_mol.pkl)

DiffSBDD (to `DiffSBDD/checkpoints`): [crossdocked_fullatom_cond.ckpt](https://zenodo.org/record/8183747/files/crossdocked_fullatom_cond.ckpt?download=1)

PMDM (to `PMDM/`): [500.pt](https://zenodo.org/records/10630921)

## Usage

### Examples

You can access the `example` folder to run multiple kinds of benchmark. We encourage you to use the `example/run_multiple_targets_multiple_parameters` since this is the folder that is most often tested. Alternatively, you can use
`run_simplest_config` for very fast generation and analysis.

```bash
cd example/run_simplest_config
chmod +x ./prepare_config.sh
chmod +x ./run_scramblebench.sh

conda activate scramblebench
./run_scramblebench.sh
```

### Step-by step pipeline

#### 1. Prepare Config File

Before running the molecule generation, it is necessary to understand the appropriate input for ScrambleBench. Several examples are available in the `example` folder. Although  `p1_generate_config.py` scans all of the `yaml` keys and subkeys, we will reveal the config keys gradually throughout the documentation :).

```yaml
input:
  protein1:
    complex_path: input/5ht2c/gpcr_5ht2c_6bqh_complex_autoprepared.pdb
    pdb_path: input/5ht2c/gpcr_5ht2c_6bqh_protein_autoprepared.pdb
    sdf_path: input/5ht2c/gpcr_5ht2c_6bqh_ligand.sdf

input_dir:
  dirpath: input
```

**Input Key**

string: `input`

This key must be followed by protein name(s) (e.g., `protein1`) **MULTIPLE proteins are allowed**. Each protein name can have the following fields: 
| Field        | dtype                    | description                                                           | required | default       |
|--------------|--------------------------|-----------------------------------------------------------------------|----------|---------------|
| complex_path | str                      | directory path to protein-ligand complex pdb                          | True     | N/A           |
| pdb_path     | str                      | directory path to protein pdb                                         | True     | N/A           |
| sdf_path     | str                      | directory path to ligand sdf                                          | True     | N/A           |
| name         | str                      | protein name                                                          | True     | subheading    |
| pocket_path  | str or None              | directory path to pocket pdb                                          | False*   | autogenerated |
| pocket_coord | str or list[str] or None | pocket x,y,z center coordinate (must be comma separated or as a list) | False*   | autogenerated |

*: not required for `p1_generate_config.py` but required for other script.

string: `input_dir`
| Field   | dtype | description                          | required | default |
|---------|-------|--------------------------------------|----------|---------|
| dirpath | str   | directory path to protein folder pdb | True     | None    |

Currently, all protein target should have a 1) complex pdb, 2) protein pdb, and 3) ligand sdf, because different AI models require different inputs. If you use `input_dir` key, there **MUST** be a **SINGLE** file of protein.pdb, complex.pdb, and ligand.sdf with the string `protein`, `complex`, and `ligand`, respectively, so that the script can detect which one is which. For example, the script will detect `gpcr_ligand.sdf` as ligand file, but not `gpcr.sdf`, because there is no `ligand` string in the file.

The path directory should be the following:

```txt
input_folder
├── protein1
│   ├── complex.pdb
│   ├── protein.pdb
│   └── ligand.sdf
└── protein2
    ├── complex.pdb
    ├── protein.pdb
    └── ligand.sdf
```

In the config file, the input key should be either `input` or `input_dir`. Config file with `input_dir` key must be run with `--dirpath_input` argument in the `p1_generate_config.py` script, which can be seen in the `example/run_multiple_targets_single_dir` folder. 

```txt
usage: p1_generate_config.py [-h] [-i INPUT] [-o OUTPUT] [--dirpath_input]

Prepare config file for ScrambleBench

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        config yaml input file
  -o OUTPUT, --output OUTPUT
                        config yaml output file prefix (i.e., not path directory)
  --dirpath_input       write input key as a single directory
```
Multiple values for the input protein names and parameters (i.e., box_size and num_sample) are supported, which can be seen in the `example` folder.



In the future, we will support complex pdb only, but the file directory should still be the same (i.e., individual protein target has its own folder).

A config file and txt file containing a list of the config.yaml will be written in the folder list in value for `config[generation][input]` key.

For best practice, please write the absolute pathdir, as relative pathdir is relative to the executable path (i.e., pwd). Despite this, relative pathdir is still allowed, and ScrambleBench will resolve it in the config output.

In the output file, you can notice that the protein name in the subheading will be moved to the `name` field of the `input` key.

```yaml
input:
  complex_path: /home/Veincent/manuscript/ScrambleBench/example/run_multiple_targets_multiple_parameters/input/5ht2c/gpcr_5ht2c_6bqh_complex_autoprepared.pdb
  pdb_path: /home/Veincent/manuscript/ScrambleBench/example/run_multiple_targets_multiple_parameters/input/5ht2c/gpcr_5ht2c_6bqh_protein_autoprepared.pdb
  sdf_path: /home/Veincent/manuscript/ScrambleBench/example/run_multiple_targets_multiple_parameters/input/5ht2c/gpcr_5ht2c_6bqh_ligand.sdf
  pocket_path: /home/Veincent/manuscript/ScrambleBench/example/run_multiple_targets_multiple_parameters/input/5ht2c/gpcr_5ht2c_6bqh_complex_autopreparedcut10/gpcr_5ht2c_6bqh_complex_autopreparedcut10_pocket_withH.pdb
  pocket_coord: 38.87,31.02,56.85
  name: protein1
```

#### 2. Run Generation

```yaml
#yaml config
model: 
  diffsbdd:
    name: DiffSBDD
    dir: models/DiffSBDD
    conda_env: benchmark_diffsbdd
  pocket2mol:
    name: Pocket2Mol
    dir: models/Pocket2Mol
    conda_env: benchmark_pocket2mol

generation:
  input: output
  output: AI_Generation
  parameter:
    box_size: 16
    num_sample: 50,100 
    name: 5ht2c_prepared_0
```

**Model Key**

string: `model`

This key must be followed by a model name (e.g., `pocket2mol`, preferably lowercase) **MULTIPLE models are allowed**. Each model name can have the following fields:
| Field | dtype       | description                                           | required | default |
|-------|-------------|-------------------------------------------------------|----------|---------|
| name  | str         | model name used for output and downstream file naming | True     | N/A     |
| dir   | str or None | directory path to installed model folder              | False    | None    |
| conda | str or None | conda environment name used to run the model          | False    | None    |

If the `dir` and `conda` field is empty or filled with `non_applicable`, the `p2_execute_generation.py` won't run these models. This is reserved for models that have generated the molecules without using `ScrambleBench`.

Note that there are differences in the `name` field and the model name (e.g., `pocket2mol`). The `name` field is essential for naming the intermediates and output files for downstream analysis, while the `pocket2mol` is essential for the `script/utils/generation_utils/generation_template.sh`

**Generation Key**

string: `generation`

| Field           | dtype       | description                                        | required | default                                                |
|-----------------|-------------|----------------------------------------------------|----------|--------------------------------------------------------|
| input           | str         | typically directory path to the prepared yaml file | True     | N/A                                                    |
| output          | str         | directory path to output generation                | True     | N/A                                                    |
| script_pathfile | str or None | file path to bash script to execute generation     | False*   | `src/scramblebench/script/utils/generation_utils/generation_template.sh` |
| parameter       | dict        | model parameter used for generation                | False*   | see below                                              |

*: not required for `p1_generate_config.py` but required for other script.

**Generation Parameter Key**

string: `parameter`

| Field      | dtype                | description                                                                                            | required | default       |
|------------|----------------------|--------------------------------------------------------------------------------------------------------|----------|---------------|
| box_size   | None or float or str | box size for molecule generation (similar to docking box) in Amstrong (must be comma separated if str) | False*    | 16            |
| num_sample | None or int or str | number of requested generated molecule (must be comma separated if str)                                | False*    | 100           |
| name       | str or None          | job name                                                                                               | False    | generic_title |

*: not required for `p1_generate_config.py` but required for other script.

The generation script allows two kinds of input: a single yaml file input or a txt file containing a list of yaml pathdir. 

Before generation, this script will check whether the generation script (default: `src/script/utils/generation_template.sh`) will generate molecules by the correct AI model by their conda environment (e.g., pmdm model should have a `$model_pmdm_conda_env` and pocket2mol model should have a `$model_pocket2mol_conda_env` string in the bash script).

```txt
usage: p2_execute_generation.py [-h] -i INPUT

Run de novo molecule generation after p1_generate_config.py

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        config yaml input file or txt file containing yaml filepath
```

Expected output:
```
output_folder
├── model1
│   └── generated_model1_ligand.sdf
└── model2
│   └── generated_model2_ligand.sdf
└── summary
    ├── generated_model1_ligand.sdf
    └── generated_model2_ligand.sdf   
```


  batch_parameter: {}

#### 2.5. Parameter Key

For further downstream analysis, this key is essential as all parts will require this information one way or another. This parameter key is auto-generated from `p1_generate_config.py` and hence does not need to be prefilled.


**Parameter Key**

string: `parameter`

| Field        | dtype     | description                         | required | default |
|--------------|-----------|-------------------------------------|----------|---------|
| protein_name | str       | identical as the `input` key        | True     | N/A     |
| model_list   | list[str] | list of model `name` in `model` key | True     | N/A     |
| num_sample   | int       | identical as the `generation` key   | True     | 100     |
| batch_parameter   | dict       | `parameter` in the `generation` key that have multiple values (i.e., `num_sample` and `box_size`)  | False     | {}     |


#### 3. Prepare molecule

After generation, the molecules need to be prepared and validated.

The input **MUST** contains the `summary` folder with an `.sdf` file containing all generated ligands and named with the model name, as seen in the previous step.

Expected input:
```
summary
  ├── generated_model1_ligand.sdf
  └── generated_model2_ligand.sdf  
```
```yaml
#yaml config
post_generation:
  input: AI_Generation
  output: cheminformatics_input_prepared
  pick_last: Pocket2Mol
  pick_random: null
```

**Post Generation Key**

string: `post_generation`

| Field       | dtype       | description                                                                                | required | default |
|-------------|-------------|--------------------------------------------------------------------------------------------|----------|---------|
| input       | str         | input directory path to generated molecule                                                 | True     | N/A     |
| output      | str         | output directory path to processed generated molecule                                      | True     | N/A     |
| pick_last   | None or str | filled with model name to pick last `num_sample` molecule (must be comma separated if str) | False    | None    |
| pick_random | None or str | filled with model name to randomly pick molecule (must be comma separated if str)          | False    | None    |


The key `pick_last` and `pick_random` should be filled with the model name (case sensitive). This will only trigger when the model generate more ligand than the requested `num_sample`. For example, if the requested num_sample is `100`, but the model has generated `200` ligand, the script will trim it down to `100` ligands instead. `pick_last` will pick the last `100` ligands from the `.sdf` file, while `pick_random` will randomly pick ligands. If the model is unspecified, the script will use `pick_random` by default.


```txt
usage: p3_prepare_molecule.py [-h] -i INPUT

Prepare generated molecules for downstream analysis

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        config yaml input file or txt file containing yaml filepath
```

Expected_output:
```
output_folder
├── model1
│   └── prepared_model1_ligand.sdf
└── model2
    └── prepared_model2_ligand.sdf
```

#### 4. Analysis Key

This key is optional if the user wishes to perform the analysis of their generated molecules. In the config files, the key `analysis` must be added. The subkeys `genbench3d`, `redocking`, and `diversity` are supported.

#### 4a. GenBench3D analysis

This will run the GenBench3D analysis. Please fill in the config files in the GenBench3D repository in `config/GenAI_evaluation.yaml`. An example is also shown in the `example/run_multiple_targets_multiple_parameters/config/GenAI_evaluation.yaml`.

Expected input:
```
input_folder
├── model1
│   └── prepared_model1_ligand.sdf
└── model2
    └── prepared_model2_ligand.sdf 
```

Scramblebench script currently does not attempt to access genbench3d config, which is why there is a duplicate keys in both the `ScrambleBench` and `GenBench3D` configs (e.g., `schrodinger_dir`).

```yaml
# ScrambleBench config
analysis:
  genbench3d:
    input: cheminformatics_input_prepared
    output: cheminformatics_analysis
    genbench3d_dir: /home/Veincent/manuscript/ScrambleBench/models/genbench3d
    conda_env: benchmark_genbench3d
    schrodinger_dir: /opt/schrodinger2025-4
    genbench3d_config: /home/Veincent/manuscript/ScrambleBench/models/genbench3d/GenAI_evaluation.yaml
    do_complex_forcefield_minimisation: True
    do_docking_forcefield_minimisation: True
    skip_genbench3d_protonation: True
```


**GenBench3D Key**

string: `genbench3d`

| Field                              | dtype        | description                                                            | required | default |
|------------------------------------|--------------|------------------------------------------------------------------------|----------|---------|
| input                              | str          | input directory path to prepared generated molecule                    | True     | N/A     |
| output                             | str          | output directory path to genbench3d analysis                           | True     | N/A     |
| genbench3d_dir                     | str          | directory path to installed genbench3d folder                          | True     | N/A     |
| conda_env                          | str          | conda environment name to execute genbench3d script                    | True     | N/A     |
| genbench3d_config                  | str          | file path to the genbench3d_config (see below)                         | True     | N/A     |
| schrodinger_dir                    | None or str  | directory path to the schrodinger root dir                             | False    | None    |
| do_complex_forcefield_minimisation | None or bool | whether to perform MMFF94 minimisation before analysis                 | False    | None    |
| do_docking_forcefield_minimisation | None or bool | whether to perform mininplace docking                                  | False    | None    |
| skip_genbench3d_protonation        | None or bool | whether to skip both ligand and protein protonation by adfr and/obabel | False    | None    |



```yaml
# genbench3d config
benchmark_dirpath: GenBench3D rootdir repository
glide_working_dir: I am not sure about this
pocket_distance_from_ligand: 5.0 # Angstrom

bin:
  prepare_receptor_bin_path: 'yourpath/models/genbench3d/ADFRsuite_x86_64Linux_1.0/bin/prepare_receptor' # https://ccsb.scripps.edu/adfr/downloads/
  glide_path: '/opt/schrodinger2025-4/glide'
  structconvert_path: '/opt/schrodinger2025-4/utilities/structconvert'

data:
  ligboundconf_path: 'yourpath/models/genbench3d/S2_LigBoundConf_minimized.sdf' # https://pubs.acs.org/doi/suppl/10.1021/acs.jcim.0c01197/suppl_file/ci0c01197_si_002.zip

genbench3d:
  minimum_pattern_values: 50
  tfd_threshold: 0.2
  q_value_threshold: 0.001
  steric_clash_safety_ratio: 0.75
  maximum_ring_plane_distance: 0.1 # Angstrom
  consider_hydrogens: False
  include_torsions_in_validity3D: False
  add_minimized_docking_scores: True
  overwrite_results: True

minimization:
  distance_constraint: 1.0 # Angstrom
  n_steps: 1000

vina:
  scoring_function: 'vina'
  size_border: 35 # Angstrom
  n_cpus: 100
  seed: 2023
```

```txt
usage: p4_analyse_genbench3d.py [-h] -i INPUT

Perform genbench3d analysis after molecule preparation

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        config yaml input file or txt file containing yaml filepath
```

Because genbench3d allows for many different parameter choices, each analysis can take a long time. Below is some of the parameter that can be used:

| Parameter            | Description                                       | GenBench3D Config Field                     | Default                                                                                                |
|----------------------|---------------------------------------------------|---------------------------------------------|--------------------------------------------------------------------------------------------------------|
| Complex Minimisation | Whether to do MMFF94 minimisation before analysis | `do_complex_forcefield_minimisation`        | False (setting True will perform both analysis)                                                        |
| Docking Minimisation | Whether to do mininplace scoring                  | `do_docking_forcefield_minimisation`        | False (setting True will perform both analysis)                                                        |
| Docking Program      | Which docking program to use                      | varied (`--schrodinger_dir` to do Glide SP) | Vina (unless input key in `ScrambleBench` config is missing; see `examples/run_without_protein_input`) |

In order to track GenBench3D analysis, we have set up a `genbench3d_checkpoint.json` to prevent repeat analysis in case analysis terminated halfway.

Expected_output:
```
output_folder
├── structural_input
│   ├── protein.pdb (from input)
│   ├── ligand.sdf (from input)
│   ├── protein.pdbqt (for vina)
│   ├── protein_grid.txt (for glide)
│   └── other intermediates files
├── json_output
│   ├── prepared_model1_ligand_minimisation.json
│   ├── prepared_model1_ligand_minimisation.log
│   ├── prepared_model2_ligand_minimisation.json
│   └── prepared_model2_ligand_minimisation.log
├── Vina
├── Glide
└── genbench3d_checkpoint.json # checks if analysis is done per model per minimisation 
```

Sample of `genbench3d_checkpoint.json`:
```json
{
    "Pocket2Mol": {
        "unminimised": "COMPLETED",
        "minimised": "FAILED"
    },
    "DiffSBDD": {
        "unminimised": "FAILED",
        "minimised": "PENDING"
    }
}
```
#### 4b. Redocking analysis

For redocking, we currently support `easydock` redocking and `Glide SP` redocking. However, to run easydock, it is necessary to have some knowledge of how to run easydock to edit our configuration (mainly the easydock config and grid file). For now, my script supports vina by default, but will be able to support other docking program supported by easydock in the future. However, this is not tested yet!

Separate installation is needed to run these programs. Please refer to the easydock documentation here: [protonation](https://easydock.readthedocs.io/en/latest/usage/#protonation-options) and [docking](https://easydock.readthedocs.io/en/latest/usage/#molecular-docking)

string: `redocking`

supported subkeys: `protonation`, `docking`


Expected input:
```
input_folder
├── model1
│   └── prepared_model1_ligand.sdf
└── model2
    └── prepared_model2_ligand.sdf 
```

```yaml
analysis:
  redocking:
    protonation: 
      method: ligprep
      input: cheminformatics_input_prepared
      output: docking
      env: /opt/schrodinger2025-4

    docking:
      easydock:
        input: docking
        output: vina_output
        conda_env: benchmark_easydock
        protein_pdbqt_preparation: adfr
        protein_pdbqt_executable: /opt/veincent/software/ADFRsuite-1.0/ADFRsuite_x86_64Linux_1.0/bin/prepare_receptor 
        protein_preparation: /opt/schrodinger2025-4
        docking_program: vina
        protonation: null
        config_fname: config/easydock_config_5HT2C.yml

      glide:
        input: cheminformatics_input_prepared
        output: glide_output
        schrodinger_dir: /opt/schrodinger2025-4
        reward_intra_hbonds: null
        protonation: ligprep
        protein_preparation: protwizard
```

**Redocking Protonation Key**

string: `protonation`

| Field  | dtype | description                                                                                                       | required | default |
|--------|-------|-------------------------------------------------------------------------------------------------------------------|----------|---------|
| input  | str   | input directory path to prepared generated molecule                                                               | True     | N/A     |
| output | str   | output directory path to ligand protonation                                                                       | True     | N/A     |
| method | str   | protonation method (choose from [`molgpka`, `unipka`, `chemaxon`] via `easydock`* or `ligprep` via `schrodinger`) | True     | N/A     |
| env    | str   | conda environment name to execute `easydock`  or  directory path to the `schrodinger` root dir                    | True     | N/A     |

*: Note that these protonation method must be previously installed before running.

**Redocking Docking Easydock Key**

string: `easydock`

| Field                     | dtype       | description                                                                                                                       | required | default                                                                          |
|---------------------------|-------------|-----------------------------------------------------------------------------------------------------------------------------------|----------|----------------------------------------------------------------------------------|
| input                     | str         | input directory path to prepared generated molecule                                                                               | True     | N/A                                                                              |
| output                    | str         | output directory path to easydock docking                                                                                         | True     | N/A                                                                              |
| conda_env                 | str         | conda environment name to execute easydock script                                                                                 | True     | N/A                                                                              |
| config_fname              | None or str | path file to easydock config                                                                                                      | False*   | `src/scramblebench/script/utils/docking_utils/easydock_config.yml` for vina only |
| protein_pdbqt_preparation | str or None | method to generation pdbqt file (choose from [`obabel` or `adfr`])                                                                | False    | obabel                                                                           |
| protein_pdbqt_executable  | str         | file path to `prepare_receptor`, otherwise just `obabel`                                                                          | False    | obabel                                                                           |
| protein_preparation       | str or None | method to prepare protein (choose from [`pdbfixer`, `obabel`] or write directory path to the `schrodinger` root dir )             | False    | obabel                                                                           |
| docking_program           | str or None | method to dock ligand (choose from [`vina`, `gnina`, `smina`, `vina-gpu`, `qvina`, `server`])                                     | False    | vina                                                                             |
| protonation               | None or str | method to protonate ligand (choose from [`molgpka`, `unipka`, `chemaxon`] or write directory path to the `schrodinger` root dir ) | False    | None                                                                             |
| ncpu                      | None or int | parallel cpu to run easydock                                                                                                      | False    | None                                                                             |

*: not required for vina only

**Redocking Docking Glide Key**

string: `glide`

| Field               | dtype        | description                                                  | required | default |
|---------------------|--------------|--------------------------------------------------------------|----------|---------|
| input               | str          | input directory path to prepared generated molecule          | True     | N/A     |
| output              | str          | output directory path to glide SP docking                    | True     | N/A     |
| schrodinger_dir     | str          | directory path to the `schrodinger` root dir                 | True     | N/A     |
| reward_intra_hbonds | bool or None | whether to reward intramolecular hydrogen bonds              | False    | None    |
| protonation         | str or None  | method to protonate ligand (choose from `ligprep` or None)   | False    | None    |
| protein_preparation | None or str  | method to prepare protein (choose from `protwizard` or None) | False    | None    |



```txt
usage: p4_analyse_redocking.py [-h] -i INPUT

Redocking analysis after molecule preparation

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        config yaml input file or txt file containing yaml filepath
```

Expected output (identical structure for `protonation`, `easydock`, and `glide` redocking):
```
output_folder
├── model1
│   └── docked_model1_ligand.sdf
└── model2
    └── docked_model2_ligand.sdf 
```

#### 4c. Diversity Analysis

This pipeline outlines the codes for plotting the diversity. Shout out to the original Github and the paper discussing about this

github: https://github.com/HXYfighter/HamDiv
paper: https://jcheminf.biomedcentral.com/articles/10.1186/s13321-024-00883-4

**MOST IMPORTANTLY, USE RDKIT V2025.9.3 OR LATER, BECAUSE THE RASCAL MCES HAS SOME BUGS IN EARLIER VERSIONS**

```yaml
analysis:
  diversity:
    input: cheminformatics_input_prepared
    output: diversity
    conda_env: python_tsp
    method:
    - distance: ecfp
      diversity: hamdiv
    - distance: mces
      diversity: hamdiv
    - distance: null
      diversity: generic_bm
```

| Field     | dtype      | description                                                           | required | default                                       |
|-----------|------------|-----------------------------------------------------------------------|----------|-----------------------------------------------|
| input     | str        | input directory path to prepared generated molecule                   | True     | N/A                                           |
| output    | str        | output directory path to diversity analysis                           | True     | N/A                                           |
| conda_env | str        | conda environment name to execute diversity python_tsp script         | True     | N/A                                           |
| method    | list[dict] | method to measure diversity of a ligand set (choose from table below) | False    | [{`distance`: `ecfp`, `diversity`: `hamdiv`}] |


Supported `distance` and `diversity` combination (please refer to HamDiv paper for definition for the distance)

**Do not attempt to do HamDiv MCES for > 500 ligands unless you have the huge RAM resources!**

| distance | diversity  |
|----------|------------|
| ecfp     | hamdiv     |
| mces     | hamdiv     |
| ecfp     | average    |
| null     | richness   |
| null     | rs         |
| null     | fg         |
| null     | bm         |
| null     | generic_bm |
| null     | intdiv     |
| null     | sumdiv     |
| null     | diam       |
| null     | sumdiam    |
| null     | sumbot     |
| null     | bot        |
| null     | dpp        |

```
usage: p4_analyse_diversity.py [-h] -i INPUT

Calculate Diversity of generated compound

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        config yaml input file or txt file containing yaml filepath
```

Expected output:
```
output_folder
├── model1
│   └── diversity_output.json
└── model2
    └── diversity_output.json 
```

#### 4d. Pharmacophore-based screening

In our manuscript, we did our pharmacophore-based screening using Schrödinger Phase. Unfortunately, we currently do not have an open-source pipeline for this. However, feel free to explore Easydock PLIF option which I might use to integrate this in v0.2.0

#### 5. Compilation of data analysis

After performing the data analysis, we can collect all of the information into one single `.csv` file.

The script will detect whether each data analysis key exists in the `ScrambleBench` config file (i.e., `genbench3d`, `diversity`, `redocking`)
```
usage: p5_collect_data_analysis.py [-h] -i INPUT

Collect data of analysis done in p4_analyse.py

options:
  -h, --help            show this help message and exit
  -i INPUT, --input INPUT
                        config yaml input file or txt file containing yaml filepath
```

Expected output:
```
folder
├── yaml_list.txt # from step 1
├── protein1 folder
└── all.csv # output
```

Sample output:
```csv
,protein_name,mol_id,Model,FF_unminimised_Minimized Vina score,FF_unminimised_Vina score,FF_unminimised_Minimized Glide score,FF_unminimised_Glide score,FF_minimised_Minimized Vina score,FF_minimised_Vina score,FF_minimised_Minimized Glide score,FF_minimised_Glide score,easydock_redocking_rmsd,easydock_redocking_score,glide_redocking_rmsd,glide_redocking_score
0,GPCR_5HT2C,Pocket2Mol_0,Pocket2Mol,-2.279,-2.071,-2.124,-2.014,-2.275,-1.777,-2.119,-1.977,5.032224972117204,-3.879,5.123656916811996,-3.140124946810722
1,GPCR_5HT2C,Pocket2Mol_1,Pocket2Mol,-2.183,-2.032,-4.233,-3.825,-2.07,-1.86,-4.259,-3.819,5.54947453926946,-3.568,5.92211499618567,-4.6347445657091395
2,GPCR_5HT2C,Pocket2Mol_10,Pocket2Mol,-2.423,-2.232,-4.403,-4.07,-2.349,-2.271,-4.37,-4.213,4.195787247347035,-3.435,5.830167941496112,-4.4475081845864093
```

#### 6. Plotting

To be supported in future release.


## Data Availability
Please install the files in our [Zenodo](https://zenodo.org/records/18503149) for reproducibility.

## Reproducing Figures

Please refer to the `v0.0.1` for the codes to reproduce the figure.

The main Figures of the manuscript can be reproduced by using the `07_plot_summary.py` code with the file `output_scramblebench_data_warehouse/data_warehouse.parquet` in the Zenodo file

## FAQ

**Q: I want to use `ScrambleBench` to evaluate other Gen AI models.**

**A:** There is multiple ways to do it. A simple way is to generate the molecules outside of ScrambleBench and evaluate it similar to `examples/run_without_protein_input`. However, this limits docking analysis since docking requires protein input.

Another way to do it is to integrate the generation of molecules within the `ScrambleBench`. To do this, simply edit the bash file at `src/scramblebench/script/utils/generation_utils/generation_template.sh` and add the necessary models generation command line here. Remember to include the model names into the generated files and add them into the `summary` folder. Lastly, add the model information in the `model` config key.

If you insist on generating the molecule without using `ScrambleBench` and wanted to include docking analysis, let me know and I will add another feature for this.

**Q: I want to generate molecules using `ScrambleBench`, but I don't want to separate the protein-ligand complex to protein pdb and ligand sdf.**

**A:** This feature is likely to be added in v0.2.0 or later, but I am not sure how to efficiently separate ligand and protein. Using RDKit may not be favorable, so I would likely use oddt to extract ligand and filter **ATOM** for protein instead. Let me know if there is a much more efficient way to separate these if you wish to incorporate this quicker.

**Q: I want to add other docking/conformation/diversity program into `ScrambleBench`. How do I do that?**

**A:** Unfortunately, there is no easy way to do this. After all, the reason that I wrote dedicated `KeyConfig` in the `src/scramblebench/script/config_preparation` for every key is that I can certainly filter what is printed out for configs. Simply adding more keys will not work, because these will be filtered out.

If you are using it for your own use, it is possible to just manually edit the source code, especially below the `if __name__ == "__main__"` line of code for the `px_script.py`. However, if you wish to make a PR to integrate it, it is best to learn how the `KeyConfig` class works. In my opinion, it is fairly easy to understand them once you are familiar with one or two of them, because they are written in a very template-y manner.

For docking program, another way is to integrate your desired docking program into `easydock`, which have a dedicated page for integrating their program viewed here https://easydock.readthedocs.io/en/latest/custom_docking/.

I may release a documentation dedicated for developers for this in the future.

**Q: I want to contribute to the `ScrambleBench`.**

**A:** We are happy to have someone who are interested in developing our repository. However, since we have not added the necessary tests for our repository, it is very difficult to determine if the upcoming PR will break the codes or not. Feel free to make PR whenever you wish, but the PR will probably take a long time to review. As the standard practice, please do not upload very huge chunk of line of codes to the point where it is impossible to review hahaha.

**Q: Do I need to prepare my protein-complex and/or protein input for `ScrambleBench`?**

**A:** It depends on the purpose of using `ScrambleBench`. In my experience, if you want to generate the molecules, it is not really necessary to prepare them, because these models tend to remove the hydrogens in protein and ligand, possibly because of the encoding aspect. However, it may be a good idea to protonate them early on for analysis (e.g., `docking` and `genbench3d`) so you can save some time. Regardless, there are separate protonation/preparation in the current implementations.

**Q: I tried to install the conda environment of `ScrambleBench`/`Gen AI`/`third-party` software and it does not work.**

**A:** Unfortunately, this is an inevitable phenomenon where some dependencies are no longer supported. Please raise this issue on the `Issues` section at this repository or the repository of the model or software.

**Q: Will `ScrambleBench` be continuously improved even after its publication (if it gets published)?**

**A:** This is a very difficult question to answer. While I want to integrate more tools into the repository like `PoseBusters`, `HEAD-TED`, and other diversity tools, it depends on the current and future priorities that I have. Similar to other benchmarks and/or programs, people will move on to other projects once they have published the paper. So, for me, maintaining the repositories will depend on the interest of the community. After all, there is little benefit on improving it if no one will use it.

At the very least, in the case where I moved on to other projects, I will complete the code documentations so it can be developed by other people if they wish to fork the repository. Further, I will still answer any questions that you have at the `Issues` section.

**Q: Why is it named `ScrambleBench`?**

**A:** When developing the workflow, I realized that the workflow will depend on various `conda` environment, which as we know can be a nightmare to maintain, and this kinda `scramble` my brains out hahaha. And the fact that the word `Scramble` reminds me of the scrambled eggs led me to put the egg image in the Abstract above XD.


