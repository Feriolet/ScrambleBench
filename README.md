# README

## Introduction

![Intro graphic](asset/GenAI_Abstract.png)
Hi! Welcome to ScrambleBench, A Workflow for Comparative Assessment of Structure-based *de novo* Generative Models. This repository contains the code used for our manuscript. As this version is used for transparency in providing the workflow, it is not currently production-ready and we are currently cleaning the codes. Please look out for our version 0.1.0 soon!

## Installation

```sh
## Install Scramblebench 
conda create -n scramblebench python=3.9
pip install rdkit numpy matplotlib ptitprince seaborn pandas meeko fastparquet pyarrow
```

## Model Installation

### PMDM

We are not installing the qvina environment for our study.
The mol.yml is changed to included --extra-link-url and --find-links for pytorch
In case torch can not be imported, may need to delete the libcudart.so.11.0 in the torch/lib if libcudart.so.11.8.*.* is also present in the torch/lib.

```bash
git clone https://github.com/Layne-Huang/PMDM
cd PMDM
# change env name to benchmark_pmdm
# add - --extra-index-url https://download.pytorch.org/whl/cu118 at pip
conda env create -f mol.yml
```

### DiffSBDD
```bash
git clone https://github.com/arneschneuing/DiffSBDD
conda env create -f environment.yaml
```

### Pocket2Mol

The config of `sample_for_pdb.yaml` is changed to allow Pocket2Mol to expand its atom

### PocketFlow

```bash
conda create -n benchmark_pocketflow2 python=3.10 pymol-open-source=2.5.0 openbabel -y
pip install torch==1.13.0+cu117 torchvision==0.14.0+cu117 torchaudio==0.13.0 --extra-index-url https://download.pytorch.org/whl/cu117
pip install scipy numpy==1.23.0
pip install pyg_lib torch_scatter torch_sparse torch_cluster torch_spline_conv -f https://data.pyg.org/whl/torch-1.13.0+cu117.html --no-index
pip install rdkit tensorboard six lmdb easydict
pip install torch_geometric==2.3.1
```

### Lingo3DMol

Note: Lingo3DMol/util/pocket_code_all.py was modified to increase pocket representation length from 500 to 1000, probably because the pocket size is too big
```bash
conda create -n lingo3dmol python=3.8
conda activate lingo3dmol
conda install pytorch==1.10.1 torchvision==0.11.2 torchaudio==0.10.1 cudatoolkit=11.3 -c pytorch -c conda-forge
pip install scipy==1.7.3 pandas==1.5.1 numpy==1.20.3 rdkit==2022.09.1 psutil torch_geometric==2.3.1
#maybe add conda install pyg -c pyg if still does not work
```

### Genbench3D
```
git clone https://github.com/bbaillif/genbench3d.git
cd genbench3d
mamba env create -f environment.yml
mamba activate genbench3d
pip install ray
```

## Usage