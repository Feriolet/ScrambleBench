#!/bin/bash

## Step 0: Parsing the configuration for this script
## Before running this, please make sure that the directory is new
## That is, existing output should be removed if you want to rerun the script
## with the same protein target or input files.
## Otherwise, you may end up doubling or even tripling your compound
## number with duplicates


source parse_yaml.sh
eval $(parse_yaml $1)

eval "$(conda shell.bash hook)"

mkdir -p $output_output_dir/summary

## Process PMDM
if [ ! -z "$conda_env_pmdm" ]; then
    cat $output_output_dir/PMDM/generate_ref/*.sdf > $output_output_dir/summary/generated_$parameter_name"_ligand_PMDM_"$(date -I).sdf
fi

## Process DiffSBDD
if [ ! -z "$conda_env_diffsbdd" ]; then
    cp $output_output_dir/DiffSBDD/*.sdf $output_output_dir/summary/generated_$parameter_name"_ligand_DiffSBDD_"$(date -I).sdf
fi

## Process Pocket2Mol
if [ ! -z "$conda_env_pocket2mol" ]; then
    ./combine_sdf.sh $output_output_dir/Pocket2Mol/*/SDF #give sdf_combined.sdf as output
    cp $output_output_dir/Pocket2Mol/*/SDF/sdf_combined.sdf $output_output_dir/summary/generated_$parameter_name"_ligand_Pocket2Mol_"$(date -I).sdf
fi

## Process PocketFlow
if [ ! -z "$conda_env_pocketflow" ]; then
    cp $output_output_dir/PocketFlow/*/*/generated.sdf $output_output_dir/summary/generated_$parameter_name"_ligand_PocketFlow_"$(date -I).sdf
fi

## Process Lingo3DMol
if [ ! -z "$conda_env_lingo3dmol" ]; then
    ./combine_mol.sh $output_output_dir/Lingo3DMol_0/*/
    cp $output_output_dir/Lingo3DMol_0/*/sdf_combined.sdf $output_output_dir/summary/generated_$parameter_name"_ligand_Lingo3DMol_"$(date -I).sdf
fi