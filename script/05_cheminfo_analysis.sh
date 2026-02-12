#!/bin/bash

source parse_yaml.sh
eval $(parse_yaml $1)
eval $(parse_yaml $genbench_config)
eval "$(conda shell.bash hook)"
conda activate benchmark_genbench3d

## To do list:
## 1. Make genai_model_array not hardcoded to current models
## 2. Maybe prepare another config file so user does not need to manually change this file? Need to gauge interest for this

## File setting, please change the genbench3d according to your own directories
genbench_config=$(realpath ../models/genbench3d/config/GenAI_evaluation.yaml)
genbench_dir=$(realpath ../models/genbench3d)
cheminformatic_input_dir=$output_output_dir/cheminformatics_input_prepared
cheminformatic_output_dir=$output_output_dir/cheminformatics_analysis

## Edge case where schrodinger need to be started first
schrodinger_path_exec='/opt/schrodinger2025-4'
export SCHRODINGER=$schrodinger_path_exec
jsc_status=$($schrodinger_path_exec/jsc 'local-server-status')
if [[ "$jsc_status" == *"STOPPED"* ]]; then
    $schrodinger_path_exec/jsc local-server-start
fi

## Array of models used for generating ligands
declare -a genai_model_array=("PMDM" "DiffSBDD" "PocketFlow" "Pocket2Mol" "Lingo3DMol" "Chem42")

#Only fill 1 or 0 please
is_complex_forcefield_minimisation=1
is_complex_forcefield_unminimisation=1
is_cancel_protonation_by_obabel_or_adfr=1
###############################################################################################################################################################
declare -a minimisation=()
if [ $is_complex_forcefield_minimisation -eq 1 ]; then minimisation+=("complex_forcefield_minimised"); fi
if [ $is_complex_forcefield_unminimisation -eq 1 ]; then minimisation+=("complex_forcefield_unminimised"); fi

# sanity check at maximum one file is present in the directory, otherwise terminate the script
for genai_model_name in "${genai_model_array[@]}"
do
    sdf_regex=$cheminformatic_input_dir/*$genai_model_name*.sdf
    num_sdf_fname=$(find $cheminformatic_input_dir -type f -name *$genai_model_name*.sdf | wc -l)
    if (($num_sdf_fname > 1)); then
        echo $genai_model_name has more than 1 file. Exiting...
        exit 1
    fi
done

## Copying input file of both protein pdb and ligand sdf to the new directory, so the intermediate files won't be saved in the input directory
## Input directory should only have raw input, no intermediate file should exist there in my opinion
mkdir -p $cheminformatic_output_dir/{json_output,structure_input}
cheminformatics_input_pdb=$cheminformatic_output_dir/structure_input/$(basename $input_pdb_path)
cheminformatics_input_sdf=$cheminformatic_output_dir/structure_input/$(basename $input_sdf_path)
cp $input_pdb_path $cheminformatics_input_pdb
cp $input_sdf_path $cheminformatics_input_sdf

# now loop through the above array
for genai_model_name in "${genai_model_array[@]}"
do

    for is_minimised in "${minimisation[@]}"
    do   
        mkdir -p $cheminformatic_output_dir/{Glide,Vina}/$genai_model_name/$is_minimised
        echo generating ligand for $parameter_name $parameter_num_sample with model $genai_model_name using $is_minimised
        generated_ligand_sdf_fn=$cheminformatic_input_dir/generated_$parameter_name\_ligand*$genai_model_name*.sdf
        cheminformatic_output_json_fn=$cheminformatic_output_dir/json_output/generated_$parameter_name\_ligand_$genai_model_name\_$is_minimised\_$(date -I).json
        cheminformatic_output_stdout_fn=$cheminformatic_output_dir/json_output/generated_$parameter_name\_ligand_$genai_model_name\_$is_minimised\_$(date -I).log
        cheminformatics_output_log_fn=$cheminformatic_output_dir/json_output/generated_$parameter_name\_ligand_$genai_model_name\_$is_minimised\_$(date -I)_sb_benchmark.log

        time python $genbench_dir/sb_benchmark_mols.py \
        -c $genbench_config \
        -i $generated_ligand_sdf_fn \
        -o $cheminformatic_output_json_fn \
        -p $cheminformatics_input_pdb \
        -n $cheminformatics_input_sdf \
        --do_conf_analysis \
        --glide \
        --output_vina_dir $cheminformatic_output_dir/Vina/$genai_model_name/$is_minimised \
        --output_glide_dir $cheminformatic_output_dir/Glide/$genai_model_name/$is_minimised \
        --log_output $cheminformatics_output_log_fn \
        $(if [[ $is_minimised == 'complex_forcefield_minimised' ]]; then echo '-m'; fi) \
        $(if [ $is_complex_forcefield_unminimisation -eq 1 ]; then echo '--cancel_protonation'; fi) \
        >> $cheminformatic_output_stdout_fn
    done
done
