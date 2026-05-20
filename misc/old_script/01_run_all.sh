#!/bin/bash


## Step 0: Parsing the configuration file first
source parse_yaml.sh
eval $(parse_yaml $1)

eval "$(conda shell.bash hook)"


# Step 1: Generating ligand through PMDM model
# The PMDM requires the pocket generation through their own script
# The cutoff is 10 A, as suggested through their Github Issue #10 (closed issue)
if [ ! -z "$conda_env_pmdm" ]; then
    conda activate $conda_env_pmdm
    time python -u $models_dir_pmdm/sample_for_pdb.py \
                --ckpt $models_dir_pmdm/500.pt  \
                --num_samples $parameter_num_sample \
                --sampling_type generalized \
                --pdb_path $input_pocket_path \
                --batch_size 5 \
                --outdir $output_output_dir/PMDM
fi


# Step 2: Generating ligand through DiffSBDD model.
if [ ! -z "$conda_env_diffsbdd" ]; then
    conda activate $conda_env_diffsbdd
    # cd $exec_directory_diffsbdd
    mkdir -p $output_output_dir/DiffSBDD
    time python $models_dir_diffsbdd/generate_ligands.py \
                $models_dir_diffsbdd/checkpoints/crossdocked_fullatom_cond.ckpt \
                --n_samples $parameter_num_sample \
                --pdbfile $input_complex_path \
                --outfile $output_output_dir/DiffSBDD/generated_${parameter_name}_DiffSBDD.sdf \
                --ref_ligand $input_sdf_path \
                --batch_size 50
fi

#Step 3: Generating ligand from Pocket2Mol
if [ ! -z "$conda_env_pocket2mol" ]; then
    conda activate $conda_env_pocket2mol
    time python -u $models_dir_pocket2mol/sample_for_pdb.py --pdb_path $input_complex_path \
                --center "$parameter_pocket_coord" \
                --bbox_size $parameter_box_size \
                --outdir $output_output_dir/Pocket2Mol \
                --config $models_dir_pocket2mol/configs/sample_for_pdb.yml \
                --checkpoint $models_dir_pocket2mol/ckpt/pretrained_Pocket2Mol.pt \
                --num_samples $parameter_num_sample
fi


## Step 4: Generating ligand from PocketFlow
if [ ! -z "$conda_env_pocketflow" ]; then
    conda activate $conda_env_pocketflow
    time python $models_dir_pocketflow/main_generate.py  \
        --ckpt $models_dir_pocketflow/ckpt/ZINC-pretrained-255000.pt \
        -n $parameter_num_sample \
        -d cuda:0 -at 1.0 -bt 1.0 --max_atom_num 35 -ft 0.5 -cm 0  \
        --with_print True \
        --name $parameter_name \
        -pkt $input_pocket_path \
        --root_path $output_output_dir/PocketFlow
fi

## Step 5: Generating ligand from Lingo3DMol
if [ ! -z "$conda_env_lingo3dmol" ]; then
    conda activate $conda_env_lingo3dmol
    echo ,,$input_pocket_path > $models_dir_lingo3dmol/datasets/lingo3dmol_dataset
    time python $models_dir_lingo3dmol/inference/inference_avoid_clash.py \
    --cuda '0' --cuda_list 0 \
    --input_list $models_dir_lingo3dmol/datasets/lingo3dmol_dataset \
    --savedir $output_output_dir/Lingo3DMol_ \
    --frag_len_add 15 \
    --max_run_hours 3 \
    --gen_frag_set 40 \
    --gennums $parameter_num_sample \
    --contact_path $models_dir_lingo3dmol/checkpoint/contact.pkl \
    --caption_path $models_dir_lingo3dmol/checkpoint/gen_mol.pkl
fi