#!/bin/bash


## Step 0: Parsing the configuration file first
DIR=$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" &> /dev/null && pwd)
source $DIR/parse_yaml.sh
eval $(parse_yaml $1)

eval "$(conda shell.bash hook)"

mkdir -p $generation_output/summary

# Step 1: Generating ligand through PMDM model
# The PMDM requires the pocket generation through their own script
# The cutoff is 10 A, as suggested through their Github Issue #10 (closed issue)
if [ ! -z "$model_pmdm_conda_env" ]; then
    conda activate $model_pmdm_conda_env
    time python -u $model_pmdm_dir/sample_for_pdb.py \
                --ckpt $model_pmdm_dir/500.pt  \
                --num_samples $generation_parameter_num_sample \
                --sampling_type generalized \
                --pdb_path $input_pocket_path \
                --batch_size 5 \
                --outdir $generation_output/$model_pmdm_name

    ## Process PMDM
    cat $generation_output/$model_pmdm_name/generate_ref/*.sdf > $generation_output/summary/generated_$protein_name"_ligand_"$model_pmdm_name"_"$(date -I).sdf

fi


# Step 2: Generating ligand through DiffSBDD model.
if [ ! -z "$model_diffsbdd_conda_env" ]; then
    conda activate $model_diffsbdd_conda_env
    # cd $exec_directory_diffsbdd
    mkdir -p $generation_output/$model_diffsbdd_name
    time python $model_diffsbdd_dir/generate_ligands.py \
                $model_diffsbdd_dir/checkpoints/crossdocked_fullatom_cond.ckpt \
                --n_samples $generation_parameter_num_sample \
                --pdbfile $input_complex_path \
                --outfile $generation_output/DiffSBDD/generated_${generation_parameter_name}_DiffSBDD.sdf \
                --ref_ligand $input_sdf_path \
                --batch_size 50

    ## Process DiffSBDD
    cp $generation_output/$model_diffsbdd_name/*.sdf $generation_output/summary/generated_$protein_name"_ligand_"$model_diffsbdd_name"_"$(date -I).sdf
fi

#Step 3: Generating ligand from Pocket2Mol
if [ ! -z "$model_pocket2mol_conda_env" ]; then
    conda activate $model_pocket2mol_conda_env
    time python -u $model_pocket2mol_dir/sample_for_pdb.py --pdb_path $input_complex_path \
                --center "$input_pocket_coord" \
                --bbox_size $generation_parameter_box_size \
                --outdir $generation_output/$model_pocket2mol_name \
                --config $model_pocket2mol_dir/configs/sample_for_pdb.yml \
                --checkpoint $model_pocket2mol_dir/ckpt/pretrained_Pocket2Mol.pt \
                --num_samples $generation_parameter_num_sample


    ## Process Pocket2Mol
    $DIR/combine_sdf.sh $generation_output/$model_pocket2mol_name/*/SDF #give sdf_combined.sdf as output
    cp $generation_output/$model_pocket2mol_name/*/SDF/sdf_combined.sdf $generation_output/summary/generated_$protein_name"_ligand_"$model_pocket2mol_name"_"$(date -I).sdf
fi


## Step 4: Generating ligand from PocketFlow
if [ ! -z "$model_pocketflow_conda_env" ]; then
    conda activate $model_pocketflow_conda_env
    time python $model_pocketflow_dir/main_generate.py  \
        --ckpt $model_pocketflow_dir/ckpt/ZINC-pretrained-255000.pt \
        -n $generation_parameter_num_sample \
        -d cuda:0 -at 1.0 -bt 1.0 --max_atom_num 35 -ft 0.5 -cm 0  \
        --with_print True \
        --name $generation_parameter_name \
        -pkt $input_pocket_path \
        --root_path $generation_output/$model_pocketflow_name

    ## Process PocketFlow
    cp $generation_output/$model_pocketflow_name/*/*/generated.sdf $generation_output/summary/generated_$protein_name"_ligand_"$model_pocketflow_name"_"$(date -I).sdf

fi

## Step 5: Generating ligand from Lingo3DMol
if [ ! -z "$model_lingo3dmol_conda_env" ]; then
    conda activate $model_lingo3dmol_conda_env
    echo ,,$input_pocket_path > $model_lingo3dmol_dir/datasets/lingo3dmol_dataset
    time python $model_lingo3dmol_dir/inference/inference_avoid_clash.py \
    --cuda '0' --cuda_list 0 \
    --input_list $model_lingo3dmol_dir/datasets/lingo3dmol_dataset \
    --savedir $generation_output/$model_lingo3dmol_name"_" \
    --frag_len_add 15 \
    --max_run_hours 3 \
    --gen_frag_set 40 \
    --gennums $generation_parameter_num_sample \
    --contact_path $model_lingo3dmol_dir/checkpoint/contact.pkl \
    --caption_path $model_lingo3dmol_dir/checkpoint/gen_mol.pkl

    ## Process Lingo3DMol
    $DIR/combine_mol.sh $generation_output/Lingo3DMol_0/*/
    cp $generation_output/$model_lingo3dmol_name"_0"/*/sdf_combined.sdf $generation_output/summary/generated_$protein_name"_ligand_"$model_lingo3dmol_name"_"$(date -I).sdf

fi

