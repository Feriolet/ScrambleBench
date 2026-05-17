SCRAMBLEBENCH_SCRIPT_DIR='../../src/scramblebench/script'

python $SCRAMBLEBENCH_SCRIPT_DIR/p2_execute_generation.py -i Protein_5HT2A/config_test_run_clean_config.yml
python $SCRAMBLEBENCH_SCRIPT_DIR/p3_prepare_molecule.py -i Protein_5HT2A/config_test_run_clean_config.yml
python $SCRAMBLEBENCH_SCRIPT_DIR/p4_analyse_diversity.py -i Protein_5HT2A/config_test_run_clean_config.yml
python $SCRAMBLEBENCH_SCRIPT_DIR/p4_analyse_genbench3d.py -i Protein_5HT2A/config_test_run_clean_config.yml