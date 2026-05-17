SCRAMBLEBENCH_SCRIPT_DIR='../../src/scramblebench/script'

bash ./prepare_config.sh

python $SCRAMBLEBENCH_SCRIPT_DIR/p2_execute_generation.py -i output/yaml_list.txt
python $SCRAMBLEBENCH_SCRIPT_DIR/p3_prepare_molecule.py -i output/yaml_list.txt
python $SCRAMBLEBENCH_SCRIPT_DIR/p4_analyse_diversity.py -i output/yaml_list.txt
python $SCRAMBLEBENCH_SCRIPT_DIR/p4_analyse_genbench3d.py -i output/yaml_list.txt
python $SCRAMBLEBENCH_SCRIPT_DIR/p4_analyse_redocking.py -i output/yaml_list.txt
python $SCRAMBLEBENCH_SCRIPT_DIR/p5_collect_data_analysis.py -i output/yaml_list.txt