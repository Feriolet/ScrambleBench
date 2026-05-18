### About

This example is taken from the POKMOL dataset, which only provides the raw generated ligand of the molecules. To run these kinds of examples, the folder should be restructured to the one in the example. Specifically, all of the generated ligand must be combined into one single `.sdf` file containing the name of the model (e.g., 'GraphBP' model should have a file called `*GraphBP*.sdf`). These files should be within a single folder named `summary`, because this is the hard-coded folder in `ScrambleBench`.

Currently, there is no script to automate the config generation without protein input, so these config file and `yaml_list.txt` was generated manually. Nevertheless, this should be easily supported in future release.

Because no protein input was gathered, analysis using these input could not be supported. The config files should reflect these changes, as seen in `Protein_5HT2A/config_test_run_clean_config.yml`
### Analysis

```bash
chmod +x ./run_scramblebench.sh
./run_scramblebench.sh
```