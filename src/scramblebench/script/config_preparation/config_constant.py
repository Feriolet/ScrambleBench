from pathlib import Path

INPUT_KEY = 'input'
INPUT_DIR_KEY = 'input_dir'
MODEL_KEY = 'model'
GENERATION_KEY = 'generation'
GENERATION_PARAMETER_KEY = 'parameter'
PARAMETER_KEY = 'parameter'
POST_GENERATION_KEY = 'post_generation'
ANALYSIS_KEY = 'analysis'
ANALYSIS_GENBENCH3D_KEY = 'genbench3d'
ANALYSIS_REDOCKING_KEY = 'redocking'
ANALYSIS_REDOCKING_PROTONATION_KEY = 'protonation'
ANALYSIS_REDOCKING_DOCKING_KEY = 'docking'

ANALYSIS_DOCKING_EASYDOCK_KEY = 'easydock'
ANALYSIS_DOCKING_SCHRODINGER_KEY = 'glide'

ANALYSIS_DIVERSITY_KEY = 'diversity'
ANALYSIS_VIRTUAL_HIT_KEY = 'virtual_hit'

## Generation
GENERATION_FOLDER = 'AI_Generation'
GENERATION_TEMPLATE_SCRIPT_PATH = Path(Path(__file__).parent.parent / 'utils/generation_utils/generation_template.sh').resolve()
MODEL_CONDA_ENV_KEY = 'conda_env'
GENERATION_TEMPLATE_SCRIPT_FILE_KEY = 'script_pathfile'

CONFIG_BATCH_PARAMETERS = [{'key':[INPUT_KEY],
                            'type': 'dict'},
                           {'key':[GENERATION_KEY, GENERATION_PARAMETER_KEY, 'num_sample'],
                            'type': 'int'},
                           {'key':[GENERATION_KEY, GENERATION_PARAMETER_KEY, 'box_size'],
                            'type': 'float'}]

REPORT_KEY = 'report'