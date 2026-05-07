from pathlib import Path

INPUT_KEY = 'input'
INPUT_DIR_KEY = 'input_dir'
MODEL_KEY = 'model'
GENERATION_KEY = 'generation'
GENERATION_PARAMETER_KEY = 'parameter'
POST_GENERATION_KEY = 'post_generation'
ANALYSIS_KEY = 'analysis'
ANALYSIS_GENBENCH3D_KEY = 'genbench3d'

## Generation
GENERATION_FOLDER = 'AI_Generation'
GENERATION_TEMPLATE_SCRIPT_PATH = Path(Path(__file__).parent.parent / 'utils' / 'generation_template.sh').resolve()
MODEL_CONDA_ENV_KEY = 'conda_env'
GENERATION_TEMPLATE_SCRIPT_FILE_KEY = 'script_pathfile'