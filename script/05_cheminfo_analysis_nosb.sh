#!/bin/bash
#2023.09.2
conda activate genbench3d

time python /opt/veincent/GenAI_manuscript/models/genbench3d/benchmark_mols.py \
-c /opt/veincent/GenAI_manuscript/models/genbench3d/config/test1.yaml \
-i /opt/Enamine_4M_Library/prepared_Enamine_4M_Library_single_state.sdf \
-o please_remove.json

#-i /opt/Enamine_4M_Library/prepared_Enamine_4M_Library_single_state.sdf \
# -i /opt/veincent/FYP_GenAIChem/AKT1_docking_finetuning_control/vsw_akt1_4M_Enamine-SP_OUT_1_pv.sdf \