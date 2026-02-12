


import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np

COLOR_PALLETE = ["#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#999999"]

mces_data = pd.read_csv('/opt/veincent/GenAI_manuscript/ScrambleBench/script/compare_mces_sim_gpcr.csv')
mces_colname = 'Hamiltonian Tanimoto Diversity based on MCES'
time_colname = 'Hamiltonian Tanimoto Diversity Time'

sim05_df = mces_data[mces_data['sim'] == 0.5]
sim01_df = mces_data[mces_data['sim'] == 0.1]

sim01_df = sim01_df.set_index(['Protein', 'Model', 'num_sample'])
sim05_df = sim05_df.set_index(['Protein', 'Model', 'num_sample'])

protein_list = ['GPCR_5HT2C', 'GPCR_DRD2']
fig, axes = plt.subplots(2, 2, figsize=(15, 10))
axes = axes.flatten()

i = 0
for protein in protein_list:
    filtered_data = sim01_df[mces_colname] - sim05_df[mces_colname]
    filtered_data = filtered_data.reset_index()
    filtered_data = filtered_data[filtered_data['Protein'] == protein]

    sns.barplot(data=filtered_data, x='num_sample', 
        ax=axes[i],
        y=mces_colname,
        hue='Model', 
        palette= COLOR_PALLETE)
    
    axes[i].set_title(protein)
    axes[i].set_ylabel('ΔHamDiv MCES')
    axes[i].axhline(0, color='black', linestyle='-', linewidth=1.5)

    handles, labels = axes[i].get_legend_handles_labels()
    axes[i].get_legend().remove()
    axes[i].legend(handles, labels, ncol=2, framealpha=0.5)

    i += 1

for protein in protein_list:
    filtered_data = sim01_df[time_colname] - sim05_df[time_colname]
    filtered_data = filtered_data.reset_index()
    filtered_data = filtered_data[filtered_data['Protein'] == protein]

    sns.barplot(data=filtered_data, x='num_sample', 
        ax=axes[i],
        y=time_colname,
        hue='Model', 
        palette= COLOR_PALLETE)

    axes[i].set_title(protein)
    axes[i].set_ylabel('ΔTime for HamDiv MCES / s')
    axes[i].axhline(0, color='black', linestyle='-', linewidth=1.5)

    handles, labels = axes[i].get_legend_handles_labels()
    axes[i].get_legend().remove()
    axes[i].legend(handles, labels, framealpha=0.5)

    i += 1
plt.savefig('hi.png')