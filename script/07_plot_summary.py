import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
import sys
import ptitprince as pt
import warnings

from enum import IntEnum, Enum
from pathlib import Path
from dataclasses import dataclass
from rdkit import Chem
from matplotlib.ticker import MaxNLocator
import rdkit # type: ignore
from rdkit.Chem import rdMolTransforms
from typing import Union, Iterable
COLORBLIND_PALETTE = ["#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#999999"]
GENAI_MODEL_ORDER = ['Pocket2Mol', 'PocketFlow', 'Lingo3DMol', 'DiffSBDD', 'PMDM', 'Chem42']

INPUT_MOL_BLOCK_VINA_COLUMN = 'input_mol_block_vina'
INPUT_MOL_BLOCK_GLIDE_COLUMN = 'input_mol_block_glide'
class DockingMethodAnalysis(Enum):
    VINA = 'vina'
    GLIDE = 'glide'
    NONE = 'none'

class NumberComparison(Enum):
    LESSTHAN = 'vina'
    LESSTHANOREQUAL = 'glide'
    NONE = 'none'

class SummaryDataFrameWrapper:

    def __init__(self, dataframe: Union[pd.DataFrame, pd.Series]) -> None:
        self._dataframe = dataframe
    
    @property
    def protein_name(self): 
        return self._dataframe['protein_name']

    @protein_name.setter
    def protein_name(self, value):
        self._dataframe['protein_name'] = value

    @property
    def genai_model(self):
        return self._dataframe['genai_model']

    @genai_model.setter
    def genai_model(self, value):
        self._dataframe['genai_model'] = value

    @property
    def num_sample(self):
        return self._dataframe['num_sample']

    @num_sample.setter
    def num_sample(self, value):
        self._dataframe['num_sample'] = value

    @property
    def ligand_id(self):
        return self._dataframe['ligand_id']

    @ligand_id.setter
    def ligand_id(self, value):
        self._dataframe['ligand_id'] = value
    
    @property
    def stereo_id(self):
        return self._dataframe['stereo_id']
    
    @stereo_id.setter
    def stereo_id(self, value):
        self._dataframe['stereo_id'] = value
    
    @property
    def molecular_weight(self):
        return self._dataframe.get('MW / Da')
    
    @molecular_weight.setter
    def molecular_weight(self, value):
        self._dataframe['MW / Da'] = value
    
    @property
    def log_p(self):
        return self._dataframe.get('log P')
    
    @log_p.setter
    def log_p(self, value):
        self._dataframe['log P'] = value
    
    @property
    def sa_score(self):
        return self._dataframe.get('SA_Score')
    
    @sa_score.setter
    def sa_score(self, value):
        self._dataframe['SA_Score'] = value
    
    @property
    def qed(self):
        return self._dataframe.get('QED')
    
    @qed.setter
    def qed(self, value):
        self._dataframe['QED'] = value

    @property
    def redocked_glide_docking_score(self):
        return self._dataframe.get('redocked_glide_docking_score')
    
    @redocked_glide_docking_score.setter
    def redocked_glide_docking_score(self, value):
        self._dataframe['redocked_glide_docking_score'] = value

    @property
    def glide_phase_score(self):
        return self._dataframe.get('glide_phase_score')
    
    @glide_phase_score.setter
    def glide_phase_score(self, value):
        self._dataframe['glide_phase_score'] = value

    @property
    def redocked_glide_rms(self):
        return self._dataframe.get('glide_redock-input_rms')
    
    @redocked_glide_rms.setter
    def redocked_glide_rms(self, value):
        self._dataframe['glide_redock-input_rms'] = value

    @property
    def redocked_glide_mol_block(self):
        return self._dataframe.get('redocked_glide_mol_block')
    
    @redocked_glide_mol_block.setter
    def redocked_glide_mol_block(self, value):
        self._dataframe['redocked_glide_mol_block'] = value

    @property
    def FF_unminimised_glide_inplace_docking_score(self):
        return self._dataframe.get('FF_unminimised_glide_inplace_docking_score')
    
    @FF_unminimised_glide_inplace_docking_score.setter
    def FF_unminimised_glide_inplace_docking_score(self, value):
        self._dataframe['FF_unminimised_glide_inplace_docking_score'] = value

    @property
    def FF_unminimised_glide_mininplace_docking_score(self):
        return self._dataframe.get('FF_unminimised_glide_mininplace_docking_score')
    
    @FF_unminimised_glide_mininplace_docking_score.setter
    def FF_unminimised_glide_mininplace_docking_score(self, value):
        self._dataframe['FF_unminimised_glide_mininplace_docking_score'] = value

    @property
    def FF_minimised_glide_mininplace_docking_score(self):
        return self._dataframe.get('FF_minimised_glide_mininplace_docking_score')
    
    @FF_minimised_glide_mininplace_docking_score.setter
    def FF_minimised_glide_mininplace_docking_score(self, value):
        self._dataframe['FF_minimised_glide_mininplace_docking_score'] = value

    @property
    def FF_minimised_glide_mininplace_docked_mol_block(self):
        return self._dataframe.get('FF_minimised_glide_mininplace_docked_mol_block')
    
    @FF_minimised_glide_mininplace_docked_mol_block.setter
    def FF_minimised_glide_mininplace_docked_mol_block(self, value):
        self._dataframe['FF_minimised_glide_mininplace_docked_mol_block'] = value

    @property
    def redocked_vina_docking_score(self):
        return self._dataframe.get('redocked_vina_docking_score')
    
    @redocked_vina_docking_score.setter
    def redocked_vina_docking_score(self, value):
        self._dataframe['redocked_vina_docking_score'] = value

    @property
    def redocked_vina_rms(self):
        return self._dataframe.get('vina_redock-input_rms')
    
    @redocked_vina_rms.setter
    def redocked_vina_rms(self, value):
        self._dataframe['vina_redock-input_rms'] = value

    @property
    def redocked_vina_mol_block(self):
        return self._dataframe.get('redocked_vina_mol_block')
    
    @redocked_vina_mol_block.setter
    def redocked_vina_mol_block(self, value):
        self._dataframe['redocked_vina_mol_block'] = value

    @property
    def FF_unminimised_vina_inplace_docking_score(self):
        return self._dataframe.get('FF_unminimised_vina_inplace_docking_score')
    
    @FF_unminimised_vina_inplace_docking_score.setter
    def FF_unminimised_vina_inplace_docking_score(self, value):
        self._dataframe['FF_unminimised_vina_inplace_docking_score'] = value

    @property
    def FF_minimised_vina_mininplace_docked_mol_block(self):
        return self._dataframe.get('FF_minimised_vina_mininplace_docked_mol_block')
    
    @FF_minimised_vina_mininplace_docked_mol_block.setter
    def FF_minimised_vina_mininplace_docked_mol_block(self, value):
        self._dataframe['FF_minimised_vina_mininplace_docked_mol_block'] = value

    @property
    def input_mol_block(self):
        return self._dataframe.get('input_mol_block')
    
    @input_mol_block.setter
    def input_mol_block(self, value):
        self._dataframe['input_mol_block'] = value

    @property
    def input_mol_block_vina(self):
        return self._dataframe.get('input_mol_block_vina')
    
    @input_mol_block_vina.setter
    def input_mol_block_vina(self, value):
        self._dataframe['input_mol_block_vina'] = value

    @property
    def input_mol_block_glide(self):
        return self._dataframe.get('input_mol_block_glide')
    
    @input_mol_block_glide.setter
    def input_mol_block_glide(self, value):
        self._dataframe['input_mol_block_glide'] = value

    @property
    def centroid_distance_docking_glide(self):
        return self._dataframe.get('centroid_distance_docking_glide')
    
    @centroid_distance_docking_glide.setter
    def centroid_distance_docking_glide(self, value):
        self._dataframe['centroid_distance_docking_glide'] = value

    def __getattr__(self, attr):
        return getattr(self._dataframe, attr)


    @property
    def get_generic_feature_colname(self):
        return [self.protein_name.name,
                self.genai_model.name,
                self.num_sample.name,
                self.ligand_id.name,
                self.stereo_id.name]

    @property
    def get_generic_feature(self):
        return SummaryDataFrameWrapper(self._dataframe[self.get_generic_feature_colname])
    
    @property
    def get_physicochemical_feature(self):
        
        physicochemical_property = [self.molecular_weight,
                                    self.log_p,
                                    self.sa_score,
                                    self.qed]
        matched_colname_list = [col.name for col in physicochemical_property if isinstance(col, pd.Series)]

        for col in matched_colname_list:
            self._dataframe[col] = self._dataframe[col].astype(float) 

        return matched_colname_list, SummaryDataFrameWrapper(self._dataframe[self.get_generic_feature_colname + matched_colname_list])

    @property
    def get_docking_feature(self):
        
        docking_property = [self.redocked_glide_docking_score,
                            self.redocked_vina_docking_score,
                            self.redocked_glide_rms,
                            self.redocked_vina_rms,
                            self.FF_unminimised_glide_inplace_docking_score,
                            self.FF_unminimised_glide_mininplace_docking_score,
                            self.FF_unminimised_vina_inplace_docking_score,
                            self.FF_minimised_glide_mininplace_docking_score,
                            self.redocked_glide_mol_block,
                            self.redocked_vina_mol_block,
                            self.input_mol_block,
                            self.molecular_weight]
        
        matched_colname_list = [col.name for col in docking_property if isinstance(col, pd.Series)]

        for col in matched_colname_list:
            try:
                self._dataframe[col] = self._dataframe[col].astype(float) 
            except ValueError:
                pass
        return matched_colname_list, SummaryDataFrameWrapper(self._dataframe[self.get_generic_feature_colname + matched_colname_list])

    @property
    def get_best_molecule_feature(self):

        best_molecule_glide_property = [self.glide_phase_score,
                                        self.redocked_glide_docking_score,
                                        self.redocked_glide_rms]
        best_molecule_vina_property = [self.redocked_vina_docking_score,
                                       self.redocked_vina_rms]
        best_molecule_physicochemical_property = [self.qed]

        matched_colname_list = []
        

        matched_glide_property = [col.name for col in best_molecule_glide_property if isinstance(col, pd.Series)]
        matched_vina_property = [col.name for col in best_molecule_vina_property if isinstance(col, pd.Series)]
        matched_physicochemical_property = [col.name for col in best_molecule_physicochemical_property if isinstance(col, pd.Series)]
        
        matched_colname_list = matched_colname_list + matched_glide_property + matched_vina_property + matched_physicochemical_property

        for col in matched_colname_list:
            try:
                self._dataframe[col] = self._dataframe[col].astype(float) 
            except ValueError:
                pass
        return matched_colname_list, SummaryDataFrameWrapper(self._dataframe[self.get_generic_feature_colname + matched_colname_list])

    def get_feature(self, columns: Iterable[Union[None, pd.Series, pd.Categorical]]):

        matched_colname_list = [col.name for col in columns if isinstance(col, pd.Series)]

        for col in matched_colname_list:
            try:
                self._dataframe[col] = self._dataframe[col].astype(float) 
            except ValueError:
                continue
        return matched_colname_list, SummaryDataFrameWrapper(self._dataframe[self.get_generic_feature_colname + matched_colname_list])

    
    def get_feature_as_dataframe(self, columns: list[Union[None, pd.Series, pd.Categorical]]):
        matched_colname_list = [col.name for col in columns if isinstance(col, pd.Series)]
        return self._dataframe[matched_colname_list]

    def filter_exact_feature(self, series, filters: list):
        return SummaryDataFrameWrapper(self._dataframe[series.isin(filters)].copy())

    def filter_stereo_id_based_on_value(self, column_list: Iterable[Union[None, pd.Series]], sorted_column: pd.Series, ascending:bool=True):
        
        _, wrapped_selected_df = self.get_feature(column_list)
        selected_df = wrapped_selected_df.to_dataframe()
        
        generic_colname_without_stereo_id = list(set(wrapped_selected_df.get_generic_feature_colname) - set([wrapped_selected_df.stereo_id.name]))
        if isinstance(sorted_column, pd.Series):
            if sorted_column.name in selected_df.columns and isinstance(sorted_column.name, str):
                return selected_df.sort_values(by=sorted_column.name, ascending=ascending).drop_duplicates(subset=generic_colname_without_stereo_id, keep='first')

    def fetch_mol(self, matched_row):
        summary_column_list = self.to_dataframe().columns
        summary_column_list_without_mol_block = [column_name for column_name in summary_column_list if 'mol_block' not in column_name]

        possible_mol_block_list = [self.redocked_glide_mol_block,
                            self.redocked_vina_mol_block,
                            self.FF_minimised_glide_mininplace_docked_mol_block,
                            self.FF_minimised_vina_mininplace_docked_mol_block,
                            self.input_mol_block]
        
        matched_mol_block_property = [col.name for col in possible_mol_block_list if isinstance(col, pd.Series)][0]
        
        matched_wrapped_df = self.filter_exact_feature(self.ligand_id, [matched_row[self.ligand_id.name]])
        matched_wrapped_df = matched_wrapped_df.filter_exact_feature(matched_wrapped_df.protein_name, [matched_row[self.protein_name.name]])
        matched_wrapped_df = matched_wrapped_df.filter_exact_feature(matched_wrapped_df.stereo_id, [matched_row[self.stereo_id.name]])

        matched_df = matched_wrapped_df.to_dataframe()

        mol = Chem.MolFromMolBlock(matched_df[matched_mol_block_property].item(), removeHs=False)
        mol.SetProp('_Name', str(matched_df[matched_wrapped_df.ligand_id.name].item()))
        for column in summary_column_list_without_mol_block:
            mol.SetProp(column, str(matched_df[column].item()))
        return mol

    def copy(self):
        assert isinstance(self._dataframe, pd.DataFrame)
        return SummaryDataFrameWrapper(self._dataframe.copy())
    
    def to_dataframe(self):
        assert isinstance(self._dataframe, pd.DataFrame)
        return self._dataframe


    def calculate_redocking_centroid_distance(self, 
                                                    reference_mol_block_col: pd.Series,
                                                    redocking_mol_block_col: pd.Series) -> float:
        if not isinstance(reference_mol_block_col, pd.Series):
            raise KeyError(f'The dataframe {reference_mol_block_col} column is not a pd.Series')
        if not isinstance(redocking_mol_block_col, pd.Series):
            raise KeyError(f'The dataframe {redocking_mol_block_col} column is not a pd.Series')
        
        com_list = []
        for index, row in self.to_dataframe().iterrows():
            try:
                reference_mol = Chem.MolFromMolBlock(row[reference_mol_block_col.name])
                redocking_mol = Chem.MolFromMolBlock(row[redocking_mol_block_col.name])

                reference_centroid =  np.array(list(rdMolTransforms.ComputeCentroid(reference_mol.GetConformer(), ignoreHs=True)))
                redocking_centroid =  np.array(list(rdMolTransforms.ComputeCentroid(redocking_mol.GetConformer(), ignoreHs=True)))
                com_list.append(np.linalg.norm(reference_centroid - redocking_centroid))
            except TypeError:
                print(f'skipping {row["ligand_id"]}')
                com_list.append(None)
        self.centroid_distance_docking_glide = com_list

def plot_redocking_score(wrapped_dataframe: SummaryDataFrameWrapper, property_list, output_prefix, reference_score=None):
    fig, axs = plt.subplots(1, len(property_list), figsize=(38,6))
    sns.set(style="ticks", rc={"lines.linewidth": 0.7})

    plt.rcParams['figure.dpi'] = 300
    if reference_score:
        reference_array  = [list(row) for row in zip(*reference_score)]
    for i, (docking_property, ax) in enumerate(zip(property_list, np.ravel(axs)[:4])):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ax = pt.RainCloud(hue=wrapped_dataframe.genai_model.name, # type: ignore
                            y=docking_property, 
                            x=wrapped_dataframe.protein_name.name, # type: ignore
                            palette=COLORBLIND_PALETTE,  
                            data=wrapped_dataframe.to_dataframe(), 
                            bw = 0.2,  
                            pointplot = True,
                            width_viol = 0.5, 
                            ax = ax, 
                            orient = 'v', 
                            move = 0.2, 
                            alpha = .7,  
                            dodge = True, 
                            box_zorder = 2, 
                            linewidth=2)   
            if reference_score:
                ax.plot(reference_array[i], marker='*', linestyle='None', markersize=20, color='black') # linestyle='None' prevents drawing lines
            handles, labels = ax.get_legend_handles_labels()
            ax.get_legend().remove()
            ax.set(xlabel=None)
            if i < 2:
                ax.yaxis.set_major_locator(MaxNLocator(integer=True))
                ax.set_ylim(-17, 4)
            else:
                ax.set_ylim(0, 31)
            ax.set_xlabel(str(wrapped_dataframe.protein_name.name).replace('_', ' '), fontsize=19)
            ax.set_ylabel(docking_property.replace('_', ' '), fontsize=19)
            ax.tick_params(axis='both', which='major', labelsize=19)
    fig.legend(handles[:6], labels[:6], loc="upper center", ncol=6, prop={'size': 19})
    #plt.title(title)
    num_sample = list(set(wrapped_dataframe.num_sample))[0]
    plt.savefig(f'{output_prefix}_numsample_{num_sample}_ReDockingScore.png', bbox_inches='tight')
    

def plot_physicochemical_property(wrapped_dataframe: SummaryDataFrameWrapper, property_list, output_prefix):
    fig, axs = plt.subplots(2,2, figsize=(17,12))
    sns.set(style="ticks", rc={"lines.linewidth": 0.7})
    plt.rcParams['figure.dpi'] = 300
    for physicochemical_property, ax in zip(property_list, np.ravel(axs)[:4]):
        with warnings.catch_warnings():
            warnings.simplefilter("ignore")
            ax = pt.RainCloud(hue=wrapped_dataframe.genai_model.name, # type: ignore
                            y=physicochemical_property, 
                            x=wrapped_dataframe.protein_name.name, # type: ignore
                            palette=COLORBLIND_PALETTE,  
                            data=wrapped_dataframe.to_dataframe(), 
                            bw = 0.2,  
                            pointplot = True,
                            width_viol = 0.5, 
                            ax = ax, 
                            orient = 'v', 
                            move = 0.2, 
                            alpha = .7,  
                            dodge = True, 
                            box_zorder = 2, 
                            linewidth=2)   
            
            handles, labels = ax.get_legend_handles_labels()
            ax.get_legend().remove()
            # ax.set(xlabel=None)
            # ax.set_ylim([-15, 10])
            ax.set_xlabel(str(wrapped_dataframe.protein_name.name).replace('_', ' '), fontsize=18)
            ax.set_ylabel(physicochemical_property.replace('_', ' '), fontsize=18)
            ax.tick_params(axis='both', which='major', labelsize=18)
    fig.legend(handles[:6], labels[:6], loc="upper center", ncol=6, prop={'size': 18}, bbox_to_anchor=(0.5, 0.95))

    #plt.title(title)
    num_sample = list(set(wrapped_dataframe.num_sample))[0]
    plt.savefig(f'{output_prefix}_numsample_{num_sample}_PhysicoChemicalProperty.svg', bbox_inches='tight')
    

def plot_glide_undocked_freq(wrapped_dataframe: SummaryDataFrameWrapper, column_name: pd.Series, output_prefix: str):
    fig, axs = plt.subplots(1, figsize=(17,12))

    GLIDE_SKIPPED_LIGAND_DOCKING_SCORE = 10000
    #protein_wrapped_df = wrapped_dataframe.filter_exact_feature(wrapped_dataframe.protein_name, [protein])
    protein_df = wrapped_dataframe.to_dataframe()

    grouped_df = protein_df.groupby([wrapped_dataframe.protein_name.name, 
                                     wrapped_dataframe.genai_model.name]).size().reset_index(name='counts')

    grouped_undocked_glide_df = protein_df[column_name == GLIDE_SKIPPED_LIGAND_DOCKING_SCORE].groupby([wrapped_dataframe.protein_name.name, 
                                                                                                    wrapped_dataframe.genai_model.name,
                                                                                                    column_name.name]).size().reset_index(name='glide undocked counts')

    grouped_undocked_glide_df = grouped_undocked_glide_df.merge(grouped_df)
    grouped_undocked_glide_df['undocked_percentage'] = grouped_undocked_glide_df['glide undocked counts'] / grouped_undocked_glide_df['counts']

    sns.barplot(data=grouped_undocked_glide_df,
                x=wrapped_dataframe.genai_model.name,
                y='undocked_percentage',
                hue=wrapped_dataframe.protein_name.name,
                ax=axs,
                palette=COLORBLIND_PALETTE)
    num_sample = list(set(wrapped_dataframe.num_sample))[0]
    axs.set_ylabel(f'Undocked ligand proportion {output_prefix.replace("_", " ")}', fontsize=25)
    axs.set_xlabel(wrapped_dataframe.genai_model.name, fontsize=25)
    axs.tick_params(axis='both', which='major', labelsize=23)

    handles, labels = axs.get_legend_handles_labels()
    axs.get_legend().remove()
    
    fig.legend(handles[:6], labels[:6], loc="upper center", ncol=6, prop={'size': 25}, bbox_to_anchor=(0.5, 1.15))
    axs.set_ylim(0, 1)
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_numsample_{num_sample}_GlideUnDockedFFUnminimisedCount.png', bbox_inches='tight')


def plot_find_flipped_orientation_freq(wrapped_dataframe: SummaryDataFrameWrapper, centroi_distance_column_name: pd.Series, 
                                       rmsd_column_name: pd.Series, output_prefix: str):
    fig, axs = plt.subplots(1, figsize=(17,12))

    RMSD_THRESHOLD_FLIPPED_POCKET = 7
    LIGAND_CENTROID_DISTANCE_DIFF_DOCKED_POCKET = 8
    #protein_wrapped_df = wrapped_dataframe.filter_exact_feature(wrapped_dataframe.protein_name, [protein])
    protein_df = wrapped_dataframe.to_dataframe()

    protein_df[rmsd_column_name.name] = protein_df[rmsd_column_name.name].astype(float)
    protein_df[centroi_distance_column_name.name] = protein_df[centroi_distance_column_name.name].astype(float)

    #protein_df = protein_df[protein_df['MW / Da'] > 350]
    print(protein_df)
    grouped_df = protein_df.groupby([wrapped_dataframe.protein_name.name, 
                                     wrapped_dataframe.genai_model.name]).size().reset_index(name='counts')

    grouped_flipped_glide_df = protein_df[(protein_df[rmsd_column_name.name] > RMSD_THRESHOLD_FLIPPED_POCKET) &
                                          (protein_df[centroi_distance_column_name.name] > LIGAND_CENTROID_DISTANCE_DIFF_DOCKED_POCKET)] \
                                        .groupby([wrapped_dataframe.protein_name.name, 
                                                  wrapped_dataframe.genai_model.name]).size().reset_index(name='glide flipped counts')

    print(protein_df[['protein_name', 'ligand_id', 'centroid_distance_docking_glide', 'glide_redock-input_rms']].to_csv('tessst.csv'))
    grouped_undocked_glide_df = grouped_flipped_glide_df.merge(grouped_df)
    grouped_undocked_glide_df['flipped_percentage'] = grouped_undocked_glide_df['glide flipped counts'] / grouped_undocked_glide_df['counts']

    sns.barplot(data=grouped_undocked_glide_df,
                x=wrapped_dataframe.genai_model.name,
                y='flipped_percentage',
                hue=wrapped_dataframe.protein_name.name,
                ax=axs,
                palette=COLORBLIND_PALETTE)
    num_sample = list(set(wrapped_dataframe.num_sample))[0]
    axs.set_ylabel(f'flipped ligand proportion {output_prefix.replace("_", " ")}', fontsize=25)
    axs.set_xlabel(wrapped_dataframe.genai_model.name, fontsize=25)
    axs.tick_params(axis='both', which='major', labelsize=23)

    handles, labels = axs.get_legend_handles_labels()
    axs.get_legend().remove()
    
    fig.legend(handles[:6], labels[:6], loc="upper center", ncol=6, prop={'size': 25}, bbox_to_anchor=(0.5, 1.15))
    axs.set_ylim(0, 1)
    plt.tight_layout()
    plt.savefig(f'{output_prefix}_numsample_{num_sample}_GlideFlippedCount.png', bbox_inches='tight')


if __name__ == '__main__':
    summary_fname = '/opt/veincent/GenAI_manuscript/scramblebench/script/save_all_feature.parquet'

    is_plot_physicochemical_property = False
    is_plot_redocked_score = True
    is_find_best_molecule = False

    protein_family_dict = {'GPCR': ['5-HT2C', 'DRD2'],
                           'Kinase': ['CDK2', 'GSK3b'],
                           'Hydrolase':['3CLpro', 'NA']}
    protein_name_dict = {'5ht2c': '5-HT2C',
                        'drd2':'DRD2',
                        'cdk2':'CDK2',
                        'gsk3b':'GSK3b',
                        '3clpro': '3CLpro',
                        'na':'NA'}

    # will try to think a way to clean this, rn im too lazy zzzz
    control_redocking_score = [    [[-11.019, -12.511, 0.8755, 0.8442],
                                    [-10.504, -10.141, 0.7024, 0.6223]],
                                    [[-9.5779, -9.8834, 0.8814, 1.2211],
                                    [-6.8313, -8.6500, 0.9133, 1.4014]],
                                    [[-8.1221, -6.6220, 5.2166, 3.8787],
                                    [-8.4364, -7.0780, 0.2527, 0.3195]]]

    QED_THRESHOLD = 0.5
    DOCKING_SCORE_THRESHOLD = -8.0
    summary_path = Path(summary_fname)
    if summary_path.name.endswith('parquet'):
        summary_df = pd.read_parquet(summary_path)
    elif summary_path.name.endswith('csv'):
        summary_df = pd.read_csv(summary_path)
    else:
        TypeError(f'currently, we did not support any file besides parquet and csv. Please change your input {summary_path.name}')
        sys.exit(1)

    summary_df = summary_df.rename({'logP': 'log P',
                                    'SA_score': 'SA_Score',
                                    'MW' : 'MW / Da'}, axis='columns')
    wrapped_df = SummaryDataFrameWrapper(summary_df)
    wrapped_df.protein_name = wrapped_df.protein_name.replace(protein_name_dict)
    wrapped_df.protein_name = pd.Categorical(wrapped_df.protein_name, categories=set(wrapped_df.protein_name))
    

    
    for index, (protein_family, protein_list) in enumerate(protein_family_dict.items()):
        messy_summary_df = wrapped_df.to_dataframe()
        if is_plot_physicochemical_property:

            physicochemical_prop_list, wrapped_df_physicochem = wrapped_df.get_physicochemical_feature
            physicochem_wrapped_df_by_prot = wrapped_df_physicochem.filter_exact_feature(wrapped_df_physicochem.protein_name, protein_list)
            physicochem_wrapped_df_by_protlig = physicochem_wrapped_df_by_prot.filter_exact_feature(physicochem_wrapped_df_by_prot.stereo_id, [0, None])

            physicochem_wrapped_df_by_protlig.protein_name = pd.Categorical(physicochem_wrapped_df_by_protlig.protein_name, categories=protein_list, ordered=True)
            physicochem_wrapped_df_by_protlig.genai_model = pd.Categorical(physicochem_wrapped_df_by_protlig.genai_model, categories=GENAI_MODEL_ORDER, ordered=True)

            for unique_num_sample in set(physicochem_wrapped_df_by_protlig.num_sample):
                
                physicochem_wrapped_df_by_protlig_numsample = physicochem_wrapped_df_by_protlig.filter_exact_feature(physicochem_wrapped_df_by_protlig.num_sample, [unique_num_sample])
                if np.nan in set(physicochem_wrapped_df_by_protlig_numsample.genai_model):
                    continue
                
                plot_physicochemical_property(wrapped_dataframe=physicochem_wrapped_df_by_protlig_numsample,
                                              property_list=physicochemical_prop_list,
                                              output_prefix=f'{protein_family}')
       
        if is_plot_redocked_score:

            docking_prop_list, wrapped_df_docking = wrapped_df.get_docking_feature

            ff_unminimised_prop_list = [prop_list for prop_list in docking_prop_list if prop_list in [wrapped_df.FF_unminimised_glide_inplace_docking_score.name, 
                                                                                                      wrapped_df.FF_unminimised_glide_mininplace_docking_score.name,
                                                                                                      wrapped_df.FF_unminimised_vina_inplace_docking_score.name,
                                                                                                      wrapped_df.FF_minimised_glide_mininplace_docking_score.name]]
            docking_prop_list = list(set(docking_prop_list) - set(ff_unminimised_prop_list))

            docking_wrapped_df_by_prot = wrapped_df_docking.filter_exact_feature(wrapped_df_docking.protein_name, protein_list)

            docking_wrapped_df_by_prot.protein_name = pd.Categorical(docking_wrapped_df_by_prot.protein_name, categories=protein_list, ordered=True)
            docking_wrapped_df_by_prot.genai_model = pd.Categorical(docking_wrapped_df_by_prot.genai_model, categories=GENAI_MODEL_ORDER, ordered=True)
            docking_wrapped_df_by_prot.stereo_id = -1


            summary_generic_colname = docking_wrapped_df_by_prot.get_generic_feature_colname + [docking_wrapped_df_by_prot.molecular_weight.name]
            summary_ligand_id_colname = docking_wrapped_df_by_prot.ligand_id.name
            summary_protein_name_colname = docking_wrapped_df_by_prot.protein_name.name
            summary_generic_df = docking_wrapped_df_by_prot.to_dataframe()[summary_generic_colname].drop_duplicates()

            docking_df_by_protlig = pd.DataFrame()
            if isinstance(docking_wrapped_df_by_prot.redocked_vina_docking_score, pd.Series):

                sorted_redock_vina_df = docking_wrapped_df_by_prot.filter_stereo_id_based_on_value(
                                        column_list=[docking_wrapped_df_by_prot.redocked_vina_docking_score, 
                                                     docking_wrapped_df_by_prot.redocked_vina_rms, 
                                                     docking_wrapped_df_by_prot.FF_unminimised_vina_inplace_docking_score,
                                                     docking_wrapped_df_by_prot.redocked_vina_mol_block,
                                                     docking_wrapped_df_by_prot.input_mol_block
                                                     ], 
                                        sorted_column=docking_wrapped_df_by_prot.redocked_vina_docking_score)
            
                if isinstance(sorted_redock_vina_df, pd.DataFrame):
                    sorted_redock_vina_df = sorted_redock_vina_df.rename(columns={docking_wrapped_df_by_prot.input_mol_block.name:
                                                                                  INPUT_MOL_BLOCK_VINA_COLUMN})
                    if docking_df_by_protlig.empty:
                        docking_df_by_protlig = sorted_redock_vina_df
                    else:
                        docking_df_by_protlig = docking_df_by_protlig.merge(sorted_redock_vina_df, on=docking_wrapped_df_by_prot.get_generic_feature_colname)

            if isinstance(docking_wrapped_df_by_prot.redocked_glide_docking_score, pd.Series):

                sorted_redock_glide_df = docking_wrapped_df_by_prot.filter_stereo_id_based_on_value(
                                        column_list=[docking_wrapped_df_by_prot.redocked_glide_docking_score, 
                                                     docking_wrapped_df_by_prot.redocked_glide_rms, 
                                                     docking_wrapped_df_by_prot.FF_unminimised_glide_inplace_docking_score,
                                                     docking_wrapped_df_by_prot.FF_unminimised_glide_mininplace_docking_score,
                                                     docking_wrapped_df_by_prot.FF_minimised_glide_mininplace_docking_score,
                                                     docking_wrapped_df_by_prot.redocked_glide_mol_block,
                                                     docking_wrapped_df_by_prot.input_mol_block,
                                                     docking_wrapped_df_by_prot.molecular_weight], 
                                        sorted_column=docking_wrapped_df_by_prot.redocked_glide_docking_score)
                if isinstance(sorted_redock_glide_df, pd.DataFrame):
                    sorted_redock_glide_df = sorted_redock_glide_df.rename(columns={docking_wrapped_df_by_prot.input_mol_block.name:
                                                                                    INPUT_MOL_BLOCK_GLIDE_COLUMN})
                    if docking_df_by_protlig.empty:
                        docking_df_by_protlig = sorted_redock_glide_df
                    else:
                        docking_df_by_protlig = docking_df_by_protlig.merge(sorted_redock_glide_df, on=docking_wrapped_df_by_prot.get_generic_feature_colname)
            print(docking_df_by_protlig.columns)
            if not docking_df_by_protlig.empty:
                docking_wrapped_df_by_protlig = SummaryDataFrameWrapper(docking_df_by_protlig)
            else:
                raise TypeError(f'resulting dataframe after filtering docking parameter is not a dataframe, but {type(summary_generic_df)}')
            
            for unique_num_sample in set(docking_wrapped_df_by_protlig.num_sample):
                
                docking_wrapped_df_by_protlig_numsample = docking_wrapped_df_by_protlig.filter_exact_feature(docking_wrapped_df_by_protlig.num_sample, [unique_num_sample])
                if np.nan in set(docking_wrapped_df_by_protlig_numsample.genai_model):
                    continue
                docking_wrapped_df_by_protlig_numsample.calculate_redocking_centroid_distance(docking_wrapped_df_by_protlig_numsample.input_mol_block_glide,
                                                                                                    docking_wrapped_df_by_protlig_numsample.redocked_glide_mol_block)

                plot_find_flipped_orientation_freq(wrapped_dataframe=docking_wrapped_df_by_protlig_numsample, 
                                                   centroi_distance_column_name=docking_wrapped_df_by_protlig_numsample.centroid_distance_docking_glide,
                                                   rmsd_column_name=docking_wrapped_df_by_protlig_numsample.redocked_glide_rms,
                                                   output_prefix=f'{protein_family}_flipped_orientation')
                if unique_num_sample == 500 and protein_family == 'Kinase':
                    sys.exit(0)
                # plot_glide_undocked_freq(wrapped_dataframe=docking_wrapped_df_by_protlig_numsample, column_name = docking_wrapped_df_by_protlig_numsample.FF_unminimised_glide_inplace_docking_score,
                #                          output_prefix=f'{protein_family}_fully_unminimised')
                # plot_glide_undocked_freq(wrapped_dataframe=docking_wrapped_df_by_protlig_numsample, column_name = docking_wrapped_df_by_protlig_numsample.FF_unminimised_glide_mininplace_docking_score,
                #                          output_prefix=f'{protein_family}_unminimised_glidemininplace')
                # plot_glide_undocked_freq(wrapped_dataframe=docking_wrapped_df_by_protlig_numsample, column_name = docking_wrapped_df_by_protlig_numsample.FF_minimised_glide_mininplace_docking_score,
                #                          output_prefix=f'{protein_family}_fully_minimised')
                
                # plot_redocking_score(wrapped_dataframe=docking_wrapped_df_by_protlig_numsample,
                #                     property_list=docking_prop_list,
                #                     output_prefix=f'{protein_family}', reference_score=control_redocking_score[index])
        
        if is_find_best_molecule:
            combined_df = pd.DataFrame()
            best_prop_list, best_mol_wrapped_df = wrapped_df.get_best_molecule_feature

            for protein_name in protein_list:
                combined_df = pd.DataFrame()
                for genai_model in set(best_mol_wrapped_df.genai_model):
                    for unique_num_sample in set(best_mol_wrapped_df.num_sample):

                        best_mol_wrapped_df_by_prot = best_mol_wrapped_df.filter_exact_feature(best_mol_wrapped_df.protein_name, [protein_name])
                        best_mol_wrapped_df_by_prot_model = best_mol_wrapped_df_by_prot.filter_exact_feature(best_mol_wrapped_df_by_prot.genai_model, [genai_model])
                        best_mol_wrapped_df_by_prot_model_numsample = best_mol_wrapped_df_by_prot_model.filter_exact_feature(best_mol_wrapped_df_by_prot_model.num_sample, [unique_num_sample])
                        
                        if best_mol_wrapped_df_by_prot_model_numsample.to_dataframe().empty:
                            continue

                        filtered_best_mol_wrapped_df = best_mol_wrapped_df_by_prot_model_numsample.copy()
                        if isinstance(filtered_best_mol_wrapped_df.glide_phase_score, pd.Series):
                            glide_phase_colname = filtered_best_mol_wrapped_df.glide_phase_score.name
                            best_mol_with_phase_score_df = filtered_best_mol_wrapped_df.to_dataframe().dropna(subset=[glide_phase_colname]).copy()
                            if not best_mol_with_phase_score_df.empty:
                                filtered_best_mol_wrapped_df = SummaryDataFrameWrapper(best_mol_with_phase_score_df)    

                        if isinstance(filtered_best_mol_wrapped_df.qed, pd.Series):
                            filtered_best_mol_df = filtered_best_mol_wrapped_df.to_dataframe()
                            best_mol_with_qed_df = filtered_best_mol_df[filtered_best_mol_df[filtered_best_mol_wrapped_df.qed.name] >= QED_THRESHOLD]
                            best_mol_with_best_qed_df = filtered_best_mol_df[filtered_best_mol_df[filtered_best_mol_wrapped_df.qed.name] == filtered_best_mol_wrapped_df.qed.max()]
                            if not best_mol_with_qed_df.empty:
                                filtered_best_mol_wrapped_df = SummaryDataFrameWrapper(best_mol_with_qed_df)

                        if isinstance(filtered_best_mol_wrapped_df.redocked_glide_docking_score, pd.Series):
                            filtered_best_mol_df = filtered_best_mol_wrapped_df.to_dataframe()
                            best_mol_with_glide_df = filtered_best_mol_df[filtered_best_mol_df[filtered_best_mol_wrapped_df.redocked_glide_docking_score.name] <= DOCKING_SCORE_THRESHOLD]
                            best_mol_with_best_glide_df = filtered_best_mol_df[filtered_best_mol_df[filtered_best_mol_wrapped_df.redocked_glide_docking_score.name] == filtered_best_mol_wrapped_df.redocked_glide_docking_score.min()]
                            if not best_mol_with_glide_df.empty:
                                filtered_best_mol_df = best_mol_with_glide_df
                            else:
                                filtered_best_mol_df = best_mol_with_best_glide_df


                        elif isinstance(filtered_best_mol_wrapped_df.redocked_vina_docking_score, pd.Series):
                            filtered_best_mol_df = filtered_best_mol_wrapped_df.to_dataframe()
                            best_mol_with_vina_df = filtered_best_mol_df[filtered_best_mol_df[filtered_best_mol_wrapped_df.redocked_vina_docking_score.name] <= DOCKING_SCORE_THRESHOLD]
                            best_mol_with_best_vina_df = filtered_best_mol_df[filtered_best_mol_df[filtered_best_mol_wrapped_df.redocked_vina_docking_score.name] == filtered_best_mol_wrapped_df.redocked_vina_docking_score.min()]

                            if not best_mol_with_vina_df.empty:
                                filtered_best_mol_df = best_mol_with_vina_df 
                            else:
                                filtered_best_mol_df = best_mol_with_best_vina_df
                        
                        else:
                            filtered_best_mol_df = best_mol_with_best_qed_df

                        combined_df = pd.concat([combined_df, filtered_best_mol_df])

                mol_l = [wrapped_df.fetch_mol(row_molecule) for index, row_molecule in combined_df.iterrows()]
                combined_df.to_csv(f'final_result4_{protein_name}.csv')
                with Chem.SDWriter(f'final_compound4_{protein_name}.sdf') as writer:
                    for mol in mol_l:
                        writer.write(mol)