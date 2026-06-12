from typing import Any, Callable
from pathlib import Path

import yaml
import argparse
import logging
import sys
import subprocess
import json
from scramblebench.script.config_preparation import config_constant, config_input, config_report, config_parameter
from scramblebench.script.utils.error_handler import DirNotFoundError
from scramblebench.script.utils.process_data import read_input
from scramblebench.script.utils.process_mol import calculate_rms, calculate_physicochemical_properties, MOL_PROPERTY_CALCULATORS
import os
import pandas as pd
import tempfile
from copy import deepcopy
from collections import defaultdict
import rdkit
from rdkit import Chem
from enum import Enum, IntEnum
import seaborn as sns
import matplotlib.pyplot as plt
import numpy as np
import io
import math


from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib.enums import TA_CENTER, TA_RIGHT, TA_JUSTIFY

import importlib.metadata
from datetime import datetime

COLORBLIND_PALETTE = ["#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#999999"]


def add_space(num):
    return '&nbsp;' * num


class ScrambleBenchParameterReport:

    def __init__(self, parameter_data):
        self.parameter_data = parameter_data
        self.title_font_size = 18
        self.parameter_font_size = 12

        self.styles = getSampleStyleSheet()
        self.styles['Heading1'].alignment = TA_CENTER
        self.styles['Normal'].leading = 16


    def render(self):

        Story = []
        Story.append(Paragraph(f"<font size={self.title_font_size}>ScrambleBench Report v{importlib.metadata.version('scramblebench')}</font>", 
                               self.styles['Heading1']))
        Story.append(Paragraph(f"<font size={self.parameter_font_size}>Generation Date: {datetime.now().strftime('%d %B %Y')}</font>", 
                               self.styles['Normal']))
        Story.append(Spacer(1,4))
        parameter_txt = f"<b>Parameter</b> <br/>\
                        Protein Name {add_space(16)}: {self.parameter_data.protein_value}<br/>\
                        Model list  {add_space(23)}: {', '.join(self.parameter_data.model_list_value)} <br/>\
                        Requested Num Sample: {self.parameter_data.num_sample_value}<br/>"
        Story.append(Paragraph(f"<font size= 12>{parameter_txt}</font>", self.styles['Normal']))
        Story.append(Spacer(1,4))
        if self.parameter_data.batch_parameter_dict:
            batch_parameter_txt = 'Other parameters: <br/>'
            batch_parameter_txt += '<br/>'.join([str(key) + ' : ' + str(val) for key, val in self.parameter_data.batch_parameter_dict.items()])

        Story.append(Paragraph(f"<font size= 12>{batch_parameter_txt}</font>", self.styles['Normal']))
        Story.append(Spacer(1,12))

        return Story


class ScrambleBenchInputReport:

    def __init__(self, input_data: config_input.InputStructure):
        self.input_data = input_data
        self.font_size = 15
        self.styles = getSampleStyleSheet()
        self.styles['Heading2'].alignment = TA_CENTER
        self.styles['Normal'].leading = 16
        self.residue_font_size = 15
        self.residue_num_font_size = 8
        self.residue_caption_font_size = 12

        self.image_buffer_list = []
        self.mono_style = ParagraphStyle(name='MonoStyle',
                                        parent=self.styles['Normal'],
                                        fontName='Courier',
                                        fontSize=12,
                                        leading=14)

    def generate_pocket_paragraph(self, pocket_list):

        text = ""

        for pocket_data in pocket_list:
            text += f"<font size={self.residue_font_size}>"
            new_residue_block = True

            last_residue = pocket_data[-1]
            resnum_digit = len(str(last_residue[0]))

            for count, (res_id, resname, is_pocket) in enumerate(pocket_data, start=1):

                if new_residue_block:

                    padded_zero = '0' * (len(f"{res_id:0{resnum_digit}d}") - len(str(res_id)))
                    res_id_text = f"<font color='white'>{padded_zero}</font>{res_id}"
                    text += f"</font><font size={self.residue_num_font_size}><super>{res_id_text}</super></font><font size={self.residue_font_size}>"
                    
                    new_residue_block = False

                if is_pocket:
                    text += f"<font color='red'>{resname}</font>"
                else:
                    text += resname

                if count % 10 == 0:
                    text += ' '
                    new_residue_block = True

                if count % 40 == 0:
                    text += "<br/>"

            text += '</font>'
            text += "<br/><br/>"
  
        return text
    

    def generate_ligand_paragraph(self):
        mol = self.input_data.get_sdf_mol()
        mol.RemoveAllConformers()
        
        draw_mol = Chem.RemoveHs(mol)
        Chem.rdDepictor.Compute2DCoords(draw_mol)
        drawer = Chem.Draw.rdMolDraw2D.MolDraw2DCairo(500, 500)
        Chem.Draw.rdMolDraw2D.PrepareAndDrawMolecule(drawer, draw_mol)
        drawer.FinishDrawing()

        png_data = drawer.GetDrawingText()
        buffer = io.BytesIO(png_data)
        buffer.seek(0)
        self.image_buffer_list.append(buffer)

        mol = calculate_physicochemical_properties(mol)

        smi = Chem.MolToSmiles(Chem.RemoveHs(mol))

        smi_text = ''
        for count, character in enumerate(smi, start=1):
            if count % 30 == 0:
                smi_text += '\n'
            smi_text += character

        mol_data = Table([[mol_prop, 
                           f'{float(mol.GetProp(mol_prop)):.3f}'] for mol_prop in MOL_PROPERTY_CALCULATORS.keys()] \
                            + [['SMILES', smi_text]])
        
        mol_data.setStyle(TableStyle([('VALIGN', (0, 0), (-1, -1), 'TOP'),
                                   ('FONTNAME', (0, 0), (-1, -1), 'Courier')]))

        table = Table([[Image(buffer, 3*inch, 3*inch, hAlign='LEFT'), mol_data]],
                      rowHeights=[3*inch])
        table.setStyle(TableStyle([('VALIGN', (0, 0), (-1, -1), 'MIDDLE')]))

        return table

    def render(self):

        Story = []

        Story.append(Paragraph(f"<font size= {self.font_size}>Protein Information</font>", self.styles['Heading2']))
        Story.append(Spacer(1,5))
        pocket_information = self.input_data.fetch_residue_is_pocket()
        pocket_txt = self.generate_pocket_paragraph(pocket_information)
        Story.append(Paragraph(pocket_txt, self.mono_style))
        Story.append(Spacer(1,10))
        Story.append(Paragraph(f"<font size= {self.residue_caption_font_size}>Pocket Residues within 10 A of ligands are colored in red</font>", self.styles['Normal']))
        Story.append(Spacer(1,10))

        Story.append(Paragraph(f"<font size= {self.font_size}>Ligand Information</font>", self.styles['Heading2']))

        mol_story = self.generate_ligand_paragraph()
        Story.append(mol_story)

        return Story 
    

    def cleanup(self):

        for img in self.image_buffer_list:
            img.close()


class ScrambleBenchSummaryJSONReport:

    def __init__(self, input_data):
        self.json_data = input_data
        self.title_font_size = 18
        self.parameter_font_size = 12
        self.post_processing = False
        self.diversity = False
        self.genbench3d = False

        self.styles = getSampleStyleSheet()
        self.styles['Heading1'].alignment = TA_CENTER
        self.styles['Normal'].leading = 16


    def detect_analysis_type(self):
        for value in self.json_data.values():
            scramblebench_json_keys = value.keys()
            if 'total_mol' in scramblebench_json_keys and 'uniqueness' in scramblebench_json_keys:
                self.post_processing = True
            if 'unminimised_Validity3D' in scramblebench_json_keys or 'minimised_Validity3D' in scramblebench_json_keys:
                self.genbench3d = True
            if 'diversity' in scramblebench_json_keys:
                self.diversity = True
    

    def bold_table(self, df, index, filter='max'):
        if not isinstance(index, list):
            index = [index]
        
        for metric in index:
            if filter == 'max':
                value = df.loc[metric].max()
            elif filter == 'min':
                value = df.loc[metric].min()
            
            df.loc[metric, df.loc[metric] == value] = Paragraph(f"<b>{value}</b>")

        return df
    
    def generate_summary_txt(self):

        diversity_df = pd.DataFrame()
        non_diversity_dict = defaultdict(dict)
        for model, data in self.json_data.items():
            df = pd.DataFrame(data['diversity'])
            df['model'] = model

            if diversity_df.empty:
                diversity_df = df
            else:
                diversity_df = pd.concat([diversity_df, df])

            data.pop('diversity', None)
            non_diversity_dict[model] = data
        
        non_diversity_df = pd.DataFrame(non_diversity_dict).round(2).astype(object)


        non_diversity_df = self.bold_table(df=non_diversity_df,
                                           index=non_diversity_df.index.tolist(),
                                           filter='max')
        
        non_diversity_df = non_diversity_df.reset_index(names='Metric')

        diversity_df_list = []
        for method, data in diversity_df.groupby('method'):
            df = data.drop(columns=['method']).set_index('model').T.round(2).astype(object)
            df = self.bold_table(df=df,
                            index='score',
                            filter='max')  
            
            df = df.reset_index(names='Metric')
            diversity_df_list.append((method, df))
           
        return diversity_df_list, non_diversity_df

    def render(self):

        Story = []
        Story.append(PageBreak())

        Story.append(Paragraph(f"<font size={self.title_font_size}>ScrambleBench Analysis</font>", 
                               self.styles['Heading1']))
        
        Story.append(Spacer(1,6))
        self.detect_analysis_type()

        parameter_txt = f"<b>Ran Analysis</b> <br/>\
                        PostProcessing {add_space(16)}: {self.post_processing}<br/>\
                        GenBench3D  {add_space(19)}: {self.genbench3d} <br/>\
                        Diversity {add_space(27)}: {self.diversity}<br/>"
        
        Story.append(Paragraph(f"<font size= 12>{parameter_txt}</font>", self.styles['Normal']))
        Story.append(Spacer(1,20))

        diversity_df_list, non_diversity_df = self.generate_summary_txt()

        Story.append(Paragraph(f"<font size= 12>General Analysis Result </font>", self.styles['Normal']))
        Story.append(Table([non_diversity_df.columns[:,].values.astype(str).tolist()] + non_diversity_df.values.tolist(), hAlign='LEFT') )
        
        Story.append(Spacer(1,10))
        
        Story.append(Paragraph(f"<font size= 12>Diversity Analysis Result </font>", self.styles['Normal']))
        Story.append(Spacer(1,6))
        for diversity_tuple in diversity_df_list:
            Story.append(Paragraph(f"<font size= 10><b>{diversity_tuple[0]}</b></font>", self.styles['Normal']))
            diversity_df = diversity_tuple[-1]
            Story.append(Table([diversity_df.columns[:,].values.astype(str).tolist()] + diversity_df.values.tolist(), hAlign='LEFT'))
            Story.append(Spacer(1,6))
        # Story.append(Paragraph(f"<font size= 12>{batch_parameter_txt}</font>", self.styles['Normal']))
        # Story.append(Spacer(1,12))

        return Story
    

class ScrambleBenchPlotReport:

    def __init__(self, summary_df, columns, type='violin', title_caption=None):
        
        self.summary_df = summary_df
        self.columns = columns
        self.title_caption = title_caption
        self.plot_type = type

        self.image_buffer_list = []
        self.title_font_size = 18
        self.styles = getSampleStyleSheet()
        self.styles['Heading1'].alignment = TA_CENTER
        self.styles['Normal'].leading = 16

    def render(self):
        
        Story = []

        Story.append(PageBreak())

        if self.plot_type == 'violin':
            image_buffer = self.plot_violin(dataframe=self.summary_df, 
                                            columns=self.columns)
        elif self.plot_type == 'box':
            image_buffer = self.plot_box(dataframe=self.summary_df, 
                                            columns=self.columns)   
        elif self.plot_type == 'raincloud':
            image_buffer = self.plot_raincloud(dataframe=self.summary_df, 
                                            columns=self.columns)           
        else:
            logging.warning(f'Plot type {self.plot_type} is not recognised. Defaulting to violin plot.')
            image_buffer = self.plot_violin(dataframe=self.summary_df, 
                                            columns=self.columns)
                    
        self.image_buffer_list.append(image_buffer)

        if self.title_caption:
            Story.append(Paragraph(f"<font size= {self.title_font_size}>{self.title_caption}</font>", self.styles['Heading1']))
        
        Story.append(Image(image_buffer, 7*inch, 7*inch))

        return Story


    def plot_violin(self, dataframe, columns, hue='Model', x='protein_name'):

        if len(columns) < 3:
            fig, axs = plt.subplots(len(columns), 1, figsize=(7,7), sharex=True)

        else:
            fig, axs = plt.subplots(math.ceil(len(columns) / 2), 2, figsize=(12, 12), sharex=True)

        for column, ax in zip(columns, np.ravel(axs)[:len(columns)]):
            sns.violinplot(data=dataframe,
                                x=x,
                                y=column,
                                ax=ax,
                                palette=COLORBLIND_PALETTE,
                                hue=hue,
                                width=0.5,
                                bw_method= 0.2,
                                orient= 'v', 
                                dodge= True,
                                linewidth=2)
                
            handles, labels = ax.get_legend_handles_labels()
            ax.get_legend().remove()

            ax.set_xlabel('Protein Name', fontsize=12)
            ax.set_ylabel(' '.join([word.capitalize() for word in column.split('_')]), fontsize=12)
            ax.tick_params(axis='both', which='major', labelsize=12)

        plt.subplots_adjust(wspace=0.1, hspace=0.05)

        buf = io.BytesIO()
        plt.tight_layout(rect=[0, 0, 1, 0.93]) 

        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.98), ncol=len(handles), prop={'size': 12})

        plt.savefig(buf, format = 'png', dpi=300)
        buf.seek(0)
        plt.close()

        return buf
    

    def plot_raincloud(self, dataframe, columns, hue='Model', x='protein_name'):

        import ptitprince as pt
        if len(columns) < 3:
            fig, axs = plt.subplots(len(columns), 1, figsize=(7,7), sharex=True)

        else:
            fig, axs = plt.subplots(math.ceil(len(columns) / 2), 2, figsize=(12, 12), sharex=True)

        for column, ax in zip(columns, np.ravel(axs)[:len(columns)]):

            pt.RainCloud(hue=hue, 
                        y=column, 
                        x=x, 
                        palette=COLORBLIND_PALETTE,  
                        data=dataframe, 
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

            ax.set_xlabel('Protein Name', fontsize=12)
            ax.set_ylabel(' '.join([word.capitalize() for word in column.split('_')]), fontsize=12)
            ax.tick_params(axis='both', which='major', labelsize=12)

        plt.subplots_adjust(wspace=0.1, hspace=0.05)

        buf = io.BytesIO()
        plt.tight_layout(rect=[0, 0, 1, 0.93]) 

        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.98), ncol=len(handles), prop={'size': 12})

        plt.savefig(buf, format = 'png', dpi=300)
        buf.seek(0)
        plt.close()

        return buf



    def plot_box(self, dataframe, columns, hue='Model', x='protein_name'):
        if len(columns) < 3:
            fig, axs = plt.subplots(len(columns), 1, figsize=(7,7), sharex=True)

        else:
            fig, axs = plt.subplots(math.ceil(len(columns) / 2), 2, figsize=(12, 12), sharex=True)

        for column, ax in zip(columns, np.ravel(axs)[:len(columns)]):
            sns.boxplot(data=dataframe,
                                x=x,
                                y=column,
                                ax=ax,
                                palette=COLORBLIND_PALETTE,
                                hue=hue,
                                width=0.5,
                                orient= 'v', 
                                dodge= True,
                                linewidth=2)
                
            handles, labels = ax.get_legend_handles_labels()
            ax.get_legend().remove()

            ax.set_xlabel('Protein Name', fontsize=12)
            ax.set_ylabel(' '.join([word.capitalize() for word in column.split('_')]), fontsize=12)
            ax.tick_params(axis='both', which='major', labelsize=12)

        plt.subplots_adjust(wspace=0.1, hspace=0.05)

        buf = io.BytesIO()
        plt.tight_layout(rect=[0, 0, 1, 0.93]) 

        fig.legend(handles, labels, loc="upper center", bbox_to_anchor=(0.5, 0.98), ncol=len(handles), prop={'size': 12})

        plt.savefig(buf, format = 'png', dpi=300)
        buf.seek(0)
        plt.close()

        return buf
    

    def cleanup(self):
        for img in self.image_buffer_list:
            img.close()


class ScrambleBenchReport:

    def __init__(self, yaml_file):
        self.yaml_file = yaml_file
        self.section = []

    def add_section(self, section):
        self.section.append(section)

    def generate_report(self):

        doc = SimpleDocTemplate(str(Path(self.yaml_file).parent / 'report.pdf'), pagesize=letter,
                                rightMargin=inch/2, leftMargin=inch/2,
                                topMargin=40, bottomMargin=18)


        Story = []

        try:
            for section in self.section:
                Story += section.render()

            doc.build(Story)

        finally:

            for section in self.section:
                cleanup_func = getattr(section, "cleanup", None)
                if callable(cleanup_func):
                    cleanup_func()


def generate_report(yaml_file):

    input_csv = Path(yaml_file).parent / 'summary.csv'
    input_json = Path(yaml_file).parent / 'summary.json'

    json_data = None
    with open(input_json, 'r') as json_f:
        json_data = json.load(json_f)

    df = pd.read_csv(input_csv)

    config_data = yaml.safe_load(open(yaml_file, 'r'))

    report_data = config_report.ReportConfig(config_data)
    parameter_data = config_parameter.ParameterConfig(config_data)
 
    report = ScrambleBenchReport(yaml_file)
    report.add_section(ScrambleBenchParameterReport(parameter_data=parameter_data))


    if config_constant.INPUT_KEY in config_data:
        input_data = config_input.InputStructure(config_data[config_constant.INPUT_KEY])
        report.add_section(ScrambleBenchInputReport(input_data=input_data))

    if json_data:
        report.add_section(ScrambleBenchSummaryJSONReport(input_data=json_data))

    if report_data.docking_score_value:
        redocking_score_columns = [col for col in df.columns if 'redocking_score' in col.lower()]
        report.add_section(ScrambleBenchPlotReport(summary_df=df,
                                                columns=redocking_score_columns,
                                                type=report_data.plot_value,
                                                title_caption='Redocking Score'))

    if report_data.rmsd_value:
        rmsd_columns = [col for col in df.columns if 'rmsd' in col.lower()]
        report.add_section(ScrambleBenchPlotReport(summary_df=df,
                                                columns=rmsd_columns,
                                                type=report_data.plot_value,
                                                title_caption='RMSD Score'))

    if report_data.qed_value:
        qed_column = [col for col in df.columns if 'qed' in col.lower()]
        report.add_section(ScrambleBenchPlotReport(summary_df=df,
                                                columns=qed_column,
                                                type=report_data.plot_value,
                                                title_caption='QED'))
    
    if report_data.validity3d_value:
        validity3d_column = [col for col in df.columns if 'validity3d' in col.lower()]
        report.add_section(ScrambleBenchPlotReport(summary_df=df,
                                                columns=validity3d_column,
                                                type=report_data.plot_value,
                                                title_caption='Validity3D Metric'))

    report.generate_report()

    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate report from p5_collect_data_analysis.py")

    parser.add_argument("-i", "--input", help="config yaml input file or txt file containing yaml filepath", required=True, type=str)
    args = parser.parse_args()

    logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='%(asctime)s - %(module)s: - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
    
    logging.info('Running p6_generate_report.py')
    logging.info('Reading the config filename :)')

    yaml_file_list = read_input(args.input)

    output_analysis_fname = Path(args.input).parent / 'compiled_result.csv'
    output_df = pd.DataFrame()

    for yaml_file in yaml_file_list:
        generate_report(yaml_file)

