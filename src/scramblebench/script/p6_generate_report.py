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


class ScrambleBenchPlotReport:

    def __init__(self, summary_df, columns, title_caption=None):
        
        self.summary_df = summary_df
        self.columns = columns
        self.title_caption = title_caption

        self.image_buffer_list = []
        self.title_font_size = 18
        self.styles = getSampleStyleSheet()
        self.styles['Heading1'].alignment = TA_CENTER
        self.styles['Normal'].leading = 16

    def render(self):
        
        Story = []

        Story.append(PageBreak())
                       
        image_buffer = self.plot_violin(dataframe=self.summary_df, 
                                        columns=self.columns)
        
        self.image_buffer_list.append(image_buffer)

        if self.title_caption:
            Story.append(Paragraph(f"<font size= {self.title_font_size}>{self.title_caption}</font>", self.styles['Heading1']))
        
        Story.append(Image(image_buffer, 7*inch, 7*inch))

        return Story


    def plot_violin(self, dataframe, columns, hue='Model', x='protein_name'):
        fig, axs = plt.subplots(len(columns), 1, figsize=(7,7), sharex=True)

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

    input_csv=Path(yaml_file).parent / 'summary.csv'
    df = pd.read_csv(input_csv)

    config_data = yaml.safe_load(open(yaml_file, 'r'))

    report_data = config_report.ReportConfig(config_data)
    parameter_data = config_parameter.ParameterConfig(config_data)
 
    report = ScrambleBenchReport(yaml_file)
    report.add_section(ScrambleBenchParameterReport(parameter_data=parameter_data))


    if config_constant.INPUT_KEY in config_data:
        input_data = config_input.InputStructure(config_data[config_constant.INPUT_KEY])
        report.add_section(ScrambleBenchInputReport(input_data=input_data))

    if report_data.docking_score_value:
        redocking_score_columns = [col for col in df.columns if 'redocking_score' in col.lower()]
        report.add_section(ScrambleBenchPlotReport(summary_df=df,
                                                columns=redocking_score_columns,
                                                title_caption='Redocking Score'))

    if report_data.rmsd_value:
        rmsd_columns = [col for col in df.columns if 'rmsd' in col.lower()]
        report.add_section(ScrambleBenchPlotReport(summary_df=df,
                                                columns=rmsd_columns,
                                                title_caption='RMSD Score'))

    if report_data.qed_value:
        qed_column = [col for col in df.columns if 'qed' in col.lower()]
        report.add_section(ScrambleBenchPlotReport(summary_df=df,
                                                columns=qed_column,
                                                title_caption='QED'))
    
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

