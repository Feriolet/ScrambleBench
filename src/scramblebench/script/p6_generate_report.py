from typing import Any, Callable
from pathlib import Path

import yaml
import argparse
import logging
import sys
import subprocess
import json
from scramblebench.script.config_preparation import config_constant, config_report, config_parameter
from scramblebench.script.utils.error_handler import DirNotFoundError
from scramblebench.script.utils.process_data import read_input
from scramblebench.script.utils.process_mol import calculate_rms
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
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak
from reportlab.lib.styles import getSampleStyleSheet
from reportlab.lib.units import inch
from reportlab.lib.enums import TA_CENTER, TA_RIGHT, TA_JUSTIFY

import importlib.metadata
from datetime import datetime

COLORBLIND_PALETTE = ["#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#999999"]

def plot_violin(dataframe, columns, hue='Model', x='protein_name'):
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


def add_space(num):
    return '&nbsp;' * num

def generate_data(yaml_file):

    input_csv=Path(yaml_file).parent / 'all.csv'

    config_data = yaml.safe_load(open(yaml_file, 'r'))
    report_data = config_report.ReportConfig(config_data)
    parameter_data = config_parameter.ParameterConfig(config_data)

    df = pd.read_csv(input_csv)

    doc = SimpleDocTemplate(str(Path(yaml_file).parent / 'hi.pdf'), pagesize=letter,
                            rightMargin=inch/2, leftMargin=inch/2,
                            topMargin=40, bottomMargin=18)

    Story = []
    image_buffer_list = []

    styles = getSampleStyleSheet()
    styles['Heading1'].alignment = TA_CENTER
    styles['Heading2'].alignment = TA_CENTER
    styles['Normal'].leading = 16

    Story.append(Paragraph(f"<font size= 18>ScrambleBench Report v{importlib.metadata.version('scramblebench')}</font>", styles['Heading1']))
    Story.append(Paragraph(f"<font size= 12>Generation Date: {datetime.now().strftime('%d %B %Y')}</font>", styles['Normal']))
    Story.append(Spacer(1,4))
    parameter_txt = f"<b>Parameter</b> <br/>\
                     Protein Name {add_space(16)}: {parameter_data.protein_value}<br/>\
                     Model list  {add_space(23)}: {', '.join(parameter_data.model_list_value)} <br/>\
                     Requested Num Sample: {parameter_data.num_sample_value}<br/>"
    
    if parameter_data.batch_parameter_dict:
        parameter_txt += 'Other parameters: <br/>'
        parameter_txt += '<br/>'.join([str(key) + ' : ' + str(val) for key, val in parameter_data.batch_parameter_dict.items()])

    Story.append(Paragraph(f"<font size= 12>{parameter_txt}</font>", styles['Normal']))
    Story.append(Spacer(1,12))

    if report_data.docking_score_value:
        Story.append(PageBreak())
        redocking_score_columns = [col for col in df.columns if 'redocking_score' in col.lower()]
        
        image_buffer = plot_violin(dataframe=df, columns=redocking_score_columns)
        image_buffer_list.append(image_buffer)

        Story.append(Paragraph(f"<font size= 18>Redocking Score</font>", styles['Heading1']))
        Story.append(Image(image_buffer, 7*inch, 7*inch))

    if report_data.rmsd_value:
        Story.append(PageBreak())
        rmsd_columns = [col for col in df.columns if 'rmsd' in col.lower()]
        
        image_buffer = plot_violin(dataframe=df, columns=rmsd_columns)
        image_buffer_list.append(image_buffer)

        Story.append(Paragraph(f"<font size= 18>RMSD Score</font>", styles['Heading1']))
        Story.append(Image(image_buffer, 7*inch, 7*inch))

    if report_data.qed_value:
        Story.append(PageBreak())
        qed_column = [col for col in df.columns if 'qed' in col.lower()]
        
        image_buffer = plot_violin(dataframe=df, columns=qed_column)
        image_buffer_list.append(image_buffer)

        Story.append(Paragraph(f"<font size= 18>QED</font>", styles['Heading1']))
        Story.append(Image(image_buffer, 7*inch, 7*inch))

    doc.build(Story)

    for image_buffer in image_buffer_list:
        image_buffer.close()

    
    
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate report from p5_collect_data_analysis.py")

    parser.add_argument("-i", "--input", help="config yaml input file or txt file containing yaml filepath", required=True, type=str)
    args = parser.parse_args()

    logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='%(asctime)s - %(module)s: - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
    
    logging.info('Running p5_collect_data_analysis.py')
    logging.info('Reading the config filename :)')

    yaml_file_list = read_input(args.input)

    output_analysis_fname = Path(args.input).parent / 'compiled_result.csv'
    output_df = pd.DataFrame()

    for yaml_file in yaml_file_list:
        generate_data(yaml_file)

