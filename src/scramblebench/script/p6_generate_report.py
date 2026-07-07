"generate report as the final step of ScrambleBench"
import argparse
import logging
import sys
import json
import io
import math
import importlib.metadata

from typing import Any
from pathlib import Path
from collections import defaultdict
from datetime import datetime

import yaml
import seaborn as sns
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import rdkit
from rdkit import Chem
from matplotlib.axes import Axes


from reportlab.lib import colors
from reportlab.lib.pagesizes import letter
from reportlab.platypus import SimpleDocTemplate, Paragraph, Spacer, Image, PageBreak, Table, TableStyle
from reportlab.lib.styles import getSampleStyleSheet, ParagraphStyle
from reportlab.lib.units import inch
from reportlab.lib.enums import TA_CENTER

from scramblebench.script.config_preparation import config_constant, config_input, config_report, config_parameter
from scramblebench.script.utils.process_data import read_input
from scramblebench.script.utils.process_mol import  calculate_physicochemical_properties, MOL_PROPERTY_CALCULATORS


COLORBLIND_PALETTE = ["#E69F00","#56B4E9","#009E73","#F0E442","#0072B2","#D55E00","#CC79A7","#999999"]


def add_space(num: int) -> str:
    """generate space for reportslab. Typical whitespace such as \n or \t does not work

    Args:
        num (int): number of space

    Returns:
        str: reportslab' readable spacing
    """
    return '&nbsp;' * num


class ScrambleBenchParameterReport:
    """class containing the header parameter report section
    """

    def __init__(self, parameter_data: config_parameter.ParameterConfig):
        """initialize class

        Args:
            parameter_data (config_parameter.ParameterConfig): ParameterConfig
        """
        self.parameter_data = parameter_data
        self.title_font_size = 18
        self.subtitle_font_size = 15
        self.parameter_font_size = 12

        self.styles = getSampleStyleSheet()
        self.styles['Heading1'].alignment = TA_CENTER
        self.styles['Normal'].leading = 16

        self.table_style = TableStyle([
            # Header styling
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#2C3E50')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 12),
            ('LEFTPADDING', (0, 0), (-1, 0), 65),
            ('RIGHTPADDING', (0, 0), (-1, 0), 65),
            # Body styling
            ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
            ('FONTSIZE', (0, 1), (-1, -1), 10),
            ('ALIGN', (1, 0), (-1, -1), 'CENTER'),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),

            # Alternating row colors
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.HexColor('#F8F9F9'), colors.white]),

            # Borders
            #('LINEBELOW', (0, 0), (-1, 0), 2, colors.HexColor('#2C3E50')),
            ('GRID', (0, 1), (-1, -1), 0.5, colors.HexColor('#BDC3C7')),
            ('TOPPADDING', (0, 1), (-1, -1), 10),
            ('BOTTOMPADDING', (0, 1), (-1, -1), 10),
        ])


    def render(self) -> list[Any]:
        """render the parameter section as a story

        Returns:
            list[Any]: list of reportslab story consisting of sections of paragraphs/spacer
        """

        story = []
        story.append(Paragraph(f"<font size={self.title_font_size}>ScrambleBench Report " + \
                               f"v{importlib.metadata.version('scramblebench')}</font>",
                               self.styles['Heading1']))
        story.append(Paragraph(f"<font size={self.parameter_font_size} color='#808080'>Generation Date: " + \
                               f"{datetime.now().strftime('%d %B %Y')}</font>",
                               self.styles['Normal']))
        story.append(Spacer(1,6))


        parameter_table = Table([['Parameter Name', 'Parameter Value'],
                                 ['Protein Name', self.parameter_data.protein_value],
                                 ['Model List', ', '.join(self.parameter_data.model_list_value)],
                                  ['Requested num_sample', self.parameter_data.num_sample_value]],
                                  hAlign='LEFT')
        parameter_table.setStyle(self.table_style)
        story.append(Paragraph(f'<font size={self.parameter_font_size}><b>Parameter Used:</b></font>', self.styles['Normal']))
        story.append(Spacer(1,2))
        story.append(parameter_table)

        story.append(Spacer(1,4))
        if self.parameter_data.batch_parameter_dict:

            batch_parameter_txt = 'Other parameters: <br/>'
            for key, val in self.parameter_data.batch_parameter_dict.items():
                if key != 'num_sample':
                    batch_parameter_txt += str(key) + ' : ' + str(val) + '<br/>'
            # batch_parameter_txt += '<br/>'.join([str(key) + ' : ' + str(val) \
            #                                      for key, val in self.parameter_data.batch_parameter_dict.items()])

        story.append(Paragraph(f"<font size= 12>{batch_parameter_txt}</font>", self.styles['Normal']))
        story.append(Spacer(1,12))

        return story


class ScrambleBenchInputReport:
    """class containing the protein/ligand input section
    """

    def __init__(self, input_data: config_input.InputConfig):
        """initialize class

        Args:
            input_data (config_input.InputConfig): InputConfig
        """
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

        self.title_spacer = Spacer(1, 5)
        self.section_spacer = Spacer(1, 10)

    def generate_pocket_paragraph(self, pocket_list: list[tuple[int, str, bool]]) -> str:
        """generate paragraph of the protein pocket based on the InputConfig data.

        The envision output is a protein sequence spaced after 10 residues. Pocket residues
        should be highlighted with red color. Residue number should be tagged at the end of the
        10-residue line superscript. Residue number should reflect the actual residue number
        instead of the index starting from 0. New line should be made for every 40 printed residues.

        Example:
        RFBGDKQESD^ 10 RFBGDKQESD^ 20 RFBGDKQESD^ 30 RFBGDKQESD^ 40
        RFBGDKQESD^ 50 RFBGDKQESD^ 60 RFBGDKQESD^ 70 RFBGDKQESD^ 80
        RFBGDKQESD^ 90 RFBGDKQESD^100 RFBGDKQESD^110 RFBGDKQESD^120 <- space does not shift after 100

        Args:
            pocket_list (list[tuple[int, str, bool]]): a list of pocket information:
                                        (protein index, protein name, is_pocket (True/False))

        Returns:
            str: reportslab parseable paragraph of HTML tags
        """

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
                    text += f"</font><font size={self.residue_num_font_size}><super>{res_id_text}</super>" + \
                            f"</font><font size={self.residue_font_size}>"

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


    def generate_ligand_paragraph(self) -> Table:
        """generate ligand information. There is a table within a table to make a nice-looking
        table.

        Returns:
            Table: reportslab table detailing the ligands physicochemical information
        """
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
                           f'{float(mol.GetProp(mol_prop)):.3f}'] for mol_prop in MOL_PROPERTY_CALCULATORS] \
                            + [['SMILES', smi_text]])

        mol_data.setStyle(TableStyle([('VALIGN', (0, 0), (-1, -1), 'TOP'),
                                   ('FONTNAME', (0, 0), (-1, -1), 'Courier')]))

        table = Table([[Image(buffer, 3*inch, 3*inch, hAlign='LEFT'), mol_data]],
                      rowHeights=[3*inch])
        table.setStyle(TableStyle([('VALIGN', (0, 0), (-1, -1), 'MIDDLE')]))

        return table


    def render(self) -> list[Any]:
        """rendering the input report as a story

        Returns:
            list[Any]: story consisting of reportslab Paragraph/spacer/tables
        """

        story = []

        story.append(Paragraph(f"<font size= {self.font_size}>Protein Information</font>",
                               self.styles['Heading2']))
        story.append(self.title_spacer)

        pocket_txt = self.generate_pocket_paragraph(self.input_data.fetch_residue_is_pocket())
        story.append(Paragraph(pocket_txt, self.mono_style))
        story.append(self.section_spacer)

        story.append(Paragraph(f"<font size= {self.residue_caption_font_size}>" +
                                "Pocket Residues within 10 A of ligands are colored in red</font>",
                               self.styles['Normal']))
        story.append(self.section_spacer)

        story.append(Paragraph(f"<font size= {self.font_size}>Ligand Information</font>",
                               self.styles['Heading2']))

        story.append(self.generate_ligand_paragraph())

        return story


    def cleanup(self):
        """cleanup the image buffer, because you need to render it completely first before removing it
        """

        for img in self.image_buffer_list:
            img.close()


class ScrambleBenchSummaryJSONReport:
    """class containing the json section of the report, detailing diversity, model performance, and virtual hit.
    """

    def __init__(self, input_data: config_input.InputConfig):
        """initializing class

        Args:
            input_data (config_input.InputConfig): InputConfig
        """
        self.json_data = input_data
        self.title_font_size = 18
        self.post_processing = False
        self.diversity = False
        self.genbench3d = False
        self.virtual_hit = False

        self.styles = getSampleStyleSheet()
        self.styles['Heading1'].alignment = TA_CENTER
        self.styles['Normal'].leading = 16

        self.mini_heading_styles = getSampleStyleSheet()['Normal']
        self.mini_heading_styles.leading = 16
        self.mini_heading_styles.alignment = TA_CENTER
        self.mini_heading_styles.fontName='Helvetica-Bold'

        self.paragraph_table_style = ParagraphStyle(name='CenterText', alignment=TA_CENTER)

        self.table_style = TableStyle([
            # Header styling
            ('BACKGROUND', (0, 0), (-1, 0), colors.HexColor('#2C3E50')),
            ('TEXTCOLOR', (0, 0), (-1, 0), colors.white),
            ('FONTNAME', (0, 0), (-1, 0), 'Helvetica-Bold'),
            ('FONTSIZE', (0, 0), (-1, 0), 12),
            ('BOTTOMPADDING', (0, 0), (-1, 0), 12),

            # Body styling
            ('FONTNAME', (0, 1), (-1, -1), 'Helvetica'),
            ('FONTSIZE', (0, 1), (-1, -1), 10),
            ('ALIGN', (1, 0), (-1, -1), 'CENTER'),
            ('VALIGN', (0, 0), (-1, -1), 'MIDDLE'),

            # Alternating row colors
            ('ROWBACKGROUNDS', (0, 1), (-1, -1), [colors.HexColor('#F8F9F9'), colors.white]),

            # Borders
            #('LINEBELOW', (0, 0), (-1, 0), 2, colors.HexColor('#2C3E50')),
            ('GRID', (0, 1), (-1, -1), 0.5, colors.HexColor('#BDC3C7')),
            ('TOPPADDING', (0, 1), (-1, -1), 10),
            ('BOTTOMPADDING', (0, 1), (-1, -1), 10),
        ])

    def detect_analysis_type(self):
        """detect if types of analysis is done based on the output key in the json dictionary
        """
        for value in self.json_data.values():
            scramblebench_json_keys = value.keys()
            if 'total_mol' in scramblebench_json_keys and 'uniqueness' in scramblebench_json_keys:
                self.post_processing = True
            if 'unminimised_Validity3D' in scramblebench_json_keys or 'minimised_Validity3D' in scramblebench_json_keys:
                self.genbench3d = True
            if 'diversity' in scramblebench_json_keys:
                self.diversity = True
            if 'virtual_hit' in scramblebench_json_keys:
                self.virtual_hit = True


    def bold_table(self, df: pd.DataFrame, index: list[str] | str,
                   filter_method: str='max') -> pd.DataFrame:
        """bold the cells in the table with the best (min or max) value. Although this is a weird
        way to bold the table, I don't know how to bold using a Table directly. This is much simpler
        in my opinion. Although you cannot use Table's styling for the bolded Paragraph (because 
        Paragraph is not text, you know).

        Args:
            df (pd.DataFrame): dataframe
            index (list[str] | str): the index names of the dataframe
            filter_method (str, optional): min or max. Defaults to 'max'.

        Raises:
            ValueError: if filter is not min or max

        Returns:
            pd.DataFrame: dataframe with the best value replaced with reportslab Paragraph class
        """
        if not isinstance(index, list):
            index = [index]

        for metric in index:
            if filter_method == 'max':
                value = df.loc[metric].max()
            elif filter_method == 'min':
                value = df.loc[metric].min()
            else:
                raise ValueError(f'{value} to bold table is not supported')

            df.loc[metric, df.loc[metric] == value] = Paragraph(f"<b>{value}</b>", self.paragraph_table_style)

        return df


    def generate_summary_virtual_hit(self) -> tuple[pd.DataFrame, dict[str, Any]]:
        """generate information for virtual hit.

        Returns:
            tuple[pd.DataFrame, dict[str, Any]]: dataframe containing the hit rate of each model
            and a dictionary detailing the query, filter, and the molecule id of the virtual hit.

            I separated the molecule id from the dataframe because in the case where there are hundreds
            of virtual hit, this will cause error in reportslab where they cannot fit in 1 table. I also
            can't color the table without separating it.

        """

        virtual_hit_df = pd.DataFrame()
        virtual_hit_model_dict = defaultdict(dict)

        for model, data in self.json_data.items():

            virtual_hit_dict = data['virtual_hit']

            virtual_hit_model_dict['name'][model] = ',\n'.join(virtual_hit_dict['name'])
            virtual_hit_model_dict['query'] = virtual_hit_dict['query']
            virtual_hit_model_dict['filter'] = virtual_hit_dict['filter']

            virtual_hit_dict = {'rate' : virtual_hit_dict['rate'],
                                'model': model}

            df = pd.DataFrame(virtual_hit_dict, index=[0])

            if virtual_hit_df.empty:
                virtual_hit_df = df
            else:
                virtual_hit_df = pd.concat([virtual_hit_df, df])

        virtual_hit_df = virtual_hit_df.set_index('model').T.astype(object)
        virtual_hit_df = self.bold_table(df=virtual_hit_df,
                                index='rate',
                                filter_method='max')
        virtual_hit_df = virtual_hit_df.reset_index(names='Metric')

        return virtual_hit_df, virtual_hit_model_dict


    def generate_summary_diversity(self) -> list[tuple[str, pd.DataFrame]]:
        """generate information for library diversity

        Returns:
            list[tuple[str, pd.DataFrame]]: list of (diversity method and dataframe).
            The structure is made this way, because I want to have a method as index, and
            model name as column name to make it consistent with previous tables. However,
            because the default structure for dataframe is not like that, I need to make a
            roundabout way to structure it like this.
        """

        diversity_df = pd.DataFrame()
        diversity_df_list = []

        for model, data in self.json_data.items():
            df = pd.DataFrame(data['diversity'])
            df['model'] = model

            if diversity_df.empty:
                diversity_df = df
            else:
                diversity_df = pd.concat([diversity_df, df])

        for method, data in diversity_df.groupby('method'):
            df = data.drop(columns=['method']).set_index('model').T.round(2).astype(object)
            df = self.bold_table(df=df,
                            index='score',
                            filter_method='max')

            df = df.reset_index(names='Metric')
            diversity_df_list.append((method, df))

        return diversity_df_list


    def generate_summary_general(self) -> pd.DataFrame:
        """generate information that can be expressed as a single digit value

        Returns:
            pd.DataFrame: dataframe that details anything but virtual hit and diversity
            (usually model performance and validity3d)
        """

        single_digit_dict = defaultdict(dict)


        for model, data in self.json_data.items():
            general_data = {}
            for data_key in data:
                if data_key in ['virtual_hit', 'diversity']:
                    continue
                general_data[data_key] = data[data_key]

            single_digit_dict[model] = general_data

        single_digit_df = pd.DataFrame(single_digit_dict).round(2).astype(object)

        return self.bold_table(df=single_digit_df,
                               index=single_digit_df.index.tolist(),
                               filter_method='max').reset_index(names='Metric')


    def render(self) -> list[Any]:
        """render the json section of the report

        Returns:
            list[Any]: story as a list of paragraph/spacer/pagebreak/tables
        """

        story = []
        story.append(PageBreak())

        story.append(Paragraph(f"<font size={self.title_font_size}>ScrambleBench Analysis</font>",
                               self.styles['Heading1']))

        story.append(Spacer(1,6))
        self.detect_analysis_type()

        parameter_txt = f"<b>Ran Analysis</b> <br/>\
                        PostProcessing {add_space(16)}: {self.post_processing}<br/>\
                        GenBench3D  {add_space(19)}: {self.genbench3d} <br/>\
                        Diversity {add_space(27)}: {self.diversity}<br/>"

        story.append(Paragraph(f"<font size= 12>{parameter_txt}</font>", self.styles['Normal']))
        story.append(Spacer(1,20))

        general_df = self.generate_summary_general()

        story.append(Paragraph("<font size= 12>General Analysis Result </font>",
                               self.mini_heading_styles))
        story.append(Spacer(1,6))
        non_diversity_table = Table([general_df.columns[:,].values.astype(str).tolist()] + \
                                    general_df.values.tolist(), hAlign='LEFT')
        non_diversity_table.setStyle(self.table_style)
        story.append(non_diversity_table)

        story.append(Spacer(1,10))

        if self.virtual_hit:
            virtual_hit_df, virtual_hit_metadata = self.generate_summary_virtual_hit()

            story.append(Paragraph("<font size= 12>Virtual Hit Analysis Result </font>",
                                   self.mini_heading_styles))
            story.append(Paragraph(f"<font size= 10>Filter: {virtual_hit_metadata['filter']}, " + \
                                   f"Query: {virtual_hit_metadata['query']} </font>",
                                   self.styles['Normal']))

            virtual_hit_table = Table([virtual_hit_df.columns[:,].values.astype(str).tolist()] + \
                                       virtual_hit_df.values.tolist(), hAlign='LEFT', splitInRow=1)
            virtual_hit_table.setStyle(self.table_style)
            story.append(virtual_hit_table)
            story.append(Spacer(1,10))
            for model, hits in virtual_hit_metadata['name'].items():
                if hits:
                    story.append(Paragraph(f'{model} hits: {hits}\n'))
                    story.append(Spacer(1,10))
            story.append(Spacer(1,6))

        if self.diversity:

            story.append(Paragraph("<font size= 12>Diversity Analysis Result </font>",
                                   self.mini_heading_styles))
            story.append(Spacer(1,6))
            for diversity_tuple in self.generate_summary_diversity():
                story.append(Paragraph(f"<font size= 10><b>{diversity_tuple[0]}</b></font>",
                                       self.styles['Normal']))
                diversity_df = diversity_tuple[-1]
                diversity_table = Table([diversity_df.columns[:,].values.astype(str).tolist()] + \
                                        diversity_df.values.tolist(), hAlign='LEFT')
                diversity_table.setStyle(self.table_style)
                story.append(diversity_table)
                story.append(Spacer(1,6))

        # story.append(Paragraph(f"<font size= 12>{batch_parameter_txt}</font>", self.styles['Normal']))
        # story.append(Spacer(1,12))

        return story


class ScrambleBenchPlotReport:
    """class containing the plotting section of the report
    """

    def __init__(self, summary_df: pd.DataFrame,
                 columns: list[str], plot_type: str='violin', title_caption: str=None):
        """initialize class

        Args:
            summary_df (pd.DataFrame): dataframe from p5_collect_data_analysis.py
            columns (list[str]): columns to be plotted
            plot_type (str, optional): type of plot. Defaults to 'violin'.
            title_caption (str, optional): title caption used in matplotlib. Defaults to None.
        """

        self.summary_df = summary_df
        self.columns = columns
        self.title_caption = title_caption
        self.plot_type = plot_type

        self.image_buffer_list = []
        self.title_font_size = 18
        self.styles = getSampleStyleSheet()
        self.styles['Heading1'].alignment = TA_CENTER
        self.styles['Normal'].leading = 16


    def generate_image_render(self, dataframe: pd.DataFrame,
                              columns: list[str], hue: str='Model', x: str='protein_name') -> io.BytesIO:
        """generate the plot as an image that can be pasted to reportslab

        Args:
            dataframe (pd.DataFrame): the filtered dataframe
            columns (list[str]): columns to be plotted in y axis
            hue (str, optional): hue in seaborn. Defaults to 'Model'.
            x (str, optional): columns to be plotted in x axis. Defaults to 'protein_name'.

        Returns:
            io.BytesIO: _description_
        """

        if len(columns) < 3:
            fig, axs = plt.subplots(len(columns), 1, figsize=(7,7), sharex=True)

        else:
            fig, axs = plt.subplots(math.ceil(len(columns) / 2), 2, figsize=(12, 12), sharex=True)

        if self.plot_type == 'violin':
            plot_function = self.plot_violin
        elif self.plot_type == 'box':
            plot_function = self.plot_box
        elif self.plot_type == 'raincloud':
            plot_function = self.plot_raincloud
        else:
            logging.warning(f'Plot type {self.plot_type} is not recognised. Defaulting to violin plot.')
            plot_function = self.plot_violin

        for column, ax in zip(columns, np.ravel(axs)[:len(columns)]):
            plot_function(dataframe=dataframe, column=column, ax=ax, hue=hue, x=x)

            handles, labels = ax.get_legend_handles_labels()
            ax.get_legend().remove()

            ax.set_xlabel('Protein Name', fontsize=12)
            ax.set_ylabel(' '.join([word.capitalize() if word.upper() != word else word for word in column.split('_') ]),
                          fontsize=12)
            ax.tick_params(axis='both', which='major', labelsize=12)


        legend_cols = 3
        top_margin = 0.95 - (0.05 * (math.ceil(len(handles) / legend_cols) - 1))

        buf = io.BytesIO()

        fig.subplots_adjust(top=top_margin, bottom=0.1, left=0.1, right=0.9)

        fig.legend(
            handles,
            labels,
            loc="lower center",  # 'lower center' works great with bbox_to_anchor at the top
            bbox_to_anchor=(0.5, top_margin),
            ncol=legend_cols,
            prop={"size": 12},
            frameon=False,
        )

        plt.tight_layout(rect=[0, 0, 1, top_margin])
        plt.savefig(buf, format = 'png', dpi=300)
        buf.seek(0)
        plt.close()

        return buf


    def render(self) -> list[Any]:
        """render the plotting section of the report as a story

        Returns:
            list[Any]: story consisting of pagebreak/paragraph/image/spacer
        """

        story = []

        story.append(PageBreak())

        image_buffer = self.generate_image_render(dataframe=self.summary_df,
                                                  columns=self.columns,
                                                  hue='Model',
                                                  x='protein_name')

        self.image_buffer_list.append(image_buffer)

        if self.title_caption:
            story.append(Paragraph(f"<font size= {self.title_font_size}>{self.title_caption}</font>",
                                     self.styles['Heading1']))

        story.append(Image(image_buffer, 7*inch, 7*inch))

        return story


    def plot_violin(self, dataframe: pd.DataFrame,
                    column: str, ax: Axes, hue: str='Model', x: str='protein_name') -> Axes:
        """plot using violinplot

        Args:
            dataframe (pd.DataFrame): the filtered dataframe
            columns (list[str]): columns to be plotted in y axis
            ax (Axes): initialized ax from matplotlib
            hue (str, optional): hue in seaborn. Defaults to 'Model'.
            x (str, optional): columns to be plotted in x axis. Defaults to 'protein_name'.

        Returns:
            Axes: ax with the violin plot
        """


        return sns.violinplot(data=dataframe,
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


    def plot_raincloud(self, dataframe: pd.DataFrame,
                    column: str, ax: Axes, hue: str='Model', x: str='protein_name') -> Axes:
        """plot using raincloud

        Args:
            dataframe (pd.DataFrame): the filtered dataframe
            columns (list[str]): columns to be plotted in y axis
            ax (Axes): initialized ax from matplotlib
            hue (str, optional): hue in seaborn. Defaults to 'Model'.
            x (str, optional): columns to be plotted in x axis. Defaults to 'protein_name'.

        Returns:
            Axes: ax with the raincloud plot
        """

        import ptitprince as pt
        return pt.RainCloud(hue=hue,
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


    def plot_box(self, dataframe: pd.DataFrame,
                    column: str, ax: Axes, hue: str='Model', x: str='protein_name') -> Axes:
        """plot using boxplot

        Args:
            dataframe (pd.DataFrame): the filtered dataframe
            columns (list[str]): columns to be plotted in y axis
            ax (Axes): initialized ax from matplotlib
            hue (str, optional): hue in seaborn. Defaults to 'Model'.
            x (str, optional): columns to be plotted in x axis. Defaults to 'protein_name'.

        Returns:
            Axes: ax with the box plot
        """
        return sns.boxplot(data=dataframe,
                                x=x,
                                y=column,
                                ax=ax,
                                palette=COLORBLIND_PALETTE,
                                hue=hue,
                                width=0.5,
                                orient= 'v',
                                dodge= True,
                                linewidth=2)


    def cleanup(self):
        """cleanup the image buffer, because you need to render it completely first before removing it
        """
        for img in self.image_buffer_list:
            img.close()

ReportClass = ScrambleBenchInputReport | ScrambleBenchParameterReport | \
              ScrambleBenchPlotReport | ScrambleBenchSummaryJSONReport

class ScrambleBenchReport:
    """class containing the overall report interface
    """

    def __init__(self, yaml_filename: str):
        """initialize class

        Args:
            yaml_filename (str): user's prepared YAML filename
        """
        self.yaml_file = yaml_filename
        self.section = []

    def add_section(self, section: ReportClass):
        """add the report section to the interface

        Args:
            section (ReportClass): ScrambleBenchInputReport | ScrambleBenchParameterReport | 
                                   ScrambleBenchPlotReport | ScrambleBenchSummaryJSONReport
        """
        self.section.append(section)

    def generate_report(self):
        """generate the report after adding the section together
        """

        doc = SimpleDocTemplate(str(Path(self.yaml_file).parent / 'report.pdf'), pagesize=letter,
                                rightMargin=inch/2, leftMargin=inch/2,
                                topMargin=40, bottomMargin=18)


        story = []

        try:
            for section in self.section:
                story += section.render()

            doc.build(story)

        finally:

            for section in self.section:
                cleanup_func = getattr(section, "cleanup", None)
                if callable(cleanup_func):
                    cleanup_func()


def generate_report(yaml_filename: str):
    """main function of p6_generate_report.py. This will generate a PDF report depending on the 
    parameter given by the user in the YAML filename. The input requires the output from
    p5_collect_data_analysis.py (i.e., summary.csv and summary.json)

    Args:
        yaml_filename (str): user's prepared YAML filename
    """

    input_csv = Path(yaml_filename).parent / 'summary.csv'
    input_json = Path(yaml_filename).parent / 'summary.json'

    json_data = None
    with open(input_json, 'r', encoding='utf-8') as json_f:
        json_data = json.load(json_f)

    df = pd.read_csv(input_csv)

    with open(yaml_filename, 'r', encoding='utf-8') as yaml_f:
        config_data = yaml.safe_load(yaml_f)

    report_data = config_report.ReportConfig(config_data)
    parameter_data = config_parameter.ParameterConfig(config_data)

    report = ScrambleBenchReport(yaml_filename)
    report.add_section(ScrambleBenchParameterReport(parameter_data=parameter_data))


    if config_constant.INPUT_KEY in config_data:
        input_data = config_input.InputConfig(config_data)
        report.add_section(ScrambleBenchInputReport(input_data=input_data))

    if json_data:
        report.add_section(ScrambleBenchSummaryJSONReport(input_data=json_data))

    if report_data.docking_score_value:
        redocking_score_columns = [col for col in df.columns if 'redocking_score' in col.lower()]
        report.add_section(ScrambleBenchPlotReport(summary_df=df,
                                                columns=redocking_score_columns,
                                                plot_type=report_data.plot_value,
                                                title_caption='Redocking Score'))

    if report_data.rmsd_value:
        rmsd_columns = [col for col in df.columns if 'rmsd' in col.lower()]
        report.add_section(ScrambleBenchPlotReport(summary_df=df,
                                                columns=rmsd_columns,
                                                plot_type=report_data.plot_value,
                                                title_caption='RMSD Score'))

    if report_data.qed_value:
        qed_column = [col for col in df.columns if 'qed' in col.lower()]
        report.add_section(ScrambleBenchPlotReport(summary_df=df,
                                                columns=qed_column,
                                                plot_type=report_data.plot_value,
                                                title_caption='QED'))

    if report_data.validity3d_value:
        validity3d_column = [col for col in df.columns if 'validity3d' in col.lower()]
        report.add_section(ScrambleBenchPlotReport(summary_df=df,
                                                columns=validity3d_column,
                                                plot_type=report_data.plot_value,
                                                title_caption='Validity3D Metric'))

    report.generate_report()


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Generate report from p5_collect_data_analysis.py")

    parser.add_argument("-i", "--input",
                        help="config yaml input file or txt file containing yaml filepath",
                        required=True, type=str)

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
