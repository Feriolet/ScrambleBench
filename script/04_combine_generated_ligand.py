from glob import glob
from rdkit import Chem
from typing import Optional, Any
import os
from glob import glob
import yaml

def add_properties_to_mol(mol: Chem.Mol,
                        property: list[list[Any]]) -> None:
    """Add list of properties to a Mol object

    Args:
        mol (Chem.Mol): a molecule object
        property (list[list[Any]]): list of properties with len of two
    """

    for prop in property:
        assert len(prop) == 2
        mol.SetProp(prop[0], str(prop[1]))

def combine_all_generated_ligand(input_dir: str, 
                                 genai_models: list[str], 
                                 additional_property: Optional[list[list[Any]]]=None) -> list[Chem.Mol]:
    """Combine generated ligand from different model contained in the same folder into one list. 
    The ligands will be identified through property in the Chem.Mol object using SetProp() method. 

    Args:
        input_dir (str): The directory path containing multiple or one sdf file
        genai_models (list[str]): List of GenAI model. The str should match the one in the file.
        additional_property (Optional[list[list[Any]]], optional): _description_. Defaults to None.

    Returns:
        list[Chem.Mol]: list of generated ligands from multiple models
    """

    cheminformatic_input_fnames = glob(f'{input_dir}/*.sdf')
    all_generated_ligand_list = []

    for model in genai_models:

        detected_sdf_file = 'Unknown'
        # Find which file belongs to the corresponding model. If no output, skipped
        for generated_input_fname in cheminformatic_input_fnames:
            if model in generated_input_fname:
                detected_sdf_file = generated_input_fname
                break
        if detected_sdf_file == 'Unknown':
            continue

        sdmol = Chem.SDMolSupplier(detected_sdf_file, removeHs=False)
        mol_list = []

        for mol in sdmol:
            if mol:
                if mol.GetNumAtoms() > 0:
                    if additional_property:
                        add_properties_to_mol(mol=mol, 
                                              property=additional_property)
                    mol_list.append(mol)
       
        all_generated_ligand_list += mol_list

    return all_generated_ligand_list


if __name__ == '__main__':
    input_dirname_dirname_dirname = '/opt/veincent/GenAI_manuscript/output_test_multiple_numsample'
    input_dirname = 'cheminformatics_input_prepared'
    output_dirname = 'merged_generated_ligand'
    genai_models = ['Pocket2Mol', 'PocketFlow', 'Lingo3DMol', 'DiffSBDD', 'PMDM', 'Chem42']
    protein_list = ['GPCR_5HT2C','GPCR_DRD2','Kinase_CDK2', 'Kinase_GSK3b','Protease_3CLpro','Protease_NA']
    

    if not os.path.isdir(output_dirname):
        os.makedirs(output_dirname, exist_ok=True)

    for protein in protein_list:
        protein_specific_mol_list = []
        protein_specific_sdf_output_fname = os.path.join(input_dirname_dirname_dirname, output_dirname, f'generated_ligand_{protein}.sdf')
        step03_dirname_list = glob(f'{input_dirname_dirname_dirname}/{protein}*/{input_dirname}')
        for step03_dirname in step03_dirname_list:
            mol_list = combine_all_generated_ligand(input_dir=step03_dirname, 
                                                    genai_models=genai_models)
            protein_specific_mol_list += mol_list

        with Chem.SDWriter(protein_specific_sdf_output_fname) as w:
            for index, m in enumerate(protein_specific_mol_list):
                w.write(m)