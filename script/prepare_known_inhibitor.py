
from rdkit import Chem
import os

def prepare_known_inhibitor(input_fname, protein_name, output_fname):
    if input_fname[-4:] == '.smi':
        with open(input_fname, 'r') as known_inhibitor_fname:
            known_inhibitor_data = known_inhibitor_fname.readlines()
        
        known_inhibitor_smi_l = [known_inhibitor.split()[0] for known_inhibitor in known_inhibitor_data]
        known_inhibitor_name_l = [known_inhibitor.split('\t',1)[-1].strip() for known_inhibitor in known_inhibitor_data]

        known_inhibitor_mol_l = [Chem.MolFromSmiles(smi) for smi in known_inhibitor_smi_l]
        known_inhibitor_mol_l = [mol for mol in known_inhibitor_mol_l]
        
        for index, (known_inhibitor_name, known_inhibitor_mol) in enumerate(zip(known_inhibitor_name_l, known_inhibitor_mol_l)):
            if known_inhibitor_mol and known_inhibitor_mol != '':
                print(known_inhibitor_mol)
                known_inhibitor_mol.SetProp('original_name', known_inhibitor_name)
                known_inhibitor_mol.SetProp('_Name', f'control_{index}')
                known_inhibitor_mol.SetProp('Protein', protein_name)
                known_inhibitor_mol.SetProp('GenAI_Model', 'Control')
                print(known_inhibitor_mol.GetProp('original_name'))
    elif input_fname[-4:] == '.sdf':
        known_inhibitor_mol_l = [mol for mol in Chem.SDMolSupplier(input_fname)]
        for index, known_inhibitor_mol in enumerate(known_inhibitor_mol_l):
            if known_inhibitor_mol and known_inhibitor_mol != '':
                try:
                    known_inhibitor_mol.SetProp('original_name', known_inhibitor_mol.GetProp('_Name'))
                except KeyError:
                    pass
                known_inhibitor_mol.SetProp('_Name', f'control_{index}')
                known_inhibitor_mol.SetProp('Protein', protein_name)
                known_inhibitor_mol.SetProp('GenAI_Model', 'Control')

    known_inhibitor_mol_l = [mol for mol in known_inhibitor_mol_l if mol and mol != '']
    with Chem.SDWriter(output_fname) as w:
        for m in known_inhibitor_mol_l:
            w.write(m)


if __name__ == '__main__':
    input_fname_l = ['known_inhibitor/3cl_ligand.smi']
    protein_name_l = ['Hydrolase_3CLpro']
    for input_fname, protein_name in zip(input_fname_l, protein_name_l):
        output_fname = os.path.join(input_fname.rsplit('/')[0], f'control_ligand_{protein_name}.sdf')
        prepare_known_inhibitor(input_fname, protein_name, output_fname)


