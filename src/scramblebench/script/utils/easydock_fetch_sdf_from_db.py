import sqlite3
import rdkit
from rdkit import Chem
import argparse

def get_sdf_from_db(db_fname, output_fname):
    conn = sqlite3.connect(db_fname, timeout=60)
    cur = conn.cursor()
    data = cur.execute(f"""SELECT source_mol_block_protonated FROM mols""").fetchall()

    with Chem.SDWriter(output_fname) as output_f:
        for datum in data:
            mol_block = datum[0]
            if mol_block:
                mol = Chem.MolFromMolBlock(datum[0])
                output_f.write(mol)

if __name__ == '__main__':
    parser = argparse.ArgumentParser(description="Prepare generated molecules for downstream analysis")

    parser.add_argument("-i", "--input", help="db input of easydock output", required=True, type=str)
    parser.add_argument('-o', '--output', help='sdf output of protonated mol', required=True)
    args = parser.parse_args()

    get_sdf_from_db(args.input, args.output)