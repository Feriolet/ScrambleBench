"""This util file handles the protein preparation and protonation for vina and Glide.
Modified from the GenBench3D repository. I am not sure why the class is so confusing, so
I won't be able to document it pretty well."""
import os
import logging
import sys
import subprocess

from pathlib import Path

import rdkit

from rdkit import Chem
from openmm.app import PDBFile
from pdbfixer import PDBFixer
from openbabel import pybel

logger = logging.getLogger(__name__)

class VinaProtein():
    """class to prepare protein for Vina docking
    """

    def __init__(self,
                 pdb_filepath: str,
                 prepare_receptor_bin_path: str,
                 preparation_method:str='adfr',
                 protonation_method:str='pdbfixer') -> None:
        """initialize class

        Args:
            pdb_filepath (str): .pdb protein file
            prepare_receptor_bin_path (str): executable directory to prepare_receptor
            preparation_method (str, optional): program to prepare protein. Defaults to 'adfr'.
            protonation_method (str, optional): program to protonate protein. Defaults to 'pdbfixer'.
        """
        self.protein_filepath = pdb_filepath
        self.prepare_receptor_bin_path = prepare_receptor_bin_path
        self._pdbqt_filepath = pdb_filepath.replace('.pdb',
                                                   '.pdbqt')
        self._protein_clean_filepath = pdb_filepath.replace('.pdb',
                                                           '_protein_only_clean.pdb')
        self.preparation_method = preparation_method
        self.protonation_method = protonation_method

    @property
    def pdbqt_filepath(self) -> str:
        """preparing the pdbqt filepath as a property.

        Returns:
            str: the final .pdbqt file for docking
        """
        if not os.path.exists(self._pdbqt_filepath):
            self.vina_prepare_receptor(output_pdbqt_filepath=self._pdbqt_filepath) # Using default configuration
        assert os.path.exists(self._pdbqt_filepath), \
            'Something went wrong during Vina receptor preparation'
        return self._pdbqt_filepath

    @property
    def protein_clean_filepath(self) -> str:
        """clean the .pdb using pdbfixer. I am not sure why this is done by the Genbench3d author.
        But I don't want to touch it.

        Returns:
            str: cleaned .pdb file before preparation
        """
        if not os.path.exists(self._protein_clean_filepath):
            self.clean_protein(input_pdb_filepath=self.protein_filepath,
                                 output_pdb_filepath=self._protein_clean_filepath)
        assert os.path.exists(self._protein_clean_filepath), \
            'Something went wrong during protein cleaning'
        return self._protein_clean_filepath

    @protein_clean_filepath.setter
    def protein_clean_filepath(self, new_filepath: str):
        """setter

        Args:
            new_filepath (str): The clean .pdb file
        """
        # Add any validation or logic here
        self._protein_clean_filepath = new_filepath

    def vina_prepare_receptor(self,
                              output_pdbqt_filepath: str,
                                pH: float = 7.4
                                ) -> None:
        """
        inspired from teachopencadd talktorial 15 on protein_ligand_docking.

        I modified this code to accept schrodinger as another protein protonation method.
        Otherwise, not much is changed from the original line of code.

        """

        if self.protonation_method == 'pdbfixer':
            self.clean_protein(input_pdb_filepath=self.protein_filepath,
                            output_pdb_filepath=self.protein_clean_filepath,
                            pH=pH)

        elif 'schrodinger' in self.protonation_method:

            current_dir = Path.cwd()
            os.chdir(Path(self.protein_filepath).parent)
            command = [str(Path(self.protonation_method) / 'utilities/prepwizard'),
                       str(Path(self.protein_filepath).name),
                       str(Path(self.protein_clean_filepath).name), '-fillsidechains',
                        '-disulfides', '-mse', '-assign_all_residues', '-fillloops', '-rehtreat', 
                        '-max_states', '1', '-epik_pH', '7.4', '-epik_pHt', '1.0', '-samplewater',
                        '-include_epik_states', '-propka_pH', '7.4', '-f', 'S-OPLS', '-rmsd', '0.3',
                        '-watdist','5.0','-JOBNAME','proteinprep_4','-HOST','localhost:7', '-WAIT']

            subprocess.run(command, check=True)

            os.chdir(current_dir)

        else:
            self.protein_clean_filepath = self.protein_filepath

        if self.preparation_method == 'adfr':
            self.adfr_receptor_preparation(input_pdb_filepath=self.protein_clean_filepath,
                                           output_pdbqt_filepath=output_pdbqt_filepath)
        else:
            self.pdb_to_pdbqt()


    def adfr_receptor_preparation(self,
                                  input_pdb_filepath: str,
                                  output_pdbqt_filepath: str,
                                  ) -> None:
        """
        input_pdb_filepath must be a pbd file that only contains the protein with
        hydrogens
        """
        logging.info(f'Preparing protein from {input_pdb_filepath} to {output_pdbqt_filepath}')
        arg_list = [self.prepare_receptor_bin_path,
                    f'-r {input_pdb_filepath}',
                    f'-o {output_pdbqt_filepath}']
        cmd = ' '.join(arg_list)
        os.system(cmd)


    def pdb_to_pdbqt(self,
                     ph: float = 7.4,
                     ) -> None:
        """convert pdb to pdbqt by pybel

        Args:
            ph (float, optional): pH environment. Defaults to 7.4.
        """
        molecule = list(pybel.readfile("pdb", str(self.protein_filepath)))[0]
        self.ob_mol_to_pdbqt(molecule, ph)


    def ob_mol_to_pdbqt(self,
                        molecule: pybel.Molecule,
                        ph: float = 7.4,
                        ) -> None:
        """molecule preparation to pdbqt by obabel 

        Args:
            molecule (pybel.Molecule): ligands
            ph (float, optional): pH environment. Defaults to 7.4.
        """

        # add hydrogens at given pH
        molecule.OBMol.CorrectForPH(ph)
        molecule.addh()
        # add partial charges to each atom
        for atom in molecule.atoms:
            atom.OBAtom.GetPartialCharge()

        molecule.write("pdbqt", str(self._pdbqt_filepath), overwrite=True)

        # Only keep ATOM and TER lines in pdbqt file
        with open(self.pdbqt_filepath, 'r', encoding='utf-8') as f:
            lines = [line.strip() for line in f.readlines()]
        new_lines = [line
                     for line in lines
                     if line.startswith('ATOM') or line.startswith('TER')]
        with open(self.pdbqt_filepath, 'w', encoding='utf-8') as f:
            for line in new_lines:
                f.write(line)
                f.write('\n')


    def clean_protein(self,
                      input_pdb_filepath: str,
                      output_pdb_filepath: str,
                      pH: float = 7.4,
                      ) -> None:
        """
        Adds hydrogens to given pH for a pdb protein (input_filepath)
        """
        logging.info(f'Cleaning protein from {input_pdb_filepath} to {output_pdb_filepath}')

        fixer = PDBFixer(filename=input_pdb_filepath)
        # We only use PDBFixer to reinitialize the chains and segment id combinations
        # (e.g. chain A and B, having both segment A and B, will be renamed chain A to D)
        fixer.findMissingResidues() # it cannot find them with CrossDocked because it does not contain sequence info
        fixer.findNonstandardResidues()
        fixer.replaceNonstandardResidues()
        fixer.findMissingAtoms()
        fixer.addMissingAtoms()
        fixer.addMissingHydrogens(pH=pH)
        fixer.removeHeterogens(keepWater=False)
        intermediate_filepath = output_pdb_filepath.replace('.pdb', '_no_h.pdb')

        with open(intermediate_filepath, 'w', encoding='utf-8') as intermediate_f:
            PDBFile.writeFile(fixer.topology, fixer.positions, intermediate_f)

        # pdbfixer might miss some hydrogens because of missing atoms
        molecule = list(pybel.readfile("pdb", str(intermediate_filepath)))[0]
        # import pdb;pdb.set_trace()
        molecule.OBMol.CorrectForPH(pH) # the correct functions fails for some groups: -NH3+ becomes -NH4+
        molecule.removeh() # some pdb files have hydrogens, these might mess up the next step
        molecule.addh() 
        molecule.write("pdb", str(output_pdb_filepath), overwrite=True)


class GlideProtein():
    """class containing protein preparation and protonation for Glide SP docking
    """
    def __init__(self,
                 pdb_filepath: str,
                 native_ligand: Chem.Mol,
                 grid_output_dirpath: str,
                 schrodinger_dirpath: str,
                 protein_preparation: bool) -> None:
        """initialize class

        Args:
            pdb_filepath (str): .pdb protein file
            native_ligand (Chem.Mol): native ligand as RDKit molecule
            grid_output_dirpath (str): grid output directory path
            schrodinger_dirpath (str): root directory for schrodinger
            protein_preparation (bool): whether to prepare the protein using prot wizard

        Raises:
            FileNotFoundError: if Glide executable cannot be found.
        """

        glide_path = str(Path(schrodinger_dirpath) / 'glide')
        structconvert_path = str(Path(schrodinger_dirpath) / 'utilities/structconvert')
        protwizard_path = str(Path(schrodinger_dirpath) / 'utilities/prepwizard')
        if not os.path.exists(glide_path):
            raise FileNotFoundError(f'Cannot find Glide executable at {glide_path}')

        self.pdb_filepath = pdb_filepath
        self.grid_output_dirpath = grid_output_dirpath
        if not os.path.exists(self.grid_output_dirpath):
            os.mkdir(self.grid_output_dirpath)

        self.glide_path = glide_path
        self.structconvert_path = structconvert_path
        self.protwizard_path = protwizard_path
        
        self.mae_filepath = pdb_filepath.replace('.pdb',
                                                '.mae')
        if not os.path.exists(self.mae_filepath):
            self.generate_mae_file()

        assert os.path.exists(self.mae_filepath), 'Schrodinger structconvert failed to generate .mae file.'

        if protein_preparation:
            self.prepared_mae_filepath = str(Path(self.mae_filepath).parent / f'{Path(self.mae_filepath).stem}_prepared.mae')
            if not Path(self.prepared_mae_filepath).is_file():
                self.prepare_protein()
        else:
            self.prepared_mae_filepath = self.mae_filepath

        self.grid_filepath = self.prepared_mae_filepath.replace('.mae',
                                                                '.zip')

        self.grid_center = native_ligand.GetConformer().GetPositions().mean(axis=0)

        self.glide_grid_in_filename = str(Path(self.prepared_mae_filepath).parent / 'glide_grid_generation.in')
        self.glide_grid_in_filepath = self.glide_grid_in_filename

        if not os.path.exists(self.grid_filepath):
            if os.path.exists(self.glide_grid_in_filepath):
                os.remove(self.glide_grid_in_filepath)
            self.generate_glide_grid_in_file(self.grid_center)
            assert os.path.exists(self.glide_grid_in_filepath)
            self.generate_grid_file()

        assert os.path.exists(self.grid_filepath), 'Schrodinger Glide failed to generate grid file.'


    def generate_mae_file(self):
        """generate .pdb to .mae file for compatability
        """
        logging.info(f'Converting {self.pdb_filepath} to {self.mae_filepath}')
        command = [self.structconvert_path,
                   self.pdb_filepath,
                   self.mae_filepath]
        subprocess.run(command, check=True)


    def generate_glide_grid_in_file(self,
                                    grid_center: list[float]):
        """generate grid .in for Glide. format is compatible as of 2025-4 version

        Args:
            grid_center (list[float]): pocket coordinate
        """
        # List of keywords available in the Glide documentation


        logging.info(f'Writing glide grid generation input in {self.glide_grid_in_filepath}')
        grid_center_str = [str(value) for value in grid_center]
        d = {'FORCEFIELD': 'OPLS_2005',
            'GRIDFILE': str(Path(self.grid_filepath).name),
             'RECEP_FILE': str(Path(self.prepared_mae_filepath).name),
             'GRID_CENTER': ','.join(grid_center_str)}
        with open(self.glide_grid_in_filepath, 'w', encoding='utf-8') as f:
            for param_name, value in d.items():
                f.write(f'{param_name}   {value}')
                f.write('\n')

    def generate_grid_file(self):
        """generate grid file for Glide SP. Command line is compatible as of 2025-4 version.
        """
        current_dir = Path.cwd()
        os.chdir(Path(self.prepared_mae_filepath).parent)

        logging.info(f'Generate Glide grid using {self.glide_grid_in_filepath}')
        command = [f'{self.glide_path}',
                   self.glide_grid_in_filepath,
                   '-OVERWRITE', '-HOST', 'localhost', '-TMPLAUNCHDIR', '-WAIT']
        subprocess.run(command, check=True)

        os.chdir(current_dir)

    def prepare_protein(self):
        """prepare protein using ProtWizard. Command line is compatible as of 2025-4 version.
        """
        logging.info(f'Preparing protein using ProtWizard {self.mae_filepath}')
        current_dir = Path.cwd()
        os.chdir(Path(self.prepared_mae_filepath).parent)
        prepared_mae_filename = Path(self.prepared_mae_filepath).name

        command = [self.protwizard_path, str(Path(self.mae_filepath).name), prepared_mae_filename, '-fillsidechains',
                   '-disulfides', '-mse', '-assign_all_residues', '-fillloops', '-rehtreat', 
                   '-max_states', '1', '-epik_pH', '7.4', '-epik_pHt', '1.0', '-samplewater',
                   '-include_epik_states', '-propka_pH', '7.4', '-f', 'S-OPLS', '-rmsd', '0.3',
                   '-watdist','5.0','-JOBNAME','proteinprep_4','-HOST','localhost:7', '-WAIT']

        subprocess.run(command, check=True)

        subprocess.run(['mv', prepared_mae_filename, self.prepared_mae_filepath], check=True)

        os.chdir(current_dir)

if __name__ == '__main__':

    logging.basicConfig(stream=sys.stdout,
                    level=logging.INFO,
                    format='%(asctime)s - %(module)s: - %(levelname)s - %(message)s',
                    datefmt='%Y-%m-%d %H:%M:%S')
