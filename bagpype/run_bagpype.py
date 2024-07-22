import os
import shutil
import glob
import logging
import importlib
import pathlib
import subprocess
import pandas as pd
from tqdm import tqdm
from rdkit import Chem

import create_ligand_Bagpype as create_files
from bagpype.helper_funcs import create_dir, check_dirpath


# Setup logger
logger = logging.getLogger(__name__)
# logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)s:%(asctime)s:%(name)s:%(message)s')

# Define handlers
# fh = logging.FileHandler('{}.log'.format(__name__))
fh = logging.FileHandler('batch_run_bagpype.log')
fh.setLevel(logging.INFO)
fh.setFormatter(formatter)
sh = logging.StreamHandler()
sh.setLevel(logging.WARNING)

# Add handlers to logger
logger.addHandler(fh)
logger.addHandler(sh)


class RunBagPype:
    def __init__(
            self,
            output_dir: str = 'output/',
            del_output: bool = True,
    ) -> None:

        self.output_dir: str = create_dir(check_dirpath(output_dir), replace=del_output)
        self.mymol_dict: dict = {}


    def init_from_smiles_string(self, smiles: str, mol_id: str = '1') -> None:
        """Converts input SMILES strings to MOL and MOL2 files required for BagPype. Isomeric SMILES and
        AddHs options set to True. SMI file written to disk before being converted to MOL and MOL2 files
        with the same STEM (mol_id).

        Args:
            smiles (str): SMILES string
            mol_id (str, optional): Unique ID to used as dict key. Defaults to '1'.
        """
        smi: str = Chem.MolToSmiles(Chem.AddHs(Chem.MolFromSmiles(smiles)), isomericSmiles=True)

        self.mol_files: dict[str, str] = {mol_id: self.output_dir + mol_id + '.mol'}
        self.mol2_files: dict[str, str] = {mol_id: self.output_dir + mol_id + '.mol2'}

        smi_file: str = self.output_dir + mol_id + '.smi'

        with open(smi_file, 'w') as f:
            f.writelines(smi)

        self.smiles_to_mol_mol2(
            smi_file=smi_file,
            mol_file=self.mol_files[mol_id],
            mol2_file=self.mol2_files[mol_id],
        )


    def init_from_smiles_file(self, smiles_file: str):
        pass


    def init_from_smiles_dir(self, smiles_dir: str):
        smi_dir = {get_id(file): file for file in glob.glob(check_dirpath(smiles_dir)+'*.smi')}
        self.mol_files = {get_id(file): self.output_dir + get_id(file) + '.mol' for file in smi_dir}
        self.mol2_files = {get_id(file): self.output_dir + get_id(file) + '.mol2' for file in smi_dir}

        for mol_id in smi_dir.keys():
            self.smiles_to_mol_mol2(
                smi_file=smi_dir[mol_id],
                mol_file=self.mol_files[mol_id],
                mol2_file=self.mol2_files[mol_id],
            )


    def create_mol_dir(self, mol_id):
        return create_dir(self.output_dir, f'{mol_id}/', replace=True)


    def gen_mol_mol2(self, mol_id_lst):
        self.mol_files = {mol_id: self.mol_dir[mol_id] + mol_id + '.mol' for mol_id in mol_id_lst}
        self.mol2_files = {mol_id: self.mol_dir[mol_id] + mol_id + '.mol2' for mol_id in mol_id_lst}


    def gen_smi_file(self, smi, mol_id):
        smi_file = self.mol_dir[mol_id] + mol_id + '.smi'
        with open(smi_file, 'w') as f:
            f.writelines(smi)
        return smi_file


    def init_from_csv(self, csv_file, index_col, smiles_col_name):

        smi_ser = pd.read_csv(csv_file, index_col=index_col)[smiles_col_name]
        self.mol_dir = {mol_id: self.create_mol_dir(mol_id) for mol_id in smi_ser.index}

        self.gen_mol_mol2(
            mol_id_lst=smi_ser.index,
        )

        for mol_id, smi in smi_ser.items():
            self.smiles_to_mol_mol2(
                smi_file=self.gen_smi_file(smi, mol_id),
                mol_file=self.mol_files[mol_id],
                mol2_file=self.mol2_files[mol_id],
            )


    def init_from_df(self, df):
        pass


    def smiles_to_mol_mol2(self, smi_file: str, mol_file: str, mol2_file: str) -> None:
        """Use subprocess to run obabel from the CLI to convert the SMILES file to both MOL and MOL2 files.
        Generate 3D coordiantes (--gen3d) and add Hs (-h) options used.

        Args:
            smi_file (str): SMILES file to be converted.
            mol_file (str): Filename of output MOL file.
            mol2_file (str): Filename of output MOL2 file.
        """

        command_mol = f'obabel -ismi {smi_file} -omol -O {mol_file} -h --gen3d'.split(' ')
        subprocess.run(command_mol, stdout = subprocess.PIPE)

        command_mol2 = f'obabel -ismi {smi_file} -omol2 -O {mol2_file} -h --gen3d'.split(' ')
        subprocess.run(command_mol2, stdout = subprocess.PIPE)


    def init_from_mol_file(self, mol_file, mol2_file):
        self.mol_files = {get_id(mol_file): mol_file}
        self.mol2_files = {get_id(mol2_file): mol2_file}


    def init_from_mol_dir(self, mol_dir, mol2_dir):
        self.mol_files = {get_id(file): file for file in glob.glob(check_dirpath(mol_dir)+'*.mol')}
        self.mol2_files = {get_id(file): file for file in glob.glob(check_dirpath(mol2_dir)+'*.mol2')}


    def create_bagpype_pdb(self):

        # Setup tqdm tracking bar
        pbar = tqdm(total=len(self.mol_files))

        for mol_id in self.mol_files.keys():

            pbar.set_description(f'Working on {mol_id}')

            self.run_bagpype(
                mol_id=mol_id,
                mol_file=self.mol_files[mol_id],
                mol2_file=self.mol2_files[mol_id],
            )

            pbar.update()
        
        pbar.close()


    def run_bagpype(
        self,
        mol_id,
        mol_file,
        mol2_file,
        strip={'res_name': [], 'chain': [], 'residues': [],},
        trim_H=True,
        add_H=True,
        MakeMultimer_number=None,
    ):

        import bagpype

        generate_cif_het = create_files.Ligand_generation(
            mol_file=mol_file,
            mol2_file=mol2_file,
        )

        generate_cif_het.cif_writer('./bagpype/dependencies/mmcif/NWU.cif')
        generate_cif_het.HET_writer('./bagpype/dependencies/NWU.txt')

        mol_dir = create_dir(self.output_dir, f'{mol_id}/')

        pdb_file = mol_dir + mol_id + '.pdb'
        generate_cif_het.pdb_writer(pdb_file)

        importlib.reload(bagpype)

        mymol = bagpype.molecules.Protein()
        parser = bagpype.parsing.PDBParser(pdb_file)

        parser.parse(
            mymol,
            strip=strip,
            trim_H=trim_H,
            add_H=add_H,
            MakeMultimer_number=MakeMultimer_number,
        )

        ggenerator = bagpype.construction.Graph_constructor()

        atoms_file_name = mol_dir + mol_id + '_atoms.csv'
        bonds_file_name = mol_dir + mol_id + '_bonds.csv'

        ggenerator.construct_graph(
            mymol,
            mol_id,
            atoms_file_name=atoms_file_name,
            bonds_file_name=bonds_file_name,
        )

        self.mymol_dict[mol_id] = mymol


def get_id(file):
        return pathlib.Path(file).stem.split('.')[0]
