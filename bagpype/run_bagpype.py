import os
import shutil
import glob
import logging
import importlib
import pathlib
import subprocess
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
            output_dir='output/',           # Dir for all outputs to be saved in
            del_output=True,
    ) -> None:

        # Define paths
        self.output_dir = create_dir(check_dirpath(output_dir), replace=del_output)     # Create output dir if not exists


    def init_from_smiles_string(self, smiles, mol_id='1'):
        smi = Chem.MolToSmiles(Chem.AddHs(Chem.MolFromSmiles(smiles)), True)

        self.mol_files = [self.output_dir + mol_id + '.mol']
        self.mol2_files = [self.output_dir + mol_id + '.mol2']

        self.smiles_to_mol_mol2(
            smi=smi,
            mol_id=mol_id,
        )


    def init_from_smiles_file(self, smiles_file):
        pass


    def init_from_smiles_dir(self, smiles_dir):
        pass


    def init_from_csv(self, smiles_col_name):
        pass


    def smiles_to_mol_mol2(self, smi, mol_id):
        smi_file = self.output_dir + mol_id + '.smi'

        with open(smi_file, 'w') as f:
            f.writelines(smi)

        command_mol = f'obabel -ismi {smi_file} -omol -O {self.mol_files[0]} -h --gen3d'.split(' ')
        subprocess.run(command_mol, stdout = subprocess.PIPE)

        command_mol2 = f'obabel -ismi {smi_file} -omol2 -O {self.mol2_files[0]} -h --gen3d'.split(' ')
        subprocess.run(command_mol2, stdout = subprocess.PIPE)


    def init_from_mol_file(self, mol_file, mol2_file):
        self.mol_files = [mol_file]
        self.mol2_files = [mol2_file]


    def init_from_mol_dir(self, mol_dir, mol2_dir):
        self.mol_files = glob.glob(check_dirpath(mol_dir)+'*.mol')
        self.mol2_files = glob.glob(check_dirpath(mol2_dir)+'*.mol2')


    def create_bagpype_pdb(self):

        get_id = lambda file : pathlib.Path(file).stem.split('.')[0]

        mol_dict = {get_id(file): file for file in self.mol_files}
        mol2_dict = {get_id(file): file for file in self.mol2_files}

        # Setup tqdm tracking bar
        pbar = tqdm(total=len(mol_dict))

        for mol_id in mol_dict.keys():

            pbar.set_description(f'Working on {mol_id}')

            self.run_bagpype(
                mol_id=mol_id,
                mol_file=mol_dict[mol_id],
                mol2_file=mol2_dict[mol_id],
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

        myprot = bagpype.molecules.Protein()
        parser = bagpype.parsing.PDBParser(pdb_file)

        parser.parse(
            myprot,
            strip=strip,
            trim_H=trim_H,
            add_H=add_H,
            MakeMultimer_number=MakeMultimer_number,
        )

        ggenerator = bagpype.construction.Graph_constructor()

        atoms_file_name = mol_dir + mol_id + '_atoms.csv'
        bonds_file_name = mol_dir + mol_id + '_bonds.csv'

        ggenerator.construct_graph(
            myprot,
            mol_id,
            atoms_file_name=atoms_file_name,
            bonds_file_name=bonds_file_name,
        )
