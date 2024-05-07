import os
import shutil
import glob
import logging
import importlib
import pathlib

import create_ligand_Bagpype as create_files
from bagpype.helper_funcs import create_dir, check_dirpath


# Setup logger
logger = logging.getLogger(__name__)
logger.setLevel(logging.INFO)
formatter = logging.Formatter('%(levelname)s:%(asctime)s:%(name)s:%(message)s')

# Define handlers
# fh = logging.FileHandler('{}.log'.format(__name__))
fh = logging.FileHandler('batch_run_bagpype.log')
fh.setFormatter(formatter)
sh = logging.StreamHandler()

# Add handlers to logger
logger.addHandler(fh)
logger.addHandler(sh)


def run_bagpype(
    mol_id,
    mol_file,
    mol2_file,
    output_dir='./',
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

    mol_dir = create_dir(output_dir, f'{mol_id}/')

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
        atoms_file_name=atoms_file_name,
        bonds_file_name=bonds_file_name,
    )


def batch_run_bagpype(mol_dir, mol2_dir, output_dir):

    mol_files = glob.glob(check_dirpath(mol_dir)+'*.mol')
    mol2_files = glob.glob(check_dirpath(mol2_dir)+'*.mol2')

    get_id = lambda file : pathlib.Path(file).stem.split('.')[0]

    mol_dict = {get_id(file): file for file in mol_files}
    mol2_dict = {get_id(file): file for file in mol2_files}

    for mol_id in mol_dict.keys():
        run_bagpype(
            mol_id=mol_id,
            mol_file=mol_dict[mol_id],
            mol2_file=mol2_dict[mol_id],
            output_dir=output_dir,
        )
