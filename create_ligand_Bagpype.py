from biopandas.mol2 import PandasMol2
import pandas as pd
import networkx as nx
import rdkit
import subprocess
import numpy as np

class Mol_Mol2_PDB_from_SMILES():
    def __init__(self, smiles):
        self.smiles = rdkit.Chem.MolToSmiles(rdkit.Chem.AddHs(rdkit.Chem.MolFromSmiles(smiles)), True)

    # def gen_mol_file(self, mol_file_dest):
    #     mol = rdkit.Chem.MolFromSmiles(self.smiles)
    #     mol_H = rdkit.Chem.AddHs(mol)
    #     rdkit.Chem.MolToMolFile(mol_H, mol_file_dest)
    
    # def gen_mol2_file(self, mol_file_dest, mol2_file_dest):
    #     command = f'obabel -imol {mol_file_dest} -omol2 -O {mol2_file_dest}'.split(' ')
    #     subprocess.run(command, stdout = subprocess.PIPE)
    
    def gen_mol_mol2_file(self, smi_file_dest, mol_file_dest, mol2_file_dest):
        with open(smi_file_dest, 'w') as smi_file:
            smi_file.writelines(self.smiles)
        command_mol = f'obabel -ismi {smi_file_dest} -omol -O {mol_file_dest} -h --gen3d'.split(' ')
        subprocess.run(command_mol, stdout = subprocess.PIPE)
        command_mol2 = f'obabel -ismi {smi_file_dest} -omol2 -O {mol2_file_dest} -h --gen3d'.split(' ')
        subprocess.run(command_mol2, stdout = subprocess.PIPE)
    
    def gen_pdb_file(self, mol_file_dest, pdb_file_dest):
        molecule = rdkit.Chem.rdmolfiles.MolFromMolFile(mol_file_dest, removeHs = False)
        rdkit.Chem.rdmolfiles.MolToPDBFile(molecule, pdb_file_dest)

        pdb = pdb_file_dest
        with open(pdb) as pdb_f:
            content = pdb_f.readlines()
            pdb_info = [_ for _ in content if 'HETATM' in _]
            pdb_data = [_.split() for _ in content if 'HETATM' in _]
        
        pdb_lines = []
        for p in range(len(pdb_data)):
            p_0 = 'HETATM  '
            if len(pdb_data[p][1]) == 1:
                p_1 = '  ' + pdb_data[p][1] + ' '
                p_2_num = pdb_data[p][1] + '  '
            elif len(pdb_data[p][1]) == 2:
                p_1 = ' ' + pdb_data[p][1] + ' '
                p_2_num = pdb_data[p][1] + ' '
            elif len(pdb_data[p][1]) == 3:
                p_1 = pdb_data[p][1] + ' '
            if len(pdb_data[p][-1]) == 1:
                p_2 = ' ' + pdb_data[p][-1]
            elif len(pdb_data[p][-1]) == 2:
                p_2 = pdb_data[p][-1]
            p_rest = pdb_info[p][17:]
            p_line = p_0 + p_1 + p_2 + p_2_num + p_rest
            pdb_lines.append(p_line)

        pdb_lines_updated = []
        for p_line in pdb_lines:
            pdb_lines_updated.append(p_line.replace('UNL', 'NWU'))

        with open(pdb_file_dest, 'w') as pdb_file:
            pdb_file.writelines(pdb_lines_updated)

class Ligand_generation():
    def __init__(self, mol_file, mol2_file):
        self.mol = mol_file
        self.mol2 = PandasMol2().read_mol2(mol2_file)
        self.atoms_id = [self.mol2.df['atom_type'].values[i].split('.')[0] + str(self.mol2.df['atom_id'].values[i]) for i in range(len(self.mol2.df))]
        self.type_symbol = np.array([self.mol2.df['atom_type'].values[i].split('.')[0] for i in range(len(self.mol2.df['atom_name']))])
        self.pdbx_ordinal = self.mol2.df['atom_id'].values
        self.pdbx_aromatic_flag = ['Y' if _.endswith('ar') else 'N' for _ in self.mol2.df['atom_type'].values ]
        self.Cartn_x = self.mol2.df['x'].values
        self.Cartn_y = self.mol2.df['y'].values
        self.Cartn_z = self.mol2.df['z'].values
        self.cif_dict = self.cif_dictionary()
        self.cif_df = pd.DataFrame(self.cif_dictionary(), dtype =str)
        self.bonds_data = self.bonds_info()
        self.mol_G = self.neighbouring_atoms()
        self.CONECT_dict = self.CONECT_entries()

    # This extracts information on bonds (atoms and bond order) from the .mol file of the ligand
    def bonds_info(self):
        with open(self.mol) as f:
            # 7 for MOL file from Avogadro, 4 for MOL file from RDKit and 3 for MDL file from CCDC Python API
            mol_data = []
            for line in f.read().splitlines():
                if 'M' not in line and 'Gaussian' not in line and 'Schrodinger' not in line and (len(line.split()) == 7 or len(line.split()) == 4 or len(line.split()) == 3):
                        mol_data.append(line.split()[:3])

        bonds_data = []
        bond_num = 1
        for bond in mol_data:
            bond_info = []
            atom_1 = self.cif_df[self.cif_df['pdbx_ordinal'] == bond[0]].atom_id.values[0]
            bond_info.append(atom_1)
            atom_2 = self.cif_df[self.cif_df['pdbx_ordinal'] == bond[1]].atom_id.values[0]
            bond_info.append(atom_2)
            if bond[2] == '1':
                bond_info.append('SING')
            elif bond[2] == '2':
                bond_info.append('DOUB')
            elif bond[2] == '3':
                bond_info.append('TRIP')
            if self.cif_df[self.cif_df['atom_id'] == atom_1].pdbx_aromatic_flag.values[0] == 'Y' and self.cif_df[self.cif_df['atom_id'] == atom_2].pdbx_aromatic_flag.values[0] == 'Y':
                bond_info.append('Y')
            else:
                bond_info.append('N')
            bond_info += ['N', str(bond_num)]
            bonds_data.append(bond_info)
            bond_num += 1
        return(bonds_data)

    # This combines all the information required from .mol2 file to write a cif file
    def cif_dictionary(self):
        cif_dict = {'comp_id': ['NWU'] * len(self.atoms_id),
           'atom_id': self.atoms_id,
           'alt_atom_id': self.atoms_id,
           'type_symbol': self.type_symbol,
           'charge': ['0'] * len(self.atoms_id),
           'pdbx_align': ['1'] * len(self.atoms_id),
           'pdbx_aromatic_flag': self.pdbx_aromatic_flag,
           'pdbx_leaving_atom_flag': ['N'] * len(self.atoms_id),
           'pdbx_stereo_config': ['N'] * len(self.atoms_id),
           'model_Cartn_x': self.Cartn_x,
           'model_Cartn_y': self.Cartn_y,
           'model_Cartn_z': self.Cartn_z,
           'pdbx_model_Cartn_x_ideal': self.Cartn_x,
           'pdbx_model_Cartn_y_ideal': self.Cartn_y,
           'pdbx_model_Cartn_z_ideal': self.Cartn_z,
           'pdbx_component_atom_id': self.atoms_id,
           'pdbx_component_comp_id': ['NWU'] * len(self.atoms_id),
           'pdbx_ordinal': self.pdbx_ordinal}
        return(cif_dict)

    # The following two functions put the atoms and bonds information in correct format for a .cif file
    def cif_atom_info(self):
        atom_lines = []
        for i in range(len(self.atoms_id)):
            l = [self.cif_dict[_][i] for _ in self.cif_dict.keys()]
            l_0 = l[0] + ' '
            l_1 = list(' ' * 5)
            l_1.insert(0, l[1])
            l_1 = ''.join(l_1)[:5]
            l_2 = l_1
            l_3 = list(' ' * 3)
            l_3.insert(0, l[3])
            l_3 = ''.join(l_3)[:3]
            l_4_8 = ' '.join(l[4:9]) + ' '
            l_9 = list(' ' * 8)
            l_9.insert(0, str('%.3f' % l[9]))
            l_9 = ''.join(l_9)[:8]
            l_10 = list(' ' * 8)
            l_10.insert(0, str('%.3f' % l[10]))
            l_10 = ''.join(l_10)[:8]
            l_11 = list(' ' * 8)
            l_11.insert(0, str('%.3f' % l[11]))
            l_11 = ''.join(l_11)[:8]
            l_12_14 = l_9 + l_10 + l_11
            l_15 = l_1
            l_16 = l_0
            l_17 = list(' ' * 3)
            l_17.insert(0, str(l[17]))
            l_17 = ''.join(l_17)[:3]
            atom_line = l_0 + l_1 + l_2 + l_3 + l_4_8 + l_9 + l_10 + l_11 + l_12_14 + l_15 + l_16 + l_17 + '\n'
            atom_lines.append(atom_line)
        return(atom_lines)

    def cif_bond_info(self):
        bond_lines = []
        for b in self.bonds_data:
            b_0 = 'NWU '
            b_1 = list(' ' * 5)
            b_1.insert(0, b[0])
            b_1 = ''.join(b_1)[:5]
            b_2 = list(' ' * 5)
            b_2.insert(0, b[1])
            b_2 = ''.join(b_2)[:5]
            b_3_5 = ' '.join(b[2:5]) + ' '
            b_6 = list(' ' * 3)
            b_6.insert(0, b[5])
            b_6 = ''.join(b_6)[:3]
            bond_line = b_0 + b_1 + b_2 + b_3_5 + b_6 + '\n'
            bond_lines.append(bond_line)
        return(bond_lines)

    # This writes out the .cif file required for Bagpype
    def cif_writer(self, cif_file_dest):
        with open(cif_file_dest, 'w') as cif_file:
            cif_file.writelines(['data_NWU\n', 
                                '# \n', 
                                '_chem_comp.id                                    NWU \n',
                                '_chem_comp.type                                  NON-POLYMER \n',
                                '#\nloop_\n',
                                '_chem_comp_atom.comp_id \n',
                                '_chem_comp_atom.atom_id \n',
                                '_chem_comp_atom.alt_atom_id \n',
                                '_chem_comp_atom.type_symbol \n',
                                '_chem_comp_atom.charge \n',
                                '_chem_comp_atom.pdbx_align \n',
                                '_chem_comp_atom.pdbx_aromatic_flag \n',
                                '_chem_comp_atom.pdbx_leaving_atom_flag \n',
                                '_chem_comp_atom.pdbx_stereo_config \n',
                                '_chem_comp_atom.model_Cartn_x \n',
                                '_chem_comp_atom.model_Cartn_y \n',
                                '_chem_comp_atom.model_Cartn_z \n',
                                '_chem_comp_atom.pdbx_model_Cartn_x_ideal \n',
                                '_chem_comp_atom.pdbx_model_Cartn_y_ideal \n',
                                '_chem_comp_atom.pdbx_model_Cartn_z_ideal \n',
                                '_chem_comp_atom.pdbx_component_atom_id \n',
                                '_chem_comp_atom.pdbx_component_comp_id \n',
                                '_chem_comp_atom.pdbx_ordinal \n'])
            
            cif_file.writelines(self.cif_atom_info())
            
            cif_file.writelines(['# \n', 
                                '#\nloop_\n',
                                '_chem_comp_bond.comp_id \n',
                                '_chem_comp_bond.atom_id_1 \n',
                                '_chem_comp_bond.atom_id_2 \n',
                                '_chem_comp_bond.value_order \n',
                                '_chem_comp_bond.pdbx_aromatic_flag \n',
                                '_chem_comp_bond.pdbx_stereo_config \n',
                                '_chem_comp_bond.pdbx_ordinal \n'])
            
            cif_file.writelines(self.cif_bond_info())

    # This generates the CONECT entries required for Reduce
    def CONECT_entries(self):
        CONECT_dict = []
        for atom in self.atoms_id:
            atom_neighbour = list(self.mol_G.neighbors(atom))
            atom_neighbour_updated = []
            for conect_atom in atom_neighbour:
                if len(conect_atom) == 1:
                    atom_neighbour_updated.append(' ' + conect_atom + '   ')
                elif len(conect_atom) == 2:
                    atom_neighbour_updated.append(' ' + conect_atom + '  ')
                elif len(conect_atom) == 3:
                    atom_neighbour_updated.append(' ' + conect_atom + ' ')
                elif len(conect_atom) == 4:
                    atom_neighbour_updated.append(conect_atom)
            if len(atom) == 1:
                atom_updated = '  ' + atom + ' '
            elif len(atom) == 2:
                atom_updated = ' ' + atom + ' '
            elif len(atom) == 3:
                atom_updated = ' ' + atom
            elif len(atom) == 4:
                atom_updated = atom  
            CONECT_dict.append(([atom_updated] + [str(len(atom_neighbour))] + atom_neighbour_updated))
        return(CONECT_dict)

    # This creates the correct format for HET file used by Reduce
    def CONECT_lines_info(self):
        CONECT_lines = []
        for c in self.CONECT_dict:
            c_0 = 'CONECT     '
            c_1 = c[0] + '    '
            c_2 = c[1]
            c_3_on = ''.join(c[2:])
            c_line = c_0 + c_1 + c_2 + c_3_on + '\n'
            CONECT_lines.append(c_line)
        return(CONECT_lines)

    def HET_writer(self, HET_file_dest):
        with open(HET_file_dest, 'w') as reduce_file:
            reduce_file.writelines([f'RESIDUE   NWU     {len(self.CONECT_dict)}\n'] + self.CONECT_lines_info() + ['END'])

    # This constructs the molecular graph for CONECT entries for Reduce
    def neighbouring_atoms(self):
        mol_G = nx.Graph()
        mol_G.add_nodes_from(self.atoms_id)
        mol_G.add_edges_from([_[:2] for _ in self.bonds_data])
        return(mol_G)

    # This generates the .pdb file of the molecule for Bagpype graph construction
    def pdb_writer(self, pdb_file_dest):
        molecule = rdkit.Chem.rdmolfiles.MolFromMolFile(self.mol, removeHs = False)
        rdkit.Chem.rdmolfiles.MolToPDBFile(molecule, pdb_file_dest)

        pdb = pdb_file_dest
        with open(pdb) as pdb_f:
            content = pdb_f.readlines()
            pdb_info = [_ for _ in content if 'HETATM' in _]
            pdb_data = [_.split() for _ in content if 'HETATM' in _]
        
        pdb_lines = []
        for p in range(len(pdb_data)):
            p_0 = 'HETATM  '
            if len(pdb_data[p][1]) == 1:
                p_1 = '  ' + pdb_data[p][1] + ' '
                p_2_num = pdb_data[p][1] + '  '
            elif len(pdb_data[p][1]) == 2:
                p_1 = ' ' + pdb_data[p][1] + ' '
                p_2_num = pdb_data[p][1] + ' '
            elif len(pdb_data[p][1]) == 3:
                p_1 = pdb_data[p][1] + ' '
            if len(pdb_data[p][-1]) == 1:
                p_2 = ' ' + pdb_data[p][-1]
            elif len(pdb_data[p][-1]) == 2:
                p_2 = pdb_data[p][-1]
            else:
                if len(''.join([_ for _ in list(pdb_data[p][-1]) if _ not in '0123456789+-'])) == 1:
                    p_2 = ' ' + ''.join([_ for _ in list(pdb_data[p][-1]) if _ not in '0123456789+-'])
                else:
                    p_2 = ''.join([_ for _ in list(pdb_data[p][-1]) if _ not in '0123456789+-'])
            p_rest = pdb_info[p][17:]
            p_line = p_0 + p_1 + p_2 + p_2_num + p_rest
            pdb_lines.append(p_line)

        pdb_lines_updated = []
        for p_line in pdb_lines:
            pdb_lines_updated.append(p_line.replace('UNL', 'NWU'))

        with open(pdb_file_dest, 'w') as pdb_file:
            pdb_file.writelines(pdb_lines_updated)
    

