from biopandas.pdb import PandasPdb
import os

class Ligand():
    def __init__(self, pdb_file):
        self.pdb_file = pdb_file
        self.ligand_df = PandasPdb().read_pdb(self.pdb_file)

        self.ligand_atoms = self.ligand_df.df["HETATM"]

    def get_coords(self):
        coords = self.ligand_atoms[["x_coord", "y_coord", "z_coord"]]
        return coords

    def get_atoms_by_idx(self, list_idx):
        return self.ligand_atoms.iloc[list_idx]

    def get_name(self):
        return os.path.basename(self.pdb_file)

    def save_atoms(self, saving_file):
        self.ligand_atoms.to_csv(saving_file)
