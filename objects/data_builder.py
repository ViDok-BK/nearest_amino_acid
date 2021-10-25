from objects.receptor import Receptor
from objects.ligand import Ligand
import os

class DataBuilder():
    @classmethod
    def init(cls, data_dir):
        cls.data_dir = data_dir
        cls.receptor_dir = os.path.join(cls.data_dir, "receptors")
        cls.ligand_dir = os.path.join(cls.data_dir, "ligands")

    @classmethod
    def check_file_type(cls, filename):
        fex = filename.split(".")[-1]
        return fex in ["pdb", "pdbqt"]

    @classmethod
    def get_atom_data(cls, only_coords=False):
        list_receptors = []
        list_ligands = []

        for rfile in os.listdir(cls.receptor_dir):
            if not cls.check_file_type(rfile):
                continue
            list_receptors.append(
                Receptor(os.path.join(cls.receptor_dir, rfile))
            )

        for lfile in os.listdir(cls.ligand_dir):
            if not cls.check_file_type(lfile):
                continue
            list_ligands.append(
                Ligand(os.path.join(cls.ligand_dir, lfile))
            )

        if only_coords:
            list_receptors = list(map(lambda x: x.get_coords(), list_receptors))
            list_ligands = list(map(lambda x: x.get_coords(), list_ligands))

        return list_receptors, list_ligands
