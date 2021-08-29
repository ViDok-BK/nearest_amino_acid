import numpy as np
import pandas as pd
from biopandas.pdb import PandasPdb as pdpdb

"""
Version 1.1

We use extract_pdb function to extract coordinates and chemical infomation from the raw data.
Coordinates are defined in R3, containing a tuple (x,y,z).
Chemical information of an atom includes "name", "residue name", "chain id". Particularly, residue name of each atom in virus is the amino acid's name, and is "LIG" if 
atoms are in ligand.

When we have a matrix of distances and indices after using knn, we need to compare with the matrix of distances with
the threshold to sort out the final matrix of indices then return the final matrix of distances and indices (after sort out).

"""

def extract_pdb(virus_frame, lig_frame):
    """
    Input: pandas frame of virus, ligand.
    Output: numpy arrays containing coordinates, chemical informations, label(for feeding in knn).
    """
    frame_virus = virus_frame.df["ATOM"]
    frame_lig = lig_frame.df["HETATM"]
    train_info, x_train = frame_virus[["atom_name", "residue_name", "chain_id", "residue_number"]].to_numpy(), frame_virus[["x_coord", "y_coord", "z_coord"]].to_numpy()
    test_info, x_test = frame_lig[["atom_name", "residue_name", "chain_id"]].to_numpy(), frame_lig[["x_coord", "y_coord", "z_coord"]].to_numpy()
    """
    We just use the coordinates (x) to calculate distance; however, we need to pass (train, label) in knn, so assign a vector of ones for the virus,
    and vector of zeros for the ligand.
    """
    y_train = np.zeros(len(x_train), dtype = np.int32) 
    y_test = np.ones(len(x_test), dtype = np.int32)
    return x_train, y_train, x_test, y_test, train_info, test_info

def extract_from_threshold(distances, indices, threshold):
    """
    Input:  Initial matrix of distances and indices, distance threshold.

    """
    bool_matrix = distances <= threshold
    matrix_distance_dest = distances[bool_matrix]
    matrix_indices_dest = indices[bool_matrix]
    return matrix_distance_dest, matrix_indices_dest
