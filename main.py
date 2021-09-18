from biopandas.pdb import PandasPdb as pdpdb
import numpy as np
import pandas as pd
from util_functions import extract_pdb 
from util_functions import extract_from_threshold
from sklearn.neighbors import KNeighborsClassifier as knn
import matplotlib.pyplot as plt
import time


#Define some hyper-parameters
virus_path = "/home/nhqcs/Desktop/Github/nearest_amino_axit/Receptor_vidok.pdb"
lig_path = "/home/nhqcs/Desktop/Github/nearest_amino_axit/6lu7_Ligand.pdb"
threshold = 4.0
k_neighbors = 20
amino_group = []

if __name__ == '__main__':
    begin = time.perf_counter()

    virus_init_frame = pdpdb().read_pdb(virus_path)
    lig_init_frame = pdpdb().read_pdb(lig_path)
    x_train, y_train, x_test, y_test, train_info, test_info = extract_pdb(virus_init_frame, lig_init_frame)
   
    #Define KNN algorithm.
    model = knn(n_neighbors= k_neighbors, algorithm = 'kd_tree', metric = "euclidean")
    model.fit(x_train, y_train)
    neighs = model.kneighbors(x_test, return_distance = True)
    amino_group = []
    for i in range(len(x_test)):

    #Extract distances and indices respect to the threshold
        matrix_distance_dest, matrix_indices_dest = extract_from_threshold(neighs[0][i], neighs[1][i], 8.0)
        tmp = [train_info[x] for x in matrix_indices_dest]
        amino_group.extend(tmp)

    headings = ['atom_name', 'residue_name', 'chain_id', 'residue_number']
    _result = pd.DataFrame(amino_group, columns=headings).to_csv("result.csv")

    #Read the csv to write out the final result, including the amino acid's name and its number
    #Ex: GLU166, PHE170, etc.
    end = time.perf_counter()

    print("Elapsed time: "+ str(end-begin))