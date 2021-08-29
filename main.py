from biopandas.pdb import PandasPdb as pdpdb
import numpy as np
import pandas as pd
from util_functions import extract_pdb 
from util_functions import extract_from_threshold
from sklearn.neighbors import KNeighborsClassifier as knn
import matplotlib.pyplot as plt
import time


#Define some hyper-parameters
virus_path = "C:\\Users\\Admin\\OneDrive\\Desktop\\Github\\Vidok\\Receptor_vidok.pdb"
lig_path = "C:\\Users\\Admin\\OneDrive\\Desktop\\Github\\Vidok\\6lu7_Ligand.pdb"
threshold = 4.0
k_neighbors = 20
amino_group = []

if __name__ == '__main__':
    begin = time.time()
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
    """
    path = 

    df = pd.read_csv(path)

    frame = df.groupby(["residue_name", "residue_number"])
    with open("result.txt", 'w') as f:
        for item in frame:
            f.write(item[0][0] + "")
            f.write(str(item[0][1])+ '\n')
    end = time.time()
    print(end-begin)  
    """