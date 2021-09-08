import numpy as np

def find_candidate_amino_atom(receptor_coords, ligand_coords, model):
    list_result = []
    for i, receptor_coord in enumerate(receptor_coords):
        list_result.append([])
        model.fit(receptor_coord)

        for ligand_coord in ligand_coords:
            # list_distance = np.zeros((10,))
            distances, nearest_ligand_atom_idx = model.radius_neighbors(ligand_coord)
            list_result[i].append(nearest_ligand_atom_idx)

            # for distance in distances:
            #     list_distance += np.array(list(sorted(distance))[:10])
            # print(list_distance / len(distances))

    return list_result
