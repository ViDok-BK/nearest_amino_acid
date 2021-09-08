
def build_la_ra_mapping_by_idx(list_receptors, list_ligands, nearest_amino_atom_idx):
    mapping_la_ra = []
    for j, ligand in enumerate(list_ligands):
        mapping_la_ra.append([])
        for i, receptor in enumerate(list_receptors):
            mapping_la_ra[j].append([])

            for list_idx_each_ligand in nearest_amino_atom_idx[i][j]:
                list_receptor_atom = list_receptors[i].get_atoms_by_idx(list_idx_each_ligand)
                mapping_la_ra[j][i].append(list_receptor_atom)
    return mapping_la_ra
