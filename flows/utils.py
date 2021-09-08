from configs import *
import os

def ensure_path(path_name):
    if not os.path.exists(path_name):
        os.mkdir(path_name)

def save_results(result_type=None, **kwargs):
    if not result_type:
        return

    saving_path = os.path.join(results_path, result_type)
    ensure_path(saving_path)

    if result_type == "nearest_atoms":
        for k, v in kwargs.items():
            element_path = os.path.join(saving_path, k)
            ensure_path(element_path)

            if k != "lr_mapping":
                for molecule in v:
                    saving_file = os.path.join(element_path, molecule.get_name() + ".csv")
                    molecule.save_atoms(saving_file)

            else:
                for ligand_idx, per_ligand in enumerate(v):
                    ligand_name = kwargs["ligand_atoms"][ligand_idx].get_name()

                    for receptor_idx, per_receptor in enumerate(per_ligand):
                        receptor_name = kwargs["receptor_atoms"][receptor_idx].get_name()
                        mapping_name = ligand_name + "-" + receptor_name + ".csv"
                        saving_file = os.path.join(element_path, mapping_name)

                        with open(saving_file, "w", encoding="utf-8") as f:
                            for la_idx, latom in enumerate(per_receptor):
                                f.write(str(la_idx)+",")
                                ratoms = [str(x) for x in latom.index]
                                f.write(",".join(ratoms) + "\n")
