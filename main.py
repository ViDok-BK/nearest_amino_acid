from flows.candidate_link import find_candidate_amino_atom
from flows.post_process import build_la_ra_mapping_by_idx
from flows.utils import save_results
from objects.data_builder import DataBuilder
from models.configs import model_configs
from models import create_model
from utils import *
import configs
import time
import sys

@record_timestamp
def main(FLAGS=None):
    if FLAGS is None:
        FLAGS, _ = configs.parser.parse_known_args(args=sys.argv[1:])

    # Build data
    print("BUILDING DATA...")
    DataBuilder().init(configs.data_path)
    list_receptors, list_ligands = DataBuilder().get_atom_data()
    receptor_coords, ligand_coords = DataBuilder().get_atom_data(only_coords=True)

    # Initialize models
    print("CREATING MODELS...")
    model = create_model(FLAGS.model, **model_configs[FLAGS.model])

    # Run program
    # start_time = time.time()
    print("RUNNING ALGORITHM...")
    nearest_amino_atom_idx = find_candidate_amino_atom(receptor_coords, ligand_coords, model)

    # Post-process result
    print("POST-PROCESSING RESULTS...")
    la_ra_mapping = build_la_ra_mapping_by_idx(list_receptors, list_ligands, nearest_amino_atom_idx)
    # print(time.time() - start_time)

    # Save result
    print("SAVING RESULTS...")
    save_results("nearest_atoms", receptor_atoms=list_receptors,
                                  ligand_atoms=list_ligands,
                                  lr_mapping=la_ra_mapping)

    # Visualize result
    count_ratoms_per_latom = []
    for ligand_idx, per_ligand in enumerate(la_ra_mapping):
        for receptor_idx, per_receptor in enumerate(per_ligand):
            for la_idx, latom in enumerate(per_receptor):
                count_ratoms_per_latom.append(len(latom))

    print("Average: %.02f" % (sum(count_ratoms_per_latom)/len(count_ratoms_per_latom)))
    print("Total: %d" % sum(count_ratoms_per_latom))

    print("DONE!")

if __name__ == '__main__':
    _, runtime = main()
    print("Total runtine: %.03f" % runtime)
