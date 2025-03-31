import sys
import os
sys.path.append(os.getcwd())
from geneformer import InSilicoPerturber
from geneformer import InSilicoPerturberStats
from geneformer import EmbExtractor
import pickle

storage_dir = '/mnt/vstor/SOM_PATH_DKB50/members/rxr456/mtb'
output_prefix="cell_state_shift_2000"
vanilla_model = "/home/rxr456/Geneformer/gf-12L-95M-i4096"
model = f"{storage_dir}/250324_geneformer_cellClassifier_mtb/ksplit1/"

with open(f"{storage_dir}/state_emb.pkl", 'rb') as file:
    state_embs_dict = pickle.load(file)

print(state_embs_dict)

cell_states_to_model = {
    "state_key": "expansion", 
    "start_state": "less", 
    "goal_state": "more",
}
print(cell_states_to_model)

isp = InSilicoPerturber(perturb_type="overexpress",
                        genes_to_perturb="all",
                        combos=0,
                        anchor_gene=None,
                        model_type="CellClassifier",
                        num_classes=2,
                        cell_states_to_model=cell_states_to_model,
                        state_embs_dict=state_embs_dict,
                        emb_mode="cls",
                        max_ncells=2000,
                        emb_layer=-1,
                        forward_batch_size=128,
                        nproc=80)

print("perturbation start")

isp.perturb_data(
    model,
    f"{storage_dir}/tokenized.dataset",
    f"{storage_dir}/{output_prefix}/",
    "cell_state_shift_2000"
)