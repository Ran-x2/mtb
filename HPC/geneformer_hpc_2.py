import multiprocess as mp
mp.set_start_method('fork', force=True)
from multiprocessing import freeze_support

import sys
import os
# pip install protobuf --target=/home/rxr456/miniconda3/envs/geneformer/lib/python3.10/site-packages
# sys.path.append(os.getcwd())
from geneformer import InSilicoPerturber
from geneformer import InSilicoPerturberStats
from geneformer import EmbExtractor
import pickle

print("loading modules complete")

storage_dir = '/mnt/vstor/SOM_PATH_DKB50/members/rxr456/mtb_250409'
output_prefix="cell_state_shift_1goal2alt"
vanilla_model = "/users/PDS0353/rxr456/Geneformer/gf-12L-95M-i4096"
os.chdir(storage_dir)
model = f"{storage_dir}/250409_geneformer_cellClassifier_{output_prefix}/ksplit1/"

cell_states_to_model = {
    "state_key": "identity", 
    "start_state": "less_expanded", 
    "goal_state": "more_expanded_mtb_specific",
    "alt_states":['more_expanded_bystander','more_expanded_mtb_specific_2']
}

with open(f"{storage_dir}/state_emb_4states.pkl", 'rb') as file:
    state_embs_dict = pickle.load(file)

print("open state emb complete")

# embex = EmbExtractor(model_type="CellClassifier",
#                      num_classes=4, 
#                      max_ncells=50000,
#                      emb_layer=-1, 
#                      summary_stat="exact_mean",  # I don't want this stat
#                      forward_batch_size=16,
#                      nproc=80)

# state_embs_dict = embex.get_state_embs(
#     cell_states_to_model,
#     model,
#     f"{storage_dir}/tokenized.dataset",
#     f"{storage_dir}/",
#     "state_emb_4states"
# )

# isp = InSilicoPerturber(perturb_type="overexpress",
#                         genes_to_perturb="all",
#                         combos=0,
#                         anchor_gene=None,
#                         model_type="CellClassifier",
#                         num_classes=4,
#                         emb_mode="cls",
#                         cell_states_to_model=cell_states_to_model,
#                         state_embs_dict=state_embs_dict,
#                         max_ncells=2000,
#                         emb_layer=-1,
#                         forward_batch_size=1,
#                         nproc=60)

# print("start running")
# isp.perturb_data(
#     model,
#     f"{storage_dir}/tokenized.dataset",
#     f"{storage_dir}/",
#     "state_emb_4states"
# )


isp = InSilicoPerturber(perturb_type="overexpress",
                    genes_to_perturb="all",
                    combos=0,
                    anchor_gene=None,
                    model_type="CellClassifier",
                    num_classes=4,
                    emb_mode="cls",
                    cell_states_to_model=cell_states_to_model,
                    state_embs_dict=state_embs_dict,
                    max_ncells=5000,
                    emb_layer=-1,
                    forward_batch_size=1,
                    nproc=60)
print("start running")
isp.perturb_data(
    model,
    f"{storage_dir}/tokenized.dataset",
    f"{storage_dir}/",
    "state_emb_4states"
)