import sys
import os
# pip install protobuf --target=/home/rxr456/miniconda3/envs/geneformer/lib/python3.10/site-packages
sys.path.append(os.getcwd())
from geneformer import InSilicoPerturber
from geneformer import InSilicoPerturberStats
from geneformer import EmbExtractor
import pickle

storage_dir = '/mnt/vstor/SOM_PATH_DKB50/members/rxr456/mtb_new/'
output_prefix="cell_state_shift_6000"
vanilla_model = "/home/rxr456/Geneformer/gf-12L-95M-i4096"
model = f"{storage_dir}/250402_geneformer_cellClassifier_{output_prefix}/ksplit1/"

isp = InSilicoPerturber(perturb_type="overexpress",
                        genes_to_perturb="all",
                        combos=0,
                        anchor_gene=None,
                        model_type="CellClassifier",
                        num_classes=2,
                        emb_mode="cls_and_gene",
                        max_ncells=1000,
                        emb_layer=-1,
                        forward_batch_size=128,
                        nproc=80)

isp.perturb_data(
    model,
    f"{storage_dir}/tokenized.dataset",
    f"{storage_dir}/",
    "gene_perturb_1000_128"
)

# with open(f"{storage_dir}/state_emb.pkl", 'rb') as file:
#     state_embs_dict = pickle.load(file)

# print(state_embs_dict)

# cell_states_to_model = {
#     "state_key": "expansion", 
#     "start_state": "less", 
#     "goal_state": "more",
# }
# print(cell_states_to_model)

# isp = InSilicoPerturber(perturb_type="overexpress",
#                         genes_to_perturb="all",
#                         combos=0,
#                         anchor_gene=None,
#                         model_type="CellClassifier",
#                         num_classes=2,
#                         cell_states_to_model=cell_states_to_model,
#                         state_embs_dict=state_embs_dict,
#                         emb_mode="cls",
#                         max_ncells=1000,
#                         emb_layer=-1,
#                         forward_batch_size=128,
#                         nproc=80)

# print("perturbation start")

# isp.perturb_data(
#     model,
#     f"{storage_dir}/tokenized.dataset",
#     f"{storage_dir}/{output_prefix}/",
#     "cell_state_shift_1000"
# )

# ispstats = InSilicoPerturberStats(mode="goal_state_shift",
#                                   cell_states_to_model=cell_states_to_model,
#                                   genes_perturbed="all",
#                                   combos=0,
#                                   anchor_gene=None)

# ispstats.get_stats(
#     f"{storage_dir}/{output_prefix}/",
#     None,
#     f"{storage_dir}/{output_prefix}/",
#     "cell_state_shift_1000"
# )

ispstats = InSilicoPerturberStats(mode="aggregate_gene_shifts",
                                  genes_perturbed="all",
                                  combos=0,
                                  anchor_gene=None)

ispstats.get_stats(
    f"{storage_dir}/",
    None,
    f"{storage_dir}/",
    "gene_perturb_1000"
)