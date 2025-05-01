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

storage_dir = '/mnt/vstor/SOM_PATH_DKB50/members/rxr456/mtb_250409/'
cell_states_to_model = {
    "state_key": "identity", 
    "start_state": "more_expanded_mtb_specific", 
    "goal_state": "more_expanded_bystander"
}
ispstats = InSilicoPerturberStats(mode="goal_state_shift",
                                  genes_perturbed="all",
                                  combos=0,
                                  anchor_gene=None,
                                  cell_states_to_model=cell_states_to_model)

ispstats.get_stats(
    f"{storage_dir}",
    None,
    f"{storage_dir}",
    "2states_overexpression"
)