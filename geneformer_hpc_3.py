import sys
import os
# pip install protobuf --target=/home/rxr456/miniconda3/envs/geneformer/lib/python3.10/site-packages
sys.path.append(os.getcwd())
from geneformer import InSilicoPerturberStats
import pickle

storage_dir = '/mnt/vstor/SOM_PATH_DKB50/members/rxr456/mtb_new/'

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