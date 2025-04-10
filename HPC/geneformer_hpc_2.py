import sys
import os
# pip install protobuf --target=/home/rxr456/miniconda3/envs/geneformer/lib/python3.10/site-packages
sys.path.append(os.getcwd())
from geneformer import InSilicoPerturber
from geneformer import InSilicoPerturberStats
from geneformer import EmbExtractor
import pickle

storage_dir = '/mnt/vstor/SOM_PATH_DKB50/members/rxr456/mtb_250409/'
output_prefix="cell_state_shift_1goal2alt"
vanilla_model = "/home/rxr456/Geneformer/gf-12L-95M-i4096"
model = f"{storage_dir}/250409_geneformer_cellClassifier_{output_prefix}/ksplit1/"

with open(f"{storage_dir}/250409_geneformer_cellClassifier_{output_prefix}/{output_prefix}_eval_metrics_dict.pkl", 'rb') as file:
    all_metrics = pickle.load(file)

embex = EmbExtractor(model_type="CellClassifier",
                     num_classes=4, 
                     max_ncells=10000,
                     emb_layer=-1, 
                     emb_label=["identity"],
                     labels_to_plot=["identity"],
                     forward_batch_size=16,
                     nproc=80)


embs = embex.extract_embs(model,
                          f"{storage_dir}/tokenized.dataset",
                          f"{storage_dir}/",
                          output_prefix + "_embeddings_labeled")

embex.plot_embs(embs=embs,
                plot_style="heatmap",
                output_directory=f"{storage_dir}/",
                output_prefix="embeddings_heatmap")


# all_metrics_test = cc.evaluate_saved_model(
#         model_directory=model,
#         id_class_dict_file=f"{storage_dir}/{output_prefix}_id_class_dict.pkl",
#         test_data_file=f"{storage_dir}/{output_prefix}_labeled_test.dataset",
#         output_directory=f"{storage_dir}/",
#         output_prefix=output_prefix + output_prefix,
#     )

cc.plot_conf_mat(
        conf_mat_dict={"Geneformer": all_metrics["conf_matrix"]},
        output_directory=f"{storage_dir}/",
        output_prefix=output_prefix + output_prefix
)

cell_states_to_model = {
    "state_key": "identity", 
    "start_state": "less_expanded", 
    "goal_state": "more_expanded_mtb_specific",
    "alt_states":['more_expanded_bystander','more_expanded_mtb_specific_2']
}

embex = EmbExtractor(model_type="CellClassifier",
                     num_classes=4, 
                     max_ncells=10000,
                     emb_layer=-1, 
                     summary_stat="exact_mean",  # I don't want this stat
                     forward_batch_size=16,
                     nproc=80)

state_embs_dict = embex.get_state_embs(
    cell_states_to_model,
    model,
    f"{storage_dir}/{output_prefix}/tokenized.dataset",
    f"{storage_dir}/{output_prefix}",
    "state_emb_4states"
)

isp = InSilicoPerturber(perturb_type="overexpress",
                        genes_to_perturb="all",
                        combos=0,
                        anchor_gene=None,
                        model_type="CellClassifier",
                        num_classes=4,
                        emb_mode="cls",
                        max_ncells=20000,
                        emb_layer=-1,
                        forward_batch_size=16,
                        nproc=80)

isp.perturb_data(
    cell_states_to_model,
    model,
    f"{storage_dir}/tokenized.dataset",
    f"{storage_dir}/",
    "state_emb_4states"
)