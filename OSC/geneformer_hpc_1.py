import sys
import os
from geneformer import InSilicoPerturber
from geneformer import InSilicoPerturberStats
from geneformer import EmbExtractor
import pickle

storage_dir = '/users/PDS0353/rxr456/mtb_250409'
try:
    os.makedirs(storage_dir)
except:
    print('path already exist!')

output_prefix="mtb_250409"
vanilla_model = "/users/PDS0353/rxr456/Geneformer/gf-12L-95M-i4096"
os.chdir(storage_dir)

# Did that during testing so commented out
# from geneformer import TranscriptomeTokenizer
# tk = TranscriptomeTokenizer({"identity": "identity"}, nproc=45)
# tk.tokenize_data(f"{storage_dir}", 
#                  f"{storage_dir}",
#                  "tokenized", 
#                  file_format="h5ad")

from geneformer import Classifier
cc = Classifier(classifier="cell",
                cell_state_dict = {"state_key": "identity", "states": "all"},
                max_ncells=None,
                freeze_layers = 6,
                num_crossval_splits = 1,
                split_sizes = {"train": 0.6, "valid": 0.2, "test": 0.2},
                forward_batch_size=16,
                nproc=47)


cc.prepare_data(input_data_file=f"{storage_dir}/tokenized.dataset",
                output_directory=f"{storage_dir}/",
                output_prefix=output_prefix)

all_metrics = cc.validate(model_directory=vanilla_model,
                          prepared_input_data_file=f"{storage_dir}/{output_prefix}_labeled_train.dataset",
                          id_class_dict_file=f"{storage_dir}/{output_prefix}_id_class_dict.pkl",
                          output_directory=f"{storage_dir}/",
                          output_prefix=output_prefix,
                          #n_hyperopt_trials=1,
                          predict_eval=True)

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