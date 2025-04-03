#devtools::install_github('cole-trapnell-lab/monocle3')
library(monocle3)
library(Seurat)
library(SingleCellExperiment)

#setwd('/mnt/vstor/SOM_PATH_DKB50/members/rxr456/')
srat = readRDS('srat_subset_for_DEG.rds')
#srat = srat_sub
cds <- as.cell_data_set(srat)
cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(srat[["RNA"]])
  
# Run DE analysis
gene_fits <- fit_models(cds, model_formula_str = "~expansion",cores = 5)
fit_coefs <- coefficient_table(gene_fits)

# Adjust condition in the filter
terms <- fit_coefs %>% filter(term == 'DonorIDDonor AJG2309')
terms = terms %>% select(gene_short_name, term, q_value, normalized_effect)

write.csv(terms,'expansion_DEG.csv')