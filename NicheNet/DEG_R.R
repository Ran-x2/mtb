# srun --time=8:00:00 --ntasks=48 --mem=128G --pty /bin/bash
module load GDAL
module load GEOS
module load PROJ
module load R

module load gcc/12.3.0
module load gdal/3.7.3
module load geos/3.12.0
module load proj/9.2.1
module load R/4.4.0
module load curl/8.4.0
export LD_LIBRARY_PATH=/users/PDS0353/rxr456/miniconda3/lib:$LD_LIBRARY_PATH
# devtools::install_github('cole-trapnell-lab/monocle3')

# tar -xvf udunits-2.2.28.tar.gz
# cd udunits-2.2.28
# ./configure --prefix=/home/rxr456/udunits2
# make
# make install

export UDUNITS2_INCLUDE=/home/rxr456/udunits2/include
export UDUNITS2_LIBS=/home/rxr456/udunits2/lib

# export LD_LIBRARY_PATH=/home/rxr456/udunits2/lib:$LD_LIBRARY_PATH
# install.packages('sf', configure.args = '--/home/rxr456/udunits2/lib')
# install.packages('sf', configure.args = '--/home/rxr456/udunits2/include')

# install.packages("sf", 
#   configure.args = "--with-udunits2-lib=--/home/rxr456/udunits2/lib --with-udunits2-include=/home/rxr456/udunits2/include")

#remotes::install_github('satijalab/seurat-wrappers')
# if (!require("BiocManager", quietly = TRUE))
#   install.packages("BiocManager")
# 
# BiocManager::install("SingleCellExperiment")
# BiocManager::install("scran")
# remotes::install_github('satijalab/seurat-wrappers')

library(monocle3)
library(Seurat)
library(SingleCellExperiment)
library(SeuratWrappers)
# library(scran)
library(Matrix)
library(dplyr)  
# library(SeuratDisk)
# library(SeuratData)
library(rhdf5)
# setwd('E:/AAA_Labwork/capenterlab_mtb')
setwd('/mnt/vstor/SOM_PATH_DKB50/members/rxr456/mtb_DEG')


# to more expanded bystanders ---------------------------------------------
srat = readRDS('srat.rds')
metadata = h5read("SeuratProject.h5Seurat",'/meta.data')
srat@meta.data$seurat_clusters =as.character(metadata$seurat_clusters$values)

clusters_to_keep <- c("0", "3", "6", "12", "14")
srat_sub <- subset(srat, subset = seurat_clusters %in% clusters_to_keep)
srat_sub$expansion <- ifelse(srat_sub$seurat_clusters %in% c("0", "3"), "less", "more")
saveRDS(srat_sub,'srat_subset_for_DEG_bystanders.rds')

srat = readRDS('srat_subset_for_DEG.rds_bystanders')
#srat = srat_sub
cds <- as.cell_data_set(srat)
cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(srat[["RNA"]])
  
# Run DE analysis
gene_fits <- fit_models(cds, model_formula_str = "~expansion",cores = 20)
saveRDS(gene_fits,'gene_fits_bystanders')

fit_coefs <- coefficient_table(gene_fits)
saveRDS(fit_coefs,'fit_coefs_bystanders')
fit_coefs = readRDS('fit_coefs_bystanders')
# fit_coefs <- coefficient_table(gene_fits)
# Adjust condition in the filter
terms <- fit_coefs %>% filter(term == 'expansionmore')
terms = terms %>% select(gene_short_name, term, q_value, normalized_effect)

write.csv(terms,'expansion_DEG_bystanders.csv')



# To mtb specific 1 -------------------------------------------------------
srat = readRDS('srat.rds')
srat@meta.data$seurat_clusters =as.character(metadata$seurat_clusters$values)

clusters_to_keep <- c("0", "3", "4", "13")
srat_sub <- subset(srat, subset = seurat_clusters %in% clusters_to_keep)
srat_sub$expansion <- ifelse(srat_sub$seurat_clusters %in% c("0", "3"), "less", "more")
saveRDS(srat_sub,'srat_subset_for_DEG_mtb1.rds')

srat = readRDS('srat_subset_for_DEG_mtb1.rds')
#srat = srat_sub
cds <- as.cell_data_set(srat)
cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(srat[["RNA"]])

# Run DE analysis
gene_fits <- fit_models(cds, model_formula_str = "~expansion",cores = 20)
saveRDS(gene_fits,'gene_fits_mtb1')

fit_coefs <- coefficient_table(gene_fits)
saveRDS(fit_coefs,'fit_coefs_mtb1')
fit_coefs = readRDS('fit_coefs_mtb1')
# fit_coefs <- coefficient_table(gene_fits)
# Adjust condition in the filter
terms <- fit_coefs %>% filter(term == 'expansionmore')
terms = terms %>% select(gene_short_name, term, q_value, normalized_effect)

write.csv(terms,'expansion_DEG_mtb1.csv')
print('all done')

# mtb specific 1 to bystanders-------------------------------------------------------
srat = readRDS('srat.rds')
srat@meta.data$seurat_clusters =as.character(metadata$seurat_clusters$values)

clusters_to_keep <- c("4", "13","6", "12", "14")
srat_sub <- subset(srat, subset = seurat_clusters %in% clusters_to_keep)
srat_sub$expansion <- ifelse(srat_sub$seurat_clusters %in% c("4", "13"), "specific", "bystander")
saveRDS(srat_sub,'srat_subset_for_DEG_mtb1_to_bystander.rds')

srat = readRDS('srat_subset_for_DEG_mtb1_to_bystander.rds')
#srat = srat_sub
cds <- as.cell_data_set(srat)
cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(srat[["RNA"]])

# Run DE analysis
gene_fits <- fit_models(cds, model_formula_str = "~expansion",cores = 30)
# saveRDS(gene_fits,'gene_fits_mtb_bystander')

fit_coefs <- coefficient_table(gene_fits)
saveRDS(fit_coefs,'fit_coefs_mtb_bystander')
# fit_coefs = readRDS('fit_coefs_mtb_bystander')
# fit_coefs <- coefficient_table(gene_fits)
# Adjust condition in the filter
terms <- fit_coefs %>% filter(term == 'expansionbystander')
terms = terms %>% select(gene_short_name, term, q_value, normalized_effect)

write.csv(terms,'expansion_DEG_mtb_bystander.csv')
print('all done')

