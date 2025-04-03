# srun --nodes=1 --time=2:00:00 --ntasks=6 --mem=200G --pty /bin/bash
#for hpc, .libPaths(c("/home/rxr456/R/x86_64-pc-linux-gnu-library/4.2/", .libPaths()))
# install.packages("stringi", type = "source")
library(scran)
library(Matrix)
library(Seurat)
library(SingleCellExperiment)
library(SeuratWrappers)
library(dplyr)  
library(SeuratDisk)
library(SeuratData)
library(rhdf5)
setwd('E:/AAA_Labwork/capenterlab_mtb')
#setwd('/mnt/vstor/SOM_PATH_DKB50/members/rxr456/')
Convert("SeuratProject.h5Seurat", dest = "h5ad")

H5list = h5ls("SeuratProject.h5Seurat") #there is an error by just reading the
gene_names = h5read("SeuratProject.h5Seurat",'/assays/SCT/features')
write.table(data.frame(gene_names), 'var_names.csv', sep = ',')

metadata = h5read("SeuratProject.h5Seurat",'/meta.data')
cell_names = metadata$`_index`

SCT_corrected_counts = h5read("SeuratProject.h5Seurat",'/assays/SCT/counts')
saveRDS(SCT_corrected_counts,'SCT_corrected_counts')
SCT_corrected_counts = readRDS('SCT_corrected_counts')

n_rows <- 20395
n_cols <- 157462
nonzero_values <- as.numeric(SCT_corrected_counts$data)
row_indices <- SCT_corrected_counts$indices + 1  
sparse_mat <- sparseMatrix(
  i = row_indices,         
  p = SCT_corrected_counts$indptr,  
  x = nonzero_values,    
  dims = c(n_rows, n_cols)
)

colnames(sparse_mat) = as.character(cell_names)
row.names(sparse_mat) = as.character(gene_names)
srat_for_DEG = CreateSeuratObject(counts = sparse_mat, project = "mtb", min.cells = 0, min.features = 0, assay = "RNA")
saveRDS(srat_for_DEG,'srat_for_DEG.rds')

srat_for_DEG = readRDS('srat_for_DEG.rds')
srat_for_DEG@meta.data$seurat_clusters =as.character(metadata$seurat_clusters$values)

clusters_to_keep <- c("0", "3", "6", "12", "14")
srat_sub <- subset(srat_for_DEG, subset = seurat_clusters %in% clusters_to_keep)
srat_sub$expansion <- ifelse(srat_sub$seurat_clusters %in% c("0", "3"), "less", "more")
saveRDS(srat_sub,'srat_subset_for_DEG.rds')

