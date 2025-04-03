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
Convert("SeuratProject.h5Seurat", dest = "h5ad")

H5list = h5ls("SeuratProject.h5Seurat") #there is an error by just reading the
gene_names = h5read("SeuratProject.h5Seurat",'/assays/SCT/features')
write.table(data.frame(gene_names), 'var_names.csv', sep = ',')

SCT_corrected_counts = h5read("SeuratProject.h5Seurat",'/assays/SCT/counts')
saveRDS(SCT_corrected_counts,'SCT_corrected_counts')
