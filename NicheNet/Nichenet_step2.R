library(nichenetr) # Please update to v2.0.4
library(Seurat)
library(SeuratObject)
library(tidyverse)
library(ggplot2)
library(car)
library(DESeq2)
library(dplyr)
library(EnhancedVolcano)
# BiocManager::install("GenomicRanges")
# BiocManager::install("rtracklayer")
library(GenomicRanges)
library(rtracklayer)

setwd('E:/AAA_Labwork/capenterlab_mtb')
srat = readRDS('srat_subset_for_DEG.rds')
lr_network <- readRDS('F:/CRC/Nichenet/lr_network_human_21122021.rds')
ligand_target_matrix <- readRDS('F:/CRC/Nichenet/ligand_target_matrix_nsga2r_final.rds')
weighted_networks <- readRDS('F:/CRC/Nichenet/weighted_networks_nsga2r_final.rds')
# gr_network = readRDS(url("https://zenodo.org/record/7074291/files/gr_network_human_21122021.rds"))
# sig_network = readRDS(url("https://zenodo.org/record/7074291/files/signaling_network_human_21122021.rds"))

bg_filename = 'expressed_genes_in_less_expanded.csv'
deg_filename = 'expansion_DEG.csv'
    
# ligands influence d3 compare to d4
background_expressed_genes<- read.table(bg_filename, sep = ',', header = 1, row.names = 1)
potential_ligands <- lr_network %>% filter(to %in% background_expressed_genes$X0) %>% pull(from) %>% unique()
  
DE_table_receiver = read.table(deg_filename, sep = ',', header = 1, row.names = 1)
geneset_oi <- DE_table_receiver %>% filter(q_value <= 0.05 & normalized_effect >= 0.5) %>% pull(gene_short_name)
geneset_oi <- geneset_oi %>% .[. %in% rownames(ligand_target_matrix)]
ligand_activities <- predict_ligand_activities(geneset = geneset_oi,
                                                 background_expressed_genes = background_expressed_genes$X0,
                                                 ligand_target_matrix = ligand_target_matrix,
                                                 potential_ligands = potential_ligands)
  
sorted_table <- ligand_activities[order(ligand_activities$auroc,decreasing = TRUE), ]
write.table(sorted_table, paste0('nichenet_result.csv'), row.names = FALSE)
  
top_15_ligands = ligand_activities$test_ligand[1:15]
target_gene_union <- as.data.frame( matrix(NA, nrow = 5, ncol = 15))
colnames(target_gene_union) <- top_15_ligands
for (ligand in top_15_ligands) {
    top_5_targets =ligand_target_matrix[,ligand][order(ligand_target_matrix[,ligand],decreasing = TRUE) ][1:5]
    target_gene_union[ligand] = names(top_5_targets)
  }
all_targets <- unique(as.vector(as.matrix(target_gene_union)))
write.table(all_targets, 'nichenet_targets.csv', row.names = FALSE)  

power_mat_to_plot = ligand_target_matrix[all_targets,top_15_ligands]

srat <- NormalizeData(srat, normalization.method = "LogNormalize", scale.factor = 10000)
Idents(srat) = srat$expansion
subset_srat <- subset(srat, idents = 'more')

dotplot = DotPlot(subset_srat, features = all_targets, assay = 'RNA') +
  RotatedAxis() +
  theme(axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("downstream_in_more_dotplot.jpg", plot = dotplot, width = 12, height = 4,dpi = 300)

divisors <- unlist(apply(power_mat_to_plot, 2, max))
normalized_power_mat_to_plot <- sweep(power_mat_to_plot, 2, divisors, "/")

htmap = pheatmap(t(normalized_power_mat_to_plot),# Scale by row (optional, for better visualization)
                 cluster_rows = FALSE,   # Cluster rows
                 cluster_cols = FALSE,   # Cluster columns
                 color = colorRampPalette(c("blue", "white", "red"))(100),
                 show_rownames = TRUE,
                 show_colnames = TRUE)
ggsave("ligand_powermap.jpg", plot = htmap, width = 12, height = 5,dpi = 300)

#How's the downstream of these ligand expressed in the target cell type ('more')