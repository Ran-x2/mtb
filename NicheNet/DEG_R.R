# srun --nodes=1 --time=1:00:00 --ntasks=48 --mem=128G --pty /bin/bash
# module load GDAL
# module load GEOS
# module load PROJ
# devtools::install_github('cole-trapnell-lab/monocle3')

# tar -xvf udunits-2.2.28.tar.gz
# cd udunits-2.2.28
# ./configure --prefix=/home/rxr456/udunits2
# make
# make install

# export UDUNITS2_INCLUDE=/home/rxr456/udunits2/include
# export UDUNITS2_LIBS=/home/rxr456/udunits2/lib
# export LD_LIBRARY_PATH=/home/rxr456/udunits2/lib:$LD_LIBRARY_PATH
# install.packages('sf', configure.args = '--/home/rxr456/udunits2/lib')
# install.packages('sf', configure.args = '--/home/rxr456/udunits2/include')

# install.packages("sf", 
#   configure.args = "--with-udunits2-lib=--/home/rxr456/udunits2/lib --with-udunits2-include=/home/rxr456/udunits2/include")


library(monocle3)
library(Seurat)
library(SingleCellExperiment)
setwd('/mnt/vstor/SOM_PATH_DKB50/members/rxr456/')
srat = readRDS('srat_subset_for_DEG.rds')
#srat = srat_sub
cds <- as.cell_data_set(srat)
cds <- estimate_size_factors(cds)
cds@rowRanges@elementMetadata@listData[["gene_short_name"]] <- rownames(srat[["RNA"]])
  
# Run DE analysis
gene_fits <- fit_models(cds, model_formula_str = "~expansion",cores = 40)
fit_coefs <- coefficient_table(gene_fits)

# Adjust condition in the filter
terms <- fit_coefs %>% filter(term == 'DonorIDDonor AJG2309')
terms = terms %>% select(gene_short_name, term, q_value, normalized_effect)

write.csv(terms,'expansion_DEG.csv')