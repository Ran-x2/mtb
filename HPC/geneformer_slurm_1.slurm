#!/bin/bash
#SBATCH --job-name=rxr456_geneformer_mtb_train
#SBATCH --output=geneformer_mtb_train_output.txt
#SBATCH --error=geneformer_mtb_train_error.txt
#SBATCH --time=10:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=24
#SBATCH --mem=128G
#SBATCH --mail-user=rxr456@case.edu
#SBATCH --mail-type=ALL

# Use the Python from your Miniconda installation
export PATH="/home/rxr456/miniconda3/bin:$PATH"
source /home/rxr456/miniconda3/etc/profile.d/conda.sh
storage_dir=/mnt/vstor/SOM_PATH_DKB50/members/rxr456/mtb_250409
conda activate geneformer
python $storage_dir/geneformer_hpc_1.py