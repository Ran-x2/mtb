#!/bin/bash
#SBATCH --account=PDS0353
#SBATCH --job-name=rxr456_geneformer_mtb_train
#SBATCH --output=geneformer_mtb_train_output.txt
#SBATCH --error=geneformer_mtb_train_error.txt
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=48
#SBATCH --time=10:00:00
#SBATCH --mem=64gb
#SBATCH --mail-user=rxr456@case.edu
#SBATCH --mail-type=ALL

# Use the Python from your Miniconda installation
export PATH="/users/PDS0353/rxr456/miniconda3/bin:$PATH"
source /users/PDS0353/rxr456/miniconda3/etc/profile.d/conda.sh
storage_dir=/users/PDS0353/rxr456/mtb_250409
conda activate geneformer312
python $storage_dir/geneformer_hpc_1.py