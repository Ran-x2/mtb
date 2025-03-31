#!/bin/bash
#SBATCH --job-name=rxr456_geneformer_mtb_2000
#SBATCH --output=geneformer_mtb_2000_output.txt
#SBATCH --error=geneformer_mtb_2000_error.txt
#SBATCH --time=150:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --cpus-per-task=24
#SBATCH --mem=512G
#SBATCH --mail-user=rxr456@case.edu
#SBATCH --mail-type=ALL

# Use the Python from your Miniconda installation
export PATH="/home/rxr456/miniconda3/bin:$PATH"
source /home/rxr456/miniconda3/etc/profile.d/conda.sh
conda activate geneformer
python /home/rxr456/geneformer_hpc.py