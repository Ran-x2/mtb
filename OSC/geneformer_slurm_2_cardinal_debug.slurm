#!/bin/bash
#SBATCH --account=PDS0353
#SBATCH --job-name=rxr456_geneformer_mtb_2_cardinal_debug
#SBATCH --output=geneformer_mtb_output_cardinal_debug.txt
#SBATCH --error=geneformer_mtb_error_cardinal_debug.txt
#SBATCH --time=20:00
#SBATCH --partition=debug
#SBATCH --gres=gpu:1
#SBATCH --mem=128gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=80
#SBATCH --mail-user=rxr456@case.edu
#SBATCH --mail-type=ALL

# Use the Python from your Miniconda installation

export PATH="/users/PDS0353/rxr456/miniconda3/bin:$PATH"
source /users/PDS0353/rxr456/miniconda3/etc/profile.d/conda.sh
storage_dir=/users/PDS0353/rxr456/mtb_250409

conda activate geneformer312
echo "ready for python!"
python /users/PDS0353/rxr456/mtb_250409/geneformer_hpc_2.py