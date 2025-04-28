#!/bin/bash
#SBATCH --job-name=rxr456_geneformer_mtb_2_2states
#SBATCH --output=geneformer_mtb_output_2_2states.txt
#SBATCH --error=geneformer_mtb_error_2_2states.txt
#SBATCH --time=7-00:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=128gb
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=65
#SBATCH --mail-user=rxr456@case.edu
#SBATCH --mail-type=ALL

# Use the Python from your Miniconda installation
export PATH="/home/rxr456/miniconda3/bin:$PATH"
export PATH=/home/rxr456/miniconda3/envs/geneformer/bin:$PATH
source /home/rxr456/miniconda3/etc/profile.d/conda.sh
export PYTHONPATH=/home/rxr456/miniconda3/envs/geneformer/lib/python3.10/site-packages:$PYTHONPATH

conda activate geneformer
echo "ready for python!"
python /mnt/vstor/SOM_PATH_DKB50/members/rxr456/mtb_250409/geneformer_hpc_2_2states.py