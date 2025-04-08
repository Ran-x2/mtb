#!/bin/bash
#SBATCH --job-name=rxr456_geneformer_mtb_3
#SBATCH --output=geneformer_mtb_output.txt
#SBATCH --error=geneformer_mtb_error.txt
#SBATCH --time=2:00:00
#SBATCH --partition=gpu
#SBATCH --gres=gpu:1
#SBATCH --mem=768G
#SBATCH --mail-user=rxr456@case.edu
#SBATCH --mail-type=ALL

# Use the Python from your Miniconda installation
export PATH="/home/rxr456/miniconda3/bin:$PATH"
export PATH=/home/rxr456/miniconda3/envs/geneformer/bin:$PATH
source /home/rxr456/miniconda3/etc/profile.d/conda.sh
export PYTHONPATH=/home/rxr456/miniconda3/envs/geneformer/lib/python3.10/site-packages:$PYTHONPATH

conda activate geneformer
echo "ready for python!"
python /home/rxr456/geneformer_hpc_3.py