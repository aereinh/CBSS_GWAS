#!/bin/bash
#SBATCH --job-name=mo
#SBATCH --output=/path/to/UKB_PSC/logs/multiomics/voxatl_agg.out
#SBATCH --error=/path/to/UKB_PSC/logs/multiomics/voxatl_agg.err
#SBATCH --mem=21G
#SBATCH --cpus-per-task=1
#SBATCH --time=00:45:00

cd /path/to/scripts/6_multiomics

module load python
python 4a_aggregate_selected_features.py

