#!/bin/bash
#SBATCH --job-name=crsom
#SBATCH --output=/path/to/UKB_PSC/logs/multiomics/crossomics.out
#SBATCH --error=/path/to/UKB_PSC/logs/multiomics/crossomics.err
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --time=06:00:00

module load python
cd /path/to/scripts/6_multiomics

python 4b_cross_omics_concordance.py
