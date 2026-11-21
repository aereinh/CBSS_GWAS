#!/bin/bash
#SBATCH --job-name=merge_abcd
#SBATCH --output=/path/to/UKB_PSC/logs/merge_abcd_%j.out
#SBATCH --error=/path/to/UKB_PSC/logs/merge_abcd_%j.err
#SBATCH --mem=35G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --time=01:00:00

echo "Starting ABCD chromosome merge"
echo "Job ID: $SLURM_JOB_ID"
echo "Running on: $(hostname)"
echo "Start time: $(date)"

module purge
module load plink/2.00

BASE="/path/to/UKB_PSC/output/genotype/ABCD"
MERGELIST="${BASE}/merge_list_abcd.txt"
OUTPREFIX="${BASE}/ABCD_201_all_chr"

plink2 \
  --bfile ${BASE}/ABCD_201_chr1_imputed \
  --pmerge-list ${MERGELIST} bfile \
  --make-bed \
  --threads ${SLURM_CPUS_PER_TASK} \
  --out ${OUTPREFIX}

echo "Finished merge"
echo "End time: $(date)"