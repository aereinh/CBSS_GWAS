#!/bin/bash
# =============================================================================
# 5c_voxelgwas_screening.slurm
#
# Run voxel-level GWAS (VGWAS) across all voxels for each fiber cluster using
# HEIG.  Outputs SNP-voxel association statistics for variants passing a
# genome-wide significance threshold (adjusted for the effective number of
# independent tests across voxels).
#
# SIG_THRESH = 5e-8 / (effective number of independent voxels)
# The effective number should be computed separately (e.g., via HEIG's
# effective-number-of-tests utility) and set below.
#
# Array index = fiber cluster ID (1-indexed).
# =============================================================================
#SBATCH --job-name=voxgwas
#SBATCH --output=logs/voxelgwas_%a.out
#SBATCH --error=logs/voxelgwas_%a.err
#SBATCH --mem=35G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --time=3:00:00
#SBATCH --array=1-430%50

# ── User-configurable paths ────────────────────────────────────────────────
MODALITY="dti_FA"
HEIG_DIR="/path/to/heig"
SUMSTATS_DIR="/path/to/output/sumstats"
FPCA_DIR="/path/to/output/fpca"
LDR_DIR="/path/to/output/ldrs"
VGWAS_OUT="/path/to/output/voxelgwas"
NLDR_FILE="/path/to/output/selLDRs.txt"

# Significance threshold = 5e-8 / effective number of independent voxels
# Adjust this value based on your data's effective test count
SIG_THRESH=9.124955e-12    # Example: 5e-8 / 5479 effective voxels
# ──────────────────────────────────────────────────────────────────────────

module purge
module load anaconda
conda activate heig

mkdir -p "${VGWAS_OUT}"

IDX=${SLURM_ARRAY_TASK_ID}
N_LDRS=$(awk -F',' -v i=${IDX} 'NR == i+1 {print $2}' "${NLDR_FILE}")

python "${HEIG_DIR}/heig.py" \
    --voxel-gwas \
    --out          "${VGWAS_OUT}/UKB_${MODALITY}_phase16_${IDX}_all_voxels_whole_genome" \
    --n-ldrs       "${N_LDRS}" \
    --ldr-sumstats "${SUMSTATS_DIR}/UKB_${MODALITY}_phase16_${IDX}" \
    --bases        "${FPCA_DIR}/UKB_${MODALITY}_phase16_${IDX}_bases_top"* \
    --ldr-cov      "${LDR_DIR}/UKB_${MODALITY}_phase16_${IDX}_ldr_cov_top"* \
    --threads      4 \
    --sig-thresh   "${SIG_THRESH}"
