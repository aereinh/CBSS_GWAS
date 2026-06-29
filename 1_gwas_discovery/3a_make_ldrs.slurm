#!/bin/bash
# =============================================================================
# 3a_make_ldrs.slurm
#
# Compute Low-Dimensional Representations (LDRs) for each fiber cluster using
# HEIG.  The number of LDRs is selected from the fPCA log file to preserve
# >=95% of image variance; this is used as the upper bound for step 3b.
#
# Covariates: age, sex (categorical), interaction terms, 40 genetic PCs,
# assessment centre (categorical).
#
# Array index = fiber cluster ID (1-indexed).
# =============================================================================
#SBATCH --job-name=make_ldrs
#SBATCH --output=logs/ldrs_%a.out
#SBATCH --error=logs/ldrs_%a.err
#SBATCH --mem=30G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --time=00:30:00
#SBATCH --array=1-430%50

# ── User-configurable paths ────────────────────────────────────────────────
MODALITY="dti_FA"
HEIG_DIR="/path/to/heig"
IMAGE_DIR="/path/to/output/images"
FPCA_DIR="/path/to/output/fpca"
LDR_OUT="/path/to/output/ldrs"
COVAR_FILE="/path/to/pheno_covar/covariates.txt"
CAT_COVAR="sex,assess_center"             # Comma-separated categorical covariates
# ──────────────────────────────────────────────────────────────────────────

module purge
module load anaconda
conda activate heig

mkdir -p "${LDR_OUT}"

IDX=${SLURM_ARRAY_TASK_ID}

# Read the recommended number of LDRs from the fPCA log (>=95% variance)
N_LDRS=$(grep -A 2 "The number of LDRs for preserving varying proportions" \
    "${FPCA_DIR}/UKB_${MODALITY}_allphases_${IDX}.log" \
    | tail -n 1 | awk '{print $6}')

python "${HEIG_DIR}/heig.py" \
    --make-ldr \
    --out          "${LDR_OUT}/UKB_${MODALITY}_allphases_${IDX}" \
    --image        "${IMAGE_DIR}/UKB_${MODALITY}_allphases_${IDX}_images.h5" \
    --n-ldrs       "${N_LDRS}" \
    --bases        "${FPCA_DIR}/UKB_${MODALITY}_allphases_${IDX}_bases_"* \
    --covar        "${COVAR_FILE}" \
    --cat-covar-list "${CAT_COVAR}" \
    --threads      4
