#!/bin/bash
# =============================================================================
# 1_make_ldrs_val.slurm
#
# Compute LDRs for the validation sample (UKB Phase 7).
# Uses the fPCA bases estimated from the discovery sample (all phases),
# ensuring LDRs are projected onto the same functional basis.
#
# Run 3b_select_num_ldrs.R (with PHASE="val") after this step to determine
# the optimal number of LDRs for the validation cohort.
#
# Array index = fiber cluster ID (1-indexed).
# =============================================================================
#SBATCH --job-name=ldrs_val
#SBATCH --output=logs/ldrs_val_%a.out
#SBATCH --error=logs/ldrs_val_%a.err
#SBATCH --mem=30G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --time=00:30:00
#SBATCH --array=1-430%50

# ── User-configurable paths ────────────────────────────────────────────────
MODALITY="dti_FA"
HEIG_DIR="/path/to/heig"
IMAGE_DIR_VAL="/path/to/output/images_val"   # Phase 7 images (from 1_gwas_discovery/1_read_images.slurm)
FPCA_DIR="/path/to/output/fpca"              # Discovery fPCA bases
LDR_OUT_VAL="/path/to/output/ldrs_val"
COVAR_FILE="/path/to/pheno_covar/covariates.txt"
CAT_COVAR="sex,assess_center"
NLDR_FILE="/path/to/output/selLDRs.txt"     # From discovery (upper bound)
# ──────────────────────────────────────────────────────────────────────────

module purge
module load anaconda
conda activate heig

mkdir -p "${LDR_OUT_VAL}"

IDX=${SLURM_ARRAY_TASK_ID}
N_LDRS=$(awk -F',' -v i=${IDX} 'NR == i+1 {print $2}' "${NLDR_FILE}")

python "${HEIG_DIR}/heig.py" \
    --make-ldr \
    --out          "${LDR_OUT_VAL}/UKB_${MODALITY}_phase7_${IDX}" \
    --image        "${IMAGE_DIR_VAL}/UKB_${MODALITY}_phase7_${IDX}_images.h5" \
    --n-ldrs       "${N_LDRS}" \
    --bases        "${FPCA_DIR}/UKB_${MODALITY}_phase16_${IDX}_bases_top"* \
    --covar        "${COVAR_FILE}" \
    --cat-covar-list "${CAT_COVAR}" \
    --threads      4
