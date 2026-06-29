#!/bin/bash
# =============================================================================
# heritability.slurm
#
# Estimate SNP-based heritability of white-matter voxel phenotypes using HEIG.
# Uses pre-computed LD matrices (HapMap3 SNPs, regularized).
#
# Array index = fiber cluster ID (1-indexed).
# =============================================================================
#SBATCH --job-name=heritability
#SBATCH --output=logs/heri_%a.out
#SBATCH --error=logs/heri_%a.err
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --time=00:10:00
#SBATCH --array=1-430%80

# ── User-configurable paths ────────────────────────────────────────────────
MODALITY="dti_FA"
HEIG_DIR="/path/to/heig"
SUMSTATS_DIR="/path/to/output/sumstats"
FPCA_DIR="/path/to/output/fpca"
LDR_DIR="/path/to/output/ldrs"
HERI_OUT="/path/to/output/heritability"
NLDR_FILE="/path/to/output/selLDRs.txt"

# LD matrices (HapMap3 SNPs, per-chromosome)
# {1:22} is expanded by HEIG to iterate over chromosomes
LD_INV="/path/to/ldmatrix/ukb_ld_inv_chr{1:22}"
LD="/path/to/ldmatrix/ukb_ld_chr{1:22}"
# ──────────────────────────────────────────────────────────────────────────

module purge
module load anaconda
conda activate heig

mkdir -p "${HERI_OUT}"

IDX=${SLURM_ARRAY_TASK_ID}
N_LDRS=$(awk -F',' -v i=${IDX} 'NR == i+1 {print $2}' "${NLDR_FILE}")

python "${HEIG_DIR}/heig.py" \
    --heri-gc \
    --out          "${HERI_OUT}/UKB_${MODALITY}_allphases_${IDX}" \
    --n-ldrs       "${N_LDRS}" \
    --ldr-sumstats "${SUMSTATS_DIR}/UKB_${MODALITY}_allphases_${IDX}" \
    --bases        "${FPCA_DIR}/UKB_${MODALITY}_allphases_${IDX}_bases_top"* \
    --ldr-cov      "${LDR_DIR}/UKB_${MODALITY}_allphases_${IDX}_ldr_cov_top"* \
    --threads      4 \
    --ld-inv       "${LD_INV}" \
    --ld           "${LD}" \
    --heri-only
