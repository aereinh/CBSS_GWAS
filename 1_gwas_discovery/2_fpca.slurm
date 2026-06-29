#!/bin/bash
# =============================================================================
# 2_fpca.slurm
#
# Run functional PCA (fPCA) on merged UKB white-matter images using HEIG.
# Computes PC bases retaining all PCs (variance decomposition used downstream
# in step 3 to choose the number of LDRs).
#
# Subjects with QC exclusions (excessive heterozygosity, mismatched sex,
# sex-chromosome aneuploidy) are removed; only white British ancestry kept.
#
# Array index = fiber cluster ID (1-indexed).
# =============================================================================
#SBATCH --job-name=fpca
#SBATCH --output=logs/fpca_%a.out
#SBATCH --error=logs/fpca_%a.err
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --time=1:00:00
#SBATCH --array=1-430%50

# ── User-configurable paths ────────────────────────────────────────────────
MODALITY="dti_FA"
HEIG_DIR="/path/to/heig"
IMAGE_DIR="/path/to/output/images"       # Output of step 1
FPCA_OUT="/path/to/output/fpca"          # Where to write fPCA outputs
QC_DIR="/path/to/qc_data"               # Directory with QC exclusion lists
KEEP_LIST="${QC_DIR}/ukb_white_ancestry.txt"         # Subjects to keep
REMOVE_LIST="${QC_DIR}/excessive_heterozygosity.txt,${QC_DIR}/mismatched_sex.txt,${QC_DIR}/sex_chromosome_aneuploidy.txt"
# ──────────────────────────────────────────────────────────────────────────

module purge
module load anaconda
conda activate heig

mkdir -p "${FPCA_OUT}"

IDX=${SLURM_ARRAY_TASK_ID}

python "${HEIG_DIR}/heig.py" \
    --fpca \
    --out     "${FPCA_OUT}/UKB_${MODALITY}_allphases_${IDX}" \
    --image   "${IMAGE_DIR}/UKB_${MODALITY}_allphases_${IDX}_images.h5" \
    --threads 4 \
    --all-pc \
    --keep    "${KEEP_LIST}" \
    --remove  "${REMOVE_LIST}"
