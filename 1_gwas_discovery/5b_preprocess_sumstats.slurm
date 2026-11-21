#!/bin/bash
# =============================================================================
# 5b_preprocess_sumstats.slurm
#
# Convert HEIG GWAS output (.parquet) to preprocessed summary statistics
# format required for voxel-GWAS and downstream analyses.
#
# On success, the raw .parquet folder is deleted to save disk space.
#
# Array index = fiber cluster ID (1-indexed).
# =============================================================================
#SBATCH --job-name=sumstats
#SBATCH --output=logs/sumstats_%a.out
#SBATCH --error=logs/sumstats_%a.err
#SBATCH --mem=30G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --time=1:00:00
#SBATCH --array=1-430%50

# ── User-configurable paths ────────────────────────────────────────────────
MODALITY="dti_FA"
HEIG_DIR="/path/to/heig"
GWAS_DIR="/path/to/output/gwas"          # Input: raw .parquet GWAS output
SUMSTATS_OUT="/path/to/output/sumstats"  # Output: preprocessed summary stats
# ──────────────────────────────────────────────────────────────────────────

module purge
module load anaconda
conda activate heig

mkdir -p "${SUMSTATS_OUT}"

IDX=${SLURM_ARRAY_TASK_ID}

python "${HEIG_DIR}/heig.py" \
    --sumstats \
    --out            "${SUMSTATS_OUT}/UKB_${MODALITY}_phase16_${IDX}" \
    --ldr-gwas-heig  "${GWAS_DIR}/UKB_${MODALITY}_phase16_${IDX}.parquet"

# Delete raw .parquet only if preprocessing succeeded
if [ $? -eq 0 ]; then
    echo "Preprocessing succeeded for cluster ${IDX}. Removing .parquet folder."
    rm -rf "${GWAS_DIR}/UKB_${MODALITY}_phase16_${IDX}.parquet/"
else
    echo "Preprocessing FAILED for cluster ${IDX}. Keeping .parquet folder." >&2
fi
