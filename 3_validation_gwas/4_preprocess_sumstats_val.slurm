#!/bin/bash
# =============================================================================
# 4_preprocess_sumstats_val.slurm
#
# Preprocess validation GWAS .parquet output into HEIG summary statistics
# format.  Mirrors 5b_preprocess_sumstats.slurm for the Phase 6 cohort.
#
# Array index = fiber cluster ID (1-indexed).
# =============================================================================
#SBATCH --job-name=sumstats_val
#SBATCH --output=logs/sumstats_val_%a.out
#SBATCH --error=logs/sumstats_val_%a.err
#SBATCH --mem=30G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --time=1:00:00
#SBATCH --array=1-430%50

# ── User-configurable paths ────────────────────────────────────────────────
MODALITY="dti_FA"
HEIG_DIR="/path/to/heig"
GWAS_DIR_VAL="/path/to/output/gwas_val"
SUMSTATS_OUT_VAL="/path/to/output/sumstats_val"
# ──────────────────────────────────────────────────────────────────────────

module purge
module load anaconda
conda activate heig

mkdir -p "${SUMSTATS_OUT_VAL}"

IDX=${SLURM_ARRAY_TASK_ID}

python "${HEIG_DIR}/heig.py" \
    --sumstats \
    --out            "${SUMSTATS_OUT_VAL}/UKB_${MODALITY}_phase6_${IDX}" \
    --ldr-gwas-heig  "${GWAS_DIR_VAL}/UKB_${MODALITY}_phase6_${IDX}.parquet"

if [ $? -eq 0 ]; then
    echo "Preprocessing succeeded for cluster ${IDX}. Removing .parquet folder."
    rm -rf "${GWAS_DIR_VAL}/UKB_${MODALITY}_phase6_${IDX}.parquet/"
else
    echo "Preprocessing FAILED for cluster ${IDX}. Keeping .parquet folder." >&2
fi
