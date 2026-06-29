#!/bin/bash
# =============================================================================
# step1_run.slurm
#
# Batch submission for step1_get_unique_snps.py across all fiber clusters.
# Array index = fiber cluster ID (1-indexed).
# =============================================================================
#SBATCH --job-name=unique_snps
#SBATCH --output=logs/unique_snps_%a.out
#SBATCH --error=logs/unique_snps_%a.err
#SBATCH --time=00:30:00
#SBATCH --mem=30G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general
#SBATCH --array=1-430%80

# ── User-configurable ─────────────────────────────────────────────────────
MODALITY="dti_FA"
MAIN_DIR="/path/to/output"
SCRIPTS_DIR="/path/to/scripts/4_post_screening_replication"
# ──────────────────────────────────────────────────────────────────────────

module purge
module load anaconda
conda activate heig

python "${SCRIPTS_DIR}/step1_get_unique_snps.py" \
    --modality "${MODALITY}" \
    --roi      "${SLURM_ARRAY_TASK_ID}" \
    --main-dir "${MAIN_DIR}"
