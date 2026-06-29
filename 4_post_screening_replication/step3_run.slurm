#!/bin/bash
# =============================================================================
# step3_run.slurm
#
# Batch submission for step3_screen_snps.py across all fiber clusters.
# Array index = fiber cluster ID (1-indexed).
# =============================================================================
#SBATCH --job-name=screen_snps
#SBATCH --output=logs/screen_snps_%a.out
#SBATCH --error=logs/screen_snps_%a.err
#SBATCH --time=00:10:00
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --array=1-430%70

# ── User-configurable ─────────────────────────────────────────────────────
MODALITY="dti_FA"
MAIN_DIR="/path/to/output"
CLUSTER_DIR="/path/to/output/cluster"
SCRIPTS_DIR="/path/to/scripts/4_post_screening_replication"
# ──────────────────────────────────────────────────────────────────────────

module purge
module load anaconda
conda activate heig

python "${SCRIPTS_DIR}/step3_screen_snps.py" \
    --modality    "${MODALITY}" \
    --roi         "${SLURM_ARRAY_TASK_ID}" \
    --main-dir    "${MAIN_DIR}" \
    --cluster-dir "${CLUSTER_DIR}"
