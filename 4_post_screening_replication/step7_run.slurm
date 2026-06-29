#!/bin/bash
# =============================================================================
# step7_run.slurm
#
# Batch submission for step7_screen_replication_vgwas.py.
# Array index = fiber cluster ID (1-indexed).
# =============================================================================
#SBATCH --job-name=rep_screen
#SBATCH --output=logs/rep_screen_%a.out
#SBATCH --error=logs/rep_screen_%a.err
#SBATCH --time=00:30:00
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --array=1-430%80

# ── User-configurable ─────────────────────────────────────────────────────
MODALITY="dti_FA"
NULL_TAG="1e4"      # Match the threshold used in step 6 and the cluster-val subdir
MAIN_DIR="/path/to/output"
SCRIPTS_DIR="/path/to/scripts/4_post_screening_replication"
# ──────────────────────────────────────────────────────────────────────────

module purge
module load anaconda
conda activate heig

python "${SCRIPTS_DIR}/step7_screen_replication_vgwas.py" \
    --modality  "${MODALITY}" \
    --cluster   "${SLURM_ARRAY_TASK_ID}" \
    --null-tag  "${NULL_TAG}" \
    --main-dir  "${MAIN_DIR}"
