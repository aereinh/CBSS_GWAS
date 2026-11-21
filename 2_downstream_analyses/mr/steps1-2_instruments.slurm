#!/bin/bash
# =============================================================================
# steps1-2_instruments.slurm
#
# Instrument selection for each fiber cluster (array index = cluster ID):
#   1. step1_screen_ldr_sumstats.py  — genome-wide significant SNPs of the LDR GWAS
#   2. step2_ld_clump.sh             — LD clumping to independent loci (r² < 0.01, 10 Mb)
#
# The MR itself (steps 3b-3d and 4) is submitted with submit_mr.sh.
# =============================================================================
#SBATCH --job-name=mr_instruments
#SBATCH --output=logs/mr_instr_%a.out
#SBATCH --error=logs/mr_instr_%a.err
#SBATCH --array=1-430
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --partition=general

# ── User-configurable paths ───────────────────────────────────────────────────
MODALITY="dti_FA"
SUMSTATS_DIR="/path/to/sumstats"               # HEIG LDR sumstats directory (prefix only, no ext)
SCREENED_DIR="/path/to/output/mr/screened"     # Step 1-2 output (SCREEN_DIR in mr_config.sh)
SCRIPTS_DIR="/path/to/scripts/2_downstream_analyses/mr"
# ─────────────────────────────────────────────────────────────────────────────

CLUSTER_ID=${SLURM_ARRAY_TASK_ID}

module purge
module load anaconda
conda activate heig

mkdir -p "${SCREENED_DIR}" logs

echo "Cluster ${CLUSTER_ID} | Job ${SLURM_JOB_ID} | $(date)"

# Step 1: screen LDR sumstats
PREFIX="${SUMSTATS_DIR}/UKB_${MODALITY}_phase16_${CLUSTER_ID}"
OUTPUT="${SCREENED_DIR}/UKB_${MODALITY}_phase16_${CLUSTER_ID}.txt"

python "${SCRIPTS_DIR}/step1_screen_ldr_sumstats.py" \
    --cluster-id "${CLUSTER_ID}" \
    --prefix     "${PREFIX}" \
    --output     "${OUTPUT}"

if [[ ! -f "${OUTPUT}" ]]; then
    echo "No significant SNPs for cluster ${CLUSTER_ID}: no instruments."
    exit 0
fi

# Step 2: LD clump
bash "${SCRIPTS_DIR}/step2_ld_clump.sh" --roi "${CLUSTER_ID}"

echo "Done: cluster ${CLUSTER_ID}"
