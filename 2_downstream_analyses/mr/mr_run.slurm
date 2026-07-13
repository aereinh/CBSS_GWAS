#!/bin/bash
# =============================================================================
# mr_run.slurm
#
# SLURM array job to run the full MR pipeline (steps 1–3) for each fiber
# cluster.  Each task runs:
#   1. step1_screen_ldr_sumstats.py  — extract genome-wide significant SNPs
#   2. step2_ld_clump.sh             — LD clumping to get independent loci
#   3. step3_run_mr.R                — two-sample MR against 14 traits
#
# After all tasks complete, run step4_agg_results.R locally to combine outputs.
#
# Set OPENGWAS_JWT as an environment variable before submission:
#   export OPENGWAS_JWT="<your_token>"
#   sbatch mr_run.slurm
# =============================================================================
#SBATCH --job-name=mr_pipeline
#SBATCH --output=logs/mr_%a.out
#SBATCH --error=logs/mr_%a.err
#SBATCH --array=1-430
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=16G
#SBATCH --time=00:30:00
#SBATCH --partition=general

# ── User-configurable paths ───────────────────────────────────────────────────
MODALITY="dti_FA"
SUMSTATS_DIR="/path/to/sumstats"               # HEIG LDR sumstats directory (prefix only, no ext)
SCREENED_DIR="/path/to/output/mr/screened"     # Step 1 output
SCRIPTS_DIR="/path/to/scripts/2_downstream_analyses/mr"
# ─────────────────────────────────────────────────────────────────────────────

CLUSTER_ID=${SLURM_ARRAY_TASK_ID}

module purge
module load anaconda
conda activate heig
module load r/4.5.0
export R_LIBS_USER="${HOME}/R/x86_64-pc-linux-gnu-library/4.3"

mkdir -p "${SCREENED_DIR}" logs

echo "Cluster ${CLUSTER_ID} | Job ${SLURM_JOB_ID} | $(date)"

# Step 1: screen LDR sumstats
PREFIX="${SUMSTATS_DIR}/UKB_${MODALITY}_allphases_${CLUSTER_ID}"
OUTPUT="${SCREENED_DIR}/UKB_${MODALITY}_allphases_${CLUSTER_ID}.txt"

python "${SCRIPTS_DIR}/step1_screen_ldr_sumstats.py" \
    --cluster-id "${CLUSTER_ID}" \
    --prefix     "${PREFIX}" \
    --output     "${OUTPUT}"

if [[ ! -f "${OUTPUT}" ]]; then
    echo "No significant SNPs for cluster ${CLUSTER_ID}. Skipping MR."
    exit 0
fi

# Step 2: LD clump
bash "${SCRIPTS_DIR}/step2_ld_clump.sh" --roi "${CLUSTER_ID}"

# Step 3: MR
Rscript "${SCRIPTS_DIR}/step3_run_mr.R" "${CLUSTER_ID}"

echo "Done: cluster ${CLUSTER_ID}"
