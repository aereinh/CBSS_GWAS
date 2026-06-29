#!/bin/bash
# =============================================================================
# step4a_ld_clump_per_cluster.slurm
#
# Two-stage LD clumping of screened discovery hits within each fiber cluster:
#   Stage 1: broad clumping at r²=0.6 / 1000 kb  →  identifies major blocks
#   Stage 2: fine clumping at r²=0.1 / 500 kb    →  retains one lead SNP
#                                                      per independent locus
#
# The stage-2 output per cluster feeds into step4b (tract-level merging).
#
# Requires PLINK 1.9 and a reference LD panel (e.g., 1000 Genomes European).
# Array index = fiber cluster ID (1-indexed).
# =============================================================================
#SBATCH --job-name=ld_clump
#SBATCH --output=logs/ld_clump_%a.out
#SBATCH --error=logs/ld_clump_%a.err
#SBATCH --time=00:30:00
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=10G
#SBATCH --array=1-430%90

# ── User-configurable paths ────────────────────────────────────────────────
MODALITY="dti_FA"
MAIN_DIR="/path/to/output/output_${MODALITY}/post_screen/disc"
LD_REF="/path/to/reference/g1000_eur"    # PLINK bfile prefix for LD reference panel
# ──────────────────────────────────────────────────────────────────────────

module purge
module add plink/1.90b6.21

ROI=${SLURM_ARRAY_TASK_ID}
FIBER_DIR="${MAIN_DIR}/fibercluster"
INPUT_FILE="${FIBER_DIR}/UKB_${MODALITY}_allphases_${ROI}_sig_vgwas_screen.txt"

if [ ! -f "${INPUT_FILE}" ]; then
    echo "Skipping cluster ${ROI}: no screened hits file found."
    exit 0
fi

echo "LD clumping for cluster ${ROI}"

# Stage 1: broad clumping
plink \
    --bfile    "${LD_REF}" \
    --clump    "${INPUT_FILE}" \
    --clump-kb 1000 \
    --clump-r2 0.6 \
    --out      "${FIBER_DIR}/UKB_${MODALITY}_allphases_${ROI}_clump_r2_0.6_step1"

# Stage 2: fine clumping (input = stage-1 clumped file)
STAGE1="${FIBER_DIR}/UKB_${MODALITY}_allphases_${ROI}_clump_r2_0.6_step1.clumped"

if [ ! -f "${STAGE1}" ]; then
    echo "Stage 1 produced no output for cluster ${ROI}. Skipping stage 2."
    exit 0
fi

plink \
    --bfile    "${LD_REF}" \
    --clump    "${STAGE1}" \
    --clump-kb 500 \
    --clump-r2 0.1 \
    --out      "${FIBER_DIR}/UKB_${MODALITY}_allphases_${ROI}_clump_r2_0.1_step2"

echo "Finished cluster ${ROI}"
