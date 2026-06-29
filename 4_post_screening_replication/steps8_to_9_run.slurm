#!/bin/bash
# =============================================================================
# steps8_to_9_run.slurm
#
# Run genetic replication (step 8) and spatial overlap (step 9) sequentially
# for each white-matter tract.  Array index maps to tract index (0-based).
#
# Edit TRACTS list to match your tract set.
# =============================================================================
#SBATCH --job-name=rep_spatial
#SBATCH --output=logs/rep_spatial_%a.out
#SBATCH --error=logs/rep_spatial_%a.err
#SBATCH --time=01:00:00
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=25G
#SBATCH --array=0-16

# ── User-configurable ─────────────────────────────────────────────────────
MODALITY="dti_FA"
NULL_TAG="1e4"
MAIN_DIR="/path/to/output"
SCRIPTS_DIR="/path/to/scripts/4_post_screening_replication"
# ──────────────────────────────────────────────────────────────────────────

TRACTS=(
    ACR
    ALIC
    BCC
    CGH
    CP
    EC
    FX
    GCC
    IFO
    PCR
    PLIC
    PTR
    RLIC
    SCC
    SCR
    SLF
    SS
)

TRACT="${TRACTS[$SLURM_ARRAY_TASK_ID]}"

echo "Tract: ${TRACT}  |  Job: ${SLURM_JOB_ID}  |  Array: ${SLURM_ARRAY_TASK_ID}"

module purge
module load anaconda
conda activate heig

cd "${SCRIPTS_DIR}"

# Step 8: genetic replication
echo "--- Step 8: Genetic replication ---"
python step8_replication_one_tract.py \
    --modality "${MODALITY}" \
    --tract    "${TRACT}" \
    --null-tag "${NULL_TAG}" \
    --main-dir "${MAIN_DIR}"

# Step 9: spatial overlap
echo "--- Step 9: Spatial overlap ---"
python step9_spatial_overlap_one_tract.py \
    --modality "${MODALITY}" \
    --tract    "${TRACT}" \
    --null-tag "${NULL_TAG}" \
    --main-dir "${MAIN_DIR}"

echo "Done: ${TRACT}"
