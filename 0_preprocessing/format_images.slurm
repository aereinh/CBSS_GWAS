#!/bin/bash
# =============================================================================
# format_images.slurm
#
# Convert the CBSS along-fiber measurements (.mat) into one text file per fiber
# cluster, for every UK Biobank imaging phase.
# Array index = position in PHASES below: Phases 1-6 (discovery) and Phase 7
# (validation) are formatted in the same way.
# =============================================================================
#SBATCH --job-name=format
#SBATCH --output=logs/format_%a.out
#SBATCH --error=logs/format_%a.err
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general
#SBATCH --time=15:00:00
#SBATCH --array=0-6

# ── User-configurable paths ────────────────────────────────────────────────
MODALITY="dti_FA"                                   # dti_FA, dti_MD, dti_L1 (AD)
CBSS_ROOT="/path/to/CBSS_atlas/Modalities"          # CBSS output: <MODALITY>/step2_<MODALITY>_UKB_<PHASE>/
CLUSTER_CSV="/path/to/Fiber_atlas_100points_FullName.csv"
OUT_ROOT="/path/to/UKB_PSC/Modalities"              # formatted text files, same sub-directory layout
SCRIPTS_DIR="/path/to/scripts/0_preprocessing"
# Folder names of the imaging phases; edit to match your data release.
PHASES=(Phase1 Phase2 Phase3 Phase4 Phase5 Phase6 Phase7)
# ──────────────────────────────────────────────────────────────────────────

module purge
module load anaconda
conda activate heig
mkdir -p logs

PHASE=${PHASES[$SLURM_ARRAY_TASK_ID]}
SUBDIR="${MODALITY}/step2_${MODALITY}_UKB_${PHASE}"

python "${SCRIPTS_DIR}/format_images.py" \
    --data-dir    "${CBSS_ROOT}/${SUBDIR}" \
    --subject-ids "${CBSS_ROOT}/${SUBDIR}/UKB_${PHASE}.txt" \
    --cluster-csv "${CLUSTER_CSV}" \
    --out-dir     "${OUT_ROOT}/${SUBDIR}"
