#!/bin/bash
# =============================================================================
# 1_read_images.slurm
#
# Read white-matter fiber-cluster images of every UK Biobank imaging phase
# with HEIG.  For each fiber cluster:
#   1. Apply outlier removal (values > THRESH set to NA) per phase.
#   2. Read the images of each phase into HEIG's HDF5 format.
#   3. Discovery: merge Phases 1-6 into one image file (tag phase16).
#   4. Validation: read Phase 7 into its own image file (tag phase7).
#   5. Remove the intermediate per-phase discovery files.
#
# Input text files come from 0_preprocessing/format_images.slurm.
# Array index = fiber cluster ID (1-indexed).
# =============================================================================
#SBATCH --job-name=readimgs
#SBATCH --output=logs/readimgs_%a.out
#SBATCH --error=logs/readimgs_%a.err
#SBATCH --mem=15G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --time=00:30:00
#SBATCH --array=1-430

# ── User-configurable paths ────────────────────────────────────────────────
MODALITY="dti_FA"                        # DTI modality (dti_FA, dti_MD, dti_L1)
HEIG_DIR="/path/to/heig"                 # HEIG installation directory
DATA_ROOT="/path/to/UKB_PSC/Modalities"  # Root containing per-modality/per-phase subdirs
COORD_DIR="/path/to/coords"              # Directory with coords_fibersClust_*.txt files
OUT_DIR="/path/to/output/images"         # Discovery images (HEIG HDF5)
OUT_DIR_VAL="/path/to/output/images_val" # Validation images (HEIG HDF5)
SCRIPTS_DIR="/path/to/scripts/0_preprocessing"   # Directory containing remove_outliers.py
TMP_DIR="/path/to/tmp"                   # Temporary directory for cleaned text files
THRESH=100                               # Outlier threshold
# Folder names of the imaging phases; edit to match your data release.
DISC_PHASES=(Phase1 Phase2 Phase3 Phase4 Phase5 Phase6)   # discovery sample
VAL_PHASE="Phase7"                                         # validation sample
# ──────────────────────────────────────────────────────────────────────────

module purge
module load anaconda
conda activate heig

mkdir -p "${OUT_DIR}" "${OUT_DIR_VAL}" "${TMP_DIR}"

IDX=${SLURM_ARRAY_TASK_ID}
TMP_TXT="${TMP_DIR}/fibersClust_${IDX}_clean.txt"

# read_phase <phase folder> <output prefix>
read_phase() {
    local PHASE="$1" OUT_PREFIX="$2"
    local RAW_TXT="${DATA_ROOT}/${MODALITY}/step2_${MODALITY}_UKB_${PHASE}/fibersClust_${IDX}.txt"

    # Remove outliers
    python "${SCRIPTS_DIR}/remove_outliers.py" \
        --in-txt    "${RAW_TXT}" \
        --out-txt   "${TMP_TXT}" \
        --threshold "${THRESH}"

    # Read images into HEIG HDF5 format
    python "${HEIG_DIR}/heig.py" \
        --read-image \
        --out       "${OUT_PREFIX}" \
        --threads   4 \
        --image-txt "${TMP_TXT}" \
        --coord-txt "${COORD_DIR}/coords_fibersClust_${IDX}.txt"

    rm -f "${TMP_TXT}"
}

# ── Discovery: Phases 1-6 ─────────────────────────────────────────────────
IMAGE_LIST=""
for PHASE in "${DISC_PHASES[@]}"; do
    TAG=$(echo "${PHASE}" | tr '[:upper:]' '[:lower:]')
    read_phase "${PHASE}" "${OUT_DIR}/UKB_${MODALITY}_${TAG}_${IDX}"
    IMAGE_LIST+="${OUT_DIR}/UKB_${MODALITY}_${TAG}_${IDX}_images.h5,"
done
IMAGE_LIST="${IMAGE_LIST%,}"   # strip trailing comma

# Merge the discovery phases into a single image file
python "${HEIG_DIR}/heig.py" \
    --read-image \
    --out        "${OUT_DIR}/UKB_${MODALITY}_phase16_${IDX}" \
    --threads    4 \
    --image-list "${IMAGE_LIST}"

# Remove the per-phase HDF5 files (keep only the merged file)
for PHASE in "${DISC_PHASES[@]}"; do
    TAG=$(echo "${PHASE}" | tr '[:upper:]' '[:lower:]')
    rm -f "${OUT_DIR}/UKB_${MODALITY}_${TAG}_${IDX}_images.h5"
    rm -f "${OUT_DIR}/UKB_${MODALITY}_${TAG}_${IDX}.log"
done

# ── Validation: Phase 7 ───────────────────────────────────────────────────
read_phase "${VAL_PHASE}" "${OUT_DIR_VAL}/UKB_${MODALITY}_phase7_${IDX}"
