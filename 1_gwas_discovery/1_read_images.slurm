#!/bin/bash
# =============================================================================
# 1_read_images.slurm
#
# Read white-matter fiber-cluster images across UKB release phases using HEIG.
# For each cluster ROI:
#   1. Apply outlier removal (values > THRESH set to NA) per phase.
#   2. Read images into HEIG's HDF5 format for each phase.
#   3. Merge all phases into a single combined image file.
#   4. Clean up intermediate per-phase files.
#
# Array index = fiber cluster ID (1-indexed).
# Adjust MODALITY and array range as needed.
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
MODALITY="dti_FA"                        # DTI modality (e.g., dti_FA, dti_MD, dti_RD)
HEIG_DIR="/path/to/heig"                 # HEIG installation directory
DATA_ROOT="/path/to/UKB_PSC/Modalities"  # Root containing per-modality/per-phase subdirs
COORD_DIR="/path/to/coords"              # Directory with coords_fibersClust_*.txt files
OUT_DIR="/path/to/output/images"         # Output directory for HEIG image HDF5 files
SCRIPTS_DIR="/path/to/scripts"           # Directory containing remove_outliers.py
TMP_DIR="/path/to/tmp"                   # Temporary directory for cleaned text files
THRESH=100                               # Outlier threshold
# ──────────────────────────────────────────────────────────────────────────

module purge
module load anaconda
conda activate heig

mkdir -p "${OUT_DIR}" "${TMP_DIR}"

IDX=${SLURM_ARRAY_TASK_ID}
TMP_TXT="${TMP_DIR}/fibersClust_${IDX}_clean.txt"

# Read and merge each UKB release phase
PHASES=(Phase12 Phase3 Phase4 Phase5 Phase6)
PHASE_TAGS=(phase12 phase3 phase4 phase5 phase6)
IMAGE_LIST=""

for i in "${!PHASES[@]}"; do
    PHASE="${PHASES[$i]}"
    TAG="${PHASE_TAGS[$i]}"
    RAW_TXT="${DATA_ROOT}/${MODALITY}/step2_${MODALITY}_UKB_${PHASE}/fibersClust_${IDX}.txt"

    # Remove outliers
    python "${SCRIPTS_DIR}/remove_outliers.py" \
        --in-txt    "${RAW_TXT}" \
        --out-txt   "${TMP_TXT}" \
        --threshold "${THRESH}"

    # Read images into HEIG HDF5 format
    python "${HEIG_DIR}/heig.py" \
        --read-image \
        --out     "${OUT_DIR}/UKB_${MODALITY}_${TAG}_${IDX}" \
        --threads 4 \
        --image-txt "${TMP_TXT}" \
        --coord-txt "${COORD_DIR}/coords_fibersClust_${IDX}.txt"

    rm -f "${TMP_TXT}"

    IMAGE_LIST+="${OUT_DIR}/UKB_${MODALITY}_${TAG}_${IDX}_images.h5,"
done

# Merge all phases into a single combined image file
IMAGE_LIST="${IMAGE_LIST%,}"   # strip trailing comma

python "${HEIG_DIR}/heig.py" \
    --read-image \
    --out        "${OUT_DIR}/UKB_${MODALITY}_allphases_${IDX}" \
    --threads    4 \
    --image-list "${IMAGE_LIST}"

# Remove individual per-phase HDF5 files (keep only the merged file)
for TAG in "${PHASE_TAGS[@]}"; do
    rm -f "${OUT_DIR}/UKB_${MODALITY}_${TAG}_${IDX}_images.h5"
    rm -f "${OUT_DIR}/UKB_${MODALITY}_${TAG}_${IDX}.log"
done
