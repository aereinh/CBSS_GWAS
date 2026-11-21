#!/bin/bash
#SBATCH --job-name=abcd_read
#SBATCH --output=/path/to/UKB_PSC/logs/abcd/%x_%A_%a.out
#SBATCH --error=/path/to/UKB_PSC/logs/abcd/%x_%A_%a.err
#SBATCH --mem=15G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --time=00:30:00
#SBATCH --array=1-430%80

# Step 1: fibersClust_<IDX>.txt -> output_ABCD_<MOD>/images/ABCD_<MOD>_<IDX>_images.h5
# Run with:  bash submit.sh <MOD> read
# Not needed for FA (images already exist and don't depend on covariates).

: "${MOD:?Set MOD, e.g. bash submit.sh MD read}"
source "${PIPE_DIR:?submit with submit.sh}/config.sh"
set_modality "$MOD" || exit 1
load_heig_env
set -eo pipefail

IDX=${SLURM_ARRAY_TASK_ID}
TAG=ABCD_${MOD}_${IDX}
IMG_TXT=${IMG_TXT_DIR}/fibersClust_${IDX}.txt
H5=${IMG_DIR}/${TAG}_images.h5

if [[ -s "${H5}" && -z "${OVERWRITE:-}" ]]; then
  echo "[SKIP] ${H5} exists (set OVERWRITE=1 to redo)"
  exit 0
fi
[[ -s "${IMG_TXT}" ]] || die "missing ${IMG_TXT} (run the format step first)"
mkdir -p "${IMG_DIR}"

# Optional outlier cleaning (|value| > OUTLIER_THRESH -> NA), as in the UKB dti_* runs
INPUT_TXT=${IMG_TXT}
if [[ -n "${OUTLIER_THRESH}" ]]; then
  INPUT_TXT=${TMP_ROOT}/${TAG}_clean_${SLURM_JOB_ID}.txt
  trap 'rm -f "${INPUT_TXT}"' EXIT
  python "${PIPE_DIR}/0.1_remove_outliers.py" \
    --in-txt "${IMG_TXT}" --out-txt "${INPUT_TXT}" --threshold "${OUTLIER_THRESH}"
fi

python heig.py \
  --read-image \
  --out "${IMG_DIR}/${TAG}" \
  --threads 4 \
  --image-txt "${INPUT_TXT}" \
  --coord-txt "${COORD_DIR}/coords_fibersClust_${IDX}.txt"

[[ -s "${H5}" ]] || die "heig --read-image did not write ${H5}"

# HEIG truncates the long NDAR_INV IDs; restore them from the txt file
python "${PIPE_DIR}/fix_h5_ids.py" --h5 "${H5}" --ids-txt "${IMG_TXT}"
