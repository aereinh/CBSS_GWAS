#!/bin/bash
#SBATCH --job-name=abcd_format
#SBATCH --output=/path/to/UKB_PSC/logs/abcd/%x_%A_%a.out
#SBATCH --error=/path/to/UKB_PSC/logs/abcd/%x_%A_%a.err
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general
#SBATCH --time=00:30:00
#SBATCH --array=1-430%90

# Step 0: ABCD <IDX>_Aligned.mat files -> Modalities/ABCD_<MOD>/fibersClust_<IDX>.txt
# Run with:  bash submit.sh <MOD> format
# Not needed for FA (the txt files already exist).

: "${MOD:?Set MOD, e.g. bash submit.sh MD format}"
source "${PIPE_DIR:?submit with submit.sh}/config.sh"
set_modality "$MOD" || exit 1
set -eo pipefail

IDX=${SLURM_ARRAY_TASK_ID}
OUT_TXT=${IMG_TXT_DIR}/fibersClust_${IDX}.txt

if [[ -s "${OUT_TXT}" && -z "${OVERWRITE:-}" ]]; then
  echo "[SKIP] ${OUT_TXT} exists (set OVERWRITE=1 to redo)"
  exit 0
fi

module purge
module load python/3.9.6
mkdir -p "${IMG_TXT_DIR}"

python "${PIPE_DIR}/0_format_images.py" \
  --cluster "${IDX}" \
  --data-dir "${RAW_DIR}" \
  --out-dir "${IMG_TXT_DIR}" \
  --clust-file "${CLUST_CSV}"

[[ -s "${OUT_TXT}" ]] || die "missing output ${OUT_TXT}"
