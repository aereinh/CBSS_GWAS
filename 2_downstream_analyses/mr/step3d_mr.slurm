#!/bin/bash
#SBATCH --job-name=mr
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general
#SBATCH --time=06:00:00

# Step 3d: one array task = one fiber cluster, for the trait given in TRAIT.
# Submitted by submit_mr.sh run (one array per trait).
: "${TRAIT:?TRAIT not set (use submit_mr.sh run)}"
source "${MR_DIR:?submit with submit_mr.sh}/mr_config.sh"
set -eo pipefail
FID=${SLURM_ARRAY_TASK_ID}
OUT=${OUT_ROOT}/${TRAIT}/${EXP_PREFIX}${FID}_${TRAIT}_MR_results.tsv
DONE=${OUT_ROOT}/${TRAIT}/${EXP_PREFIX}${FID}_${TRAIT}_MR.done

if [[ -z "${OVERWRITE:-}" && ( -s "${OUT}" || -e "${DONE}" ) ]]; then
  echo "[SKIP] ${TRAIT} fiber ${FID} already done (OVERWRITE=1 to redo)"
  exit 0
fi
load_r
Rscript "${MR_DIR}/step3d_run_mr.R" "${FID}" "${TRAIT}"
