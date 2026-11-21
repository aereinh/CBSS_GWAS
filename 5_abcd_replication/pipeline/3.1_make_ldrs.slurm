#!/bin/bash
#SBATCH --job-name=abcd_ldr
#SBATCH --output=/path/to/UKB_PSC/logs/abcd/%x_%A_%a.out
#SBATCH --error=/path/to/UKB_PSC/logs/abcd/%x_%A_%a.err
#SBATCH --mem=20G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --time=00:15:00
#SBATCH --array=1-430%80

# Step 3.1: ABCD LDRs using the UKB phase 1-6 fPCA bases, adjusted for COVAR (config.sh).
# Run with:  bash submit.sh <MOD> ldr
# Starts from the number of LDRs preserving >=95% variance in UKB (trimmed in step 3.2).

: "${MOD:?Set MOD, e.g. bash submit.sh FA ldr}"
source "${PIPE_DIR:?submit with submit.sh}/config.sh"
set_modality "$MOD" || exit 1
load_heig_env
set -eo pipefail

IDX=${SLURM_ARRAY_TASK_ID}
TAG=ABCD_${MOD}_${IDX}
H5=${IMG_DIR}/${TAG}_images.h5
FPCA_LOG=${UKB_FPCA_LOG_PREFIX}_${IDX}.log

[[ -s "${H5}" ]]       || die "missing ${H5} (run the read step first)"
[[ -s "${FPCA_LOG}" ]] || die "missing UKB fPCA log ${FPCA_LOG}"
compgen -G "${UKB_BASES_PREFIX}_${IDX}_bases_*" > /dev/null || die "no UKB bases ${UKB_BASES_PREFIX}_${IDX}_bases_*"

N_LDRS=$(grep -A 2 "The number of LDRs for preserving varying proportions" "${FPCA_LOG}" | tail -n 1 | awk '{print $6}' || true)
[[ -n "${N_LDRS}" ]] || die "could not read the number of LDRs from ${FPCA_LOG}"
echo "MOD=${MOD} IDX=${IDX} N_LDRS(start)=${N_LDRS} COVAR=${COVAR} CAT=${CAT_COVARS}"

mkdir -p "${LDR_DIR}"
python heig.py \
  --make-ldr \
  --out "${LDR_DIR}/${TAG}" \
  --image "${H5}" \
  --n-ldrs "${N_LDRS}" \
  --bases ${UKB_BASES_PREFIX}_${IDX}_bases_* \
  --covar "${COVAR}" \
  --cat-covar-list "${CAT_COVARS}" \
  --threads 4

compgen -G "${LDR_DIR}/${TAG}_ldr_top*.txt" > /dev/null || die "heig --make-ldr did not write ${LDR_DIR}/${TAG}_ldr_top*.txt"
