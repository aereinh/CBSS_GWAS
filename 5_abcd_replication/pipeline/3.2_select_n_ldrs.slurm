#!/bin/bash
#SBATCH --job-name=abcd_select
#SBATCH --output=/path/to/UKB_PSC/logs/abcd/%x_%j.out
#SBATCH --error=/path/to/UKB_PSC/logs/abcd/%x_%j.err
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general
#SBATCH --time=00:30:00

# Step 3.2: write output_ABCD_<MOD>_<COV_TAG>/selLDRs_val.txt
# Run with:  bash submit.sh <MOD> select
# For FA, REUSE_SEL_LDRS (config.sh) copies the original ABCD FA selection instead.

: "${MOD:?Set MOD, e.g. bash submit.sh MD select}"
source "${PIPE_DIR:?submit with submit.sh}/config.sh"
set_modality "$MOD" || exit 1
set -eo pipefail
mkdir -p "${OUT_DIR}"

if [[ -n "${REUSE_SEL_LDRS}" ]]; then
  [[ -s "${REUSE_SEL_LDRS}" ]] || die "REUSE_SEL_LDRS=${REUSE_SEL_LDRS} not found"
  cp -v "${REUSE_SEL_LDRS}" "${SEL_LDRS}"
  exit 0
fi

module purge
module load r
Rscript "${PIPE_DIR}/3.2_select_n_ldrs.R" \
  "$(dirname "${UKB_FPCA_LOG_PREFIX}")" "$(basename "${UKB_FPCA_LOG_PREFIX}")" \
  "${LDR_DIR}" "ABCD_${MOD}" \
  "${SEL_LDRS}"

[[ -s "${SEL_LDRS}" ]] || die "missing ${SEL_LDRS}"
