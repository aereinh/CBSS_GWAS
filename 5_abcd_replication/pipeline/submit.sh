#!/bin/bash
# =====================================================================
# Submit one step of the ABCD pipeline for one modality.
#
#   bash submit.sh <MOD> <STEP> [--array SPEC] [--after JOBID]
#
#   MOD    FA | MD | AD
#   STEP   format  -> 0_format_images.slurm        (.mat -> fibersClust_<IDX>.txt)
#          read    -> 1_read_images.slurm          (txt -> images.h5)
#          ldr     -> 3.1_make_ldrs.slurm          (LDRs with UKB bases + covariates)
#          select  -> 3.2_select_n_ldrs.slurm      (selLDRs_val.txt)
#          gwas    -> 4-5_rel_gwas_sumstats.slurm  (relatedness + GWAS + sumstats)
#   --array  override the array range, e.g. --array 7,19,200-210
#   --after  start only after job JOBID finished successfully (afterok)
#
# Prints the job ID, so steps can be chained:
#   J=$(bash submit.sh MD read); bash submit.sh MD ldr --after $J
# =====================================================================
set -euo pipefail
PIPE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export PIPE_DIR
source "${PIPE_DIR}/config.sh"

[[ $# -ge 2 ]] || { sed -n '3,20p' "$0"; exit 1; }
MOD_ARG=$1; STEP=$2; shift 2
ARRAY=""; AFTER=""
while [[ $# -gt 0 ]]; do
  case "$1" in
    --array) ARRAY="$2"; shift 2 ;;
    --after) AFTER="$2"; shift 2 ;;
    *) echo "Unknown option $1" >&2; exit 1 ;;
  esac
done

set_modality "$MOD_ARG"

case "$STEP" in
  format) SCRIPT=0_format_images.slurm ;;
  read)   SCRIPT=1_read_images.slurm ;;
  ldr)    SCRIPT=3.1_make_ldrs.slurm ;;
  select) SCRIPT=3.2_select_n_ldrs.slurm ;;
  gwas)   SCRIPT=4-5_rel_gwas_sumstats.slurm ;;
  *) echo "Unknown step '$STEP' (format|read|ldr|select|gwas)" >&2; exit 1 ;;
esac

LOGDIR=${LOG_ROOT}/${RUN_NAME}/${STEP}
mkdir -p "${LOGDIR}"

if [[ "$STEP" == "select" ]]; then
  LOGPAT="${LOGDIR}/%x_%j"
else
  LOGPAT="${LOGDIR}/%x_%a"
fi

ARGS=(--parsable
      --job-name="${STEP}_${RUN_NAME}"
      --output="${LOGPAT}.out"
      --error="${LOGPAT}.err"
      --export=ALL,MOD="${MOD}",PIPE_DIR="${PIPE_DIR}")
[[ -n "$ARRAY" ]] && ARGS+=(--array="$ARRAY")
[[ -n "$AFTER" ]] && ARGS+=(--dependency=afterok:"$AFTER" --kill-on-invalid-dep=yes)

JOBID=$(sbatch "${ARGS[@]}" "${PIPE_DIR}/${SCRIPT}")
echo "Submitted ${STEP} for ${MOD} (${RUN_NAME}): job ${JOBID}  logs: ${LOGDIR}" >&2
echo "${JOBID}"
