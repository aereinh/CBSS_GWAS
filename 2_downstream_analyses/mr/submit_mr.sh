#!/bin/bash
# =====================================================================
# Submit the MR analysis (steps 3b-3d and 4).
#
#   bash submit_mr.sh prep [trait ...]     step 3b (instrument list, once) + step 3c (harmonize traits)
#   bash submit_mr.sh harmonize trait ...  step 3c only, for the traits named (after editing trait_config.R)
#   bash submit_mr.sh run  [trait ...]     step 3d: MR for every fiber cluster x trait
#   bash submit_mr.sh agg                  step 4: merge results (run when everything is finished)
#
# Options (before the traits):
#   --max N      total MR jobs running at once (default 150)
#   --array SPEC fiber clusters to run, e.g. --array 7,19,200-210 (default 1-N_FIBERS)
#   --after ID   start only after job ID succeeded
# With no trait named, all traits in mr_config.sh are used.
# =====================================================================
set -euo pipefail
export MR_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${MR_DIR}/mr_config.sh"

[[ $# -ge 1 ]] || { sed -n '3,16p' "$0"; exit 1; }
MODE=$1; shift
MAX=150; ARRAY="1-${N_FIBERS}"; AFTER=""; SEL=()
while [[ $# -gt 0 ]]; do
  case "$1" in
    --max) MAX=$2; shift 2 ;;
    --array) ARRAY=$2; shift 2 ;;
    --after) AFTER=$2; shift 2 ;;
    -*) echo "Unknown option $1" >&2; exit 1 ;;
    *) SEL+=("$1"); shift ;;
  esac
done
[[ ${#SEL[@]} -gt 0 ]] || SEL=("${TRAITS[@]}")

# position (1-based) of each selected trait in TRAITS
idx_of() { local i; for i in "${!TRAITS[@]}"; do [[ "${TRAITS[$i]}" == "$1" ]] && { echo $((i + 1)); return; }; done
           echo "Unknown trait '$1' (known: ${TRAITS[*]})" >&2; exit 1; }
IDX=(); for t in "${SEL[@]}"; do IDX+=("$(idx_of "$t")"); done
IDX_CSV=$(IFS=,; echo "${IDX[*]}")

mkdir -p "${LOG_DIR}/prep" "${HARM_DIR}" "${OUT_ROOT}"
DEP=(); [[ -n "$AFTER" ]] && DEP=(--dependency=afterok:"$AFTER" --kill-on-invalid-dep=yes)
COMMON=(--parsable --export=ALL,MR_DIR="${MR_DIR}")

case "$MODE" in
  prep)
    J1=$(sbatch "${COMMON[@]}" ${DEP[@]+"${DEP[@]}"} -o "${LOG_DIR}/prep/instruments_%j.out" -e "${LOG_DIR}/prep/instruments_%j.err" "${MR_DIR}/step3b_instruments.slurm")
    J2=$(sbatch "${COMMON[@]}" --dependency=afterok:"$J1" --kill-on-invalid-dep=yes --array="${IDX_CSV}" \
           -o "${LOG_DIR}/prep/harmonize_%a.out" -e "${LOG_DIR}/prep/harmonize_%a.err" "${MR_DIR}/step3c_harmonize.slurm")
    echo "instruments: job $J1 ; harmonize (${SEL[*]}): job $J2 ; logs: ${LOG_DIR}/prep" >&2
    echo "$J2" ;;
  harmonize)
    J2=$(sbatch "${COMMON[@]}" ${DEP[@]+"${DEP[@]}"} --array="${IDX_CSV}" \
           -o "${LOG_DIR}/prep/harmonize_%a.out" -e "${LOG_DIR}/prep/harmonize_%a.err" "${MR_DIR}/step3c_harmonize.slurm")
    echo "harmonize (${SEL[*]}): job $J2 ; logs: ${LOG_DIR}/prep" >&2
    echo "$J2" ;;
  run)
    PER=$(( MAX / ${#SEL[@]} )); [[ $PER -ge 1 ]] || PER=1
    for t in "${SEL[@]}"; do
      [[ -n "$AFTER" || -s "${HARM_DIR}/${t}.tsv" ]] || { echo "[skip] ${t}: ${HARM_DIR}/${t}.tsv missing (harmonize it first)" >&2; continue; }
      mkdir -p "${LOG_DIR}/${t}" "${OUT_ROOT}/${t}"
      J=$(sbatch --parsable --export=ALL,MR_DIR="${MR_DIR}",TRAIT="${t}" ${DEP[@]+"${DEP[@]}"} --job-name="mr_${t}" --array="${ARRAY}%${PER}" \
            -o "${LOG_DIR}/${t}/mr_${t}_%a.out" -e "${LOG_DIR}/${t}/mr_${t}_%a.err" "${MR_DIR}/step3d_mr.slurm")
      echo "${t}: job ${J} (array ${ARRAY}, ${PER} at a time)" >&2
    done
    echo "At most $(( PER * ${#SEL[@]} )) MR jobs run at once. Results: ${OUT_ROOT}/<trait>/ ; logs: ${LOG_DIR}/<trait>/" >&2 ;;
  agg)
    load_r; Rscript "${MR_DIR}/step4_agg_res.R" "${SEL[@]}" ;;
  *) echo "Unknown mode '$MODE' (prep|harmonize|run|agg)" >&2; exit 1 ;;
esac
