#!/bin/bash
# Progress of the MR analysis: what is in place, and finished fiber clusters per trait.
#   bash check_mr.sh
MR_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "${MR_DIR}/mr_config.sh"
shopt -s nullglob
nf() { local f=( $1 ); echo ${#f[@]}; }
n_exp=$(nf "${SCREEN_DIR}/${EXP_PREFIX}*_clump_r2_0.01.clumped")
echo "fiber clusters with clumped instruments: ${n_exp}/${N_FIBERS}"
[[ -s "${HARM_DIR}/instruments.tsv" ]] && echo "instruments: $(($(wc -l < "${HARM_DIR}/instruments.tsv") - 1))" || echo "instruments: [MISS] run 'bash submit_mr.sh prep'"
printf "%-18s %-6s %-12s %-10s %-10s %s\n" trait raw harmonized results no_result failed_logs
for t in "${TRAITS[@]}"; do
  pat=$(grep -A2 "^  ${t} = list" "${MR_DIR}/trait_config.R" | grep -o 'file = "[^"]*"' | head -1 | cut -d'"' -f2)
  raw=0; for f in "${RAW_DIR}/${pat}"*; do [[ "$f" =~ \.(zip|tar|tgz)$ ]] || raw=$((raw+1)); done
  harm="-"; [[ -s "${HARM_DIR}/${t}.tsv" ]] && harm=$(($(wc -l < "${HARM_DIR}/${t}.tsv") - 1))
  res=( "${OUT_ROOT}/${t}"/*_MR_results.tsv ); don=( "${OUT_ROOT}/${t}"/*_MR.done )
  fail=$(grep -l -i -E "^Error|Execution halted|CANCELLED|oom" "${LOG_DIR}/${t}"/*.err /dev/null 2>/dev/null | wc -l)
  printf "%-18s %-6s %-12s %-10s %-10s %s\n" "$t" "$([[ $raw -ge 1 ]] && echo OK || echo MISS)" "$harm" "${#res[@]}" "${#don[@]}" "$fail"
done
echo "results + no_result should add up to ${N_FIBERS} per trait when a trait is finished."
