#!/bin/bash
# =====================================================================
# Submit the whole ABCD pipeline for one modality as a chain of Slurm jobs:
#   format -> read -> ldr -> select -> gwas
# Each step starts only after every task of the previous one succeeded.
#
#   bash run_all.sh <MD|AD|FA> [first_step]
#
#   first_step  start the chain here (format|read|ldr|select|gwas); default format.
#               Steps already done are skipped per cluster anyway, so starting
#               earlier than needed only costs queue time.
# =====================================================================
set -euo pipefail
cd "$(dirname "${BASH_SOURCE[0]}")"
MOD=${1:?usage: bash run_all.sh <MD|AD|FA> [first_step]}
FIRST=${2:-format}
STEPS=(format read ldr select gwas)

start=-1
for i in "${!STEPS[@]}"; do [[ "${STEPS[$i]}" == "$FIRST" ]] && start=$i; done
[[ $start -ge 0 ]] || { echo "Unknown step '$FIRST'" >&2; exit 1; }

PREV=""
for ((i=start; i<${#STEPS[@]}; i++)); do
  S=${STEPS[$i]}
  if [[ -z "$PREV" ]]; then J=$(bash submit.sh "$MOD" "$S")
  else                      J=$(bash submit.sh "$MOD" "$S" --after "$PREV"); fi
  echo "$S: $J"
  PREV=$J
done
echo "Chain submitted for $MOD. Monitor with: squeue -u \$USER ; bash check_inputs.sh $MOD"
