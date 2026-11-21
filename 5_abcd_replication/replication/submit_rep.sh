#!/bin/bash
# =====================================================================
# Submit the ABCD replication analysis (FA, site-adjusted).
#   bash submit_rep.sh all            1 voxel GWAS (430 tasks) -> 2-3 per tract (17 tasks) -> 4 aggregate
#   bash submit_rep.sh vgwas [--array 7,19]
#   bash submit_rep.sh tract [--array 0-16]
#   bash submit_rep.sh agg
# Each later step can be given --after JOBID.
# =====================================================================
set -euo pipefail
export REP_SCRIPTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${REP_SCRIPTS}/rep_config.sh"
[[ $# -ge 1 ]] || { sed -n '3,9p' "$0"; exit 1; }
MODE=$1; shift
ARRAY=""; AFTER=""
while [[ $# -gt 0 ]]; do case "$1" in --array) ARRAY=$2; shift 2;; --after) AFTER=$2; shift 2;; *) echo "Unknown option $1" >&2; exit 1;; esac; done
mkdir -p "${REP_LOG_DIR}" "${REP_DIR}" "${SUMMARY_DIR}"

sub() {  # sub <script> <logpattern> <after-jobid-or-empty> [array]
  local a=(--parsable --export=ALL,REP_SCRIPTS="${REP_SCRIPTS}",MOD="${MOD}" -o "${REP_LOG_DIR}/$2.out" -e "${REP_LOG_DIR}/$2.err")
  [[ -n "$3" ]] && a+=(--dependency=afterok:"$3" --kill-on-invalid-dep=yes)
  [[ -n "${4:-}" ]] && a+=(--array="$4")
  sbatch "${a[@]}" "${REP_SCRIPTS}/$1"
}
case "$MODE" in
  vgwas) J=$(sub 1_voxel_gwas.slurm "vgwas_%a" "$AFTER" "$ARRAY"); echo "vgwas: $J" ;;
  tract) J=$(sub 2-3_tract.slurm "tract_%a" "$AFTER" "$ARRAY");    echo "tract: $J" ;;
  agg)   J=$(sub 4_aggregate.slurm "aggregate_%j" "$AFTER");       echo "agg: $J" ;;
  all)
    J1=$(sub 1_voxel_gwas.slurm "vgwas_%a" "$AFTER" "$ARRAY"); echo "vgwas: $J1"
    J2=$(sub 2-3_tract.slurm "tract_%a" "$J1");                echo "tract: $J2"
    J3=$(sub 4_aggregate.slurm "aggregate_%j" "$J2");          echo "agg:   $J3" ;;
  *) echo "Unknown mode $MODE" >&2; exit 1 ;;
esac
echo "logs: ${REP_LOG_DIR} ; results: ${SUMMARY_DIR}"
