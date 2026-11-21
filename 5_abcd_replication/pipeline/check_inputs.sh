#!/bin/bash
# =====================================================================
# Check that every input and output a modality needs is where config.sh says.
# Run on the login node before submitting (takes a few seconds):
#   bash check_inputs.sh FA
# =====================================================================
PIPE_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
source "${PIPE_DIR}/config.sh"
[[ $# -eq 1 ]] || { echo "Usage: bash check_inputs.sh <FA|MD|AD>"; exit 1; }
set_modality "$1" || exit 1

ok()   { printf "  [ OK ] %s\n" "$*"; }
bad()  { printf "  [MISS] %s\n" "$*"; }
chk()  { if [[ -e "$2" ]]; then ok "$1: $2"; else bad "$1: $2"; fi; }
count() {  # count <label> <pattern with IDX>
  local label=$1 pat=$2 n=0 missing=()
  for IDX in $(seq 1 "${N_CLUSTERS}"); do
    if compgen -G "${pat//\{IDX\}/$IDX}" > /dev/null; then n=$((n+1)); else missing+=("$IDX"); fi
  done
  local first="${missing[*]:0:15}"
  if [[ $n -eq ${N_CLUSTERS} ]]; then ok "$label: ${n}/${N_CLUSTERS}"
  else printf "  [%4d] %s: %d/%d  (missing e.g. %s)\n" "$n" "$label" "$n" "${N_CLUSTERS}" "${first}"; fi
}

echo "== ${MOD}  (UKB name ${UKB_MOD}; results -> ${RUN_NAME})"
echo "-- shared"
chk "heig.py"     "${HEIG_DIR}/heig.py"
chk "spark conf"  "${SPARK_CONF}"
chk "genotype mt" "${GENO_MT}"
chk "covariates"  "${COVAR}"
if [[ -e "${COVAR}" ]]; then
  HDR=$(head -n 1 "${COVAR}")
  for c in FID IID ${CAT_COVARS//,/ }; do
    grep -qw "$c" <<< "$HDR" && ok "covariate column '$c' present" || bad "covariate column '$c' missing in header"
  done
  echo "        $(($(wc -l < "${COVAR}") - 1)) subjects, $(wc -w <<< "$HDR") columns"
fi

echo "-- step 0 (format; skip if all txt files exist)"
chk "raw data dir"     "${RAW_DIR}"
chk "subject IDs"      "${RAW_DIR}/ABCD.txt"
chk "first .mat file"  "${RAW_DIR}/1_Aligned.mat"
chk "fiber atlas"      "${CLUST_CSV}"
count "image txt   ${IMG_TXT_DIR}/fibersClust_*.txt" "${IMG_TXT_DIR}/fibersClust_{IDX}.txt"

echo "-- step 1 (read)"
count "coords      ${COORD_DIR}/coords_fibersClust_*.txt" "${COORD_DIR}/coords_fibersClust_{IDX}.txt"
count "images h5   ${IMG_DIR}/ABCD_${MOD}_*_images.h5" "${IMG_DIR}/ABCD_${MOD}_{IDX}_images.h5"

echo "-- step 3 (ldr / select)"
count "UKB bases   ${UKB_BASES_PREFIX}_*_bases_*" "${UKB_BASES_PREFIX}_{IDX}_bases_*"
count "UKB fPCA log ${UKB_FPCA_LOG_PREFIX}_*.log" "${UKB_FPCA_LOG_PREFIX}_{IDX}.log"
count "ABCD LDRs   ${LDR_DIR}/ABCD_${MOD}_*_ldr_top*.txt" "${LDR_DIR}/ABCD_${MOD}_{IDX}_ldr_top*.txt"
[[ -n "${REUSE_SEL_LDRS}" ]] && chk "selLDRs to reuse" "${REUSE_SEL_LDRS}"
chk "selLDRs_val.txt" "${SEL_LDRS}"

echo "-- steps 4-5 (gwas)"
count "LOCO preds  ${REL_DIR}/ABCD_${MOD}_*_ldr_loco_preds.h5" "${REL_DIR}/ABCD_${MOD}_{IDX}_ldr_loco_preds.h5"
count "sumstats    ${SUMSTATS_DIR}/ABCD_${MOD}_*.snpinfo" "${SUMSTATS_DIR}/ABCD_${MOD}_{IDX}.snpinfo*"
