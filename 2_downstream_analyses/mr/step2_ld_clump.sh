#!/bin/bash
# =============================================================================
# step2_ld_clump.sh
#
# LD-clump the screened LDR sumstats from step 1 using PLINK 1.9.
# Produces a .clumped file with independent loci (r²<0.01, 10 Mb window).
#
# Usage:
#   bash step2_ld_clump.sh --roi <cluster_id>
#
# Example:
#   bash step2_ld_clump.sh --roi 42
# =============================================================================

set -euo pipefail

# ── User-configurable ─────────────────────────────────────────────────────────
LD_REF_BFILE="/path/to/ld_reference/g1000_eur"   # PLINK bfile prefix for LD reference panel
SCREENED_DIR="/path/to/output/mr/screened"        # Directory of step 1 .txt outputs
MODALITY="dti_FA"
# ─────────────────────────────────────────────────────────────────────────────

while [[ $# -gt 0 ]]; do
    case "$1" in
        --roi)  roi="$2"; shift ;;
        *)      echo "Unknown argument: $1"; echo "Usage: bash step2_ld_clump.sh --roi <cluster_id>"; exit 1 ;;
    esac
    shift
done

if [[ -z "${roi:-}" ]]; then
    echo "Error: --roi is required"
    exit 1
fi

module load plink/1.90b6

SUMSTATS="${SCREENED_DIR}/UKB_${MODALITY}_allphases_${roi}.txt"
OUT_PREFIX="${SCREENED_DIR}/UKB_${MODALITY}_allphases_${roi}_clump_r2_0.01"

if [[ ! -f "${SUMSTATS}" ]]; then
    echo "No screened sumstats found: ${SUMSTATS}. Skipping."
    exit 0
fi

echo "LD clumping cluster ${roi} ..."

plink \
    --bfile "${LD_REF_BFILE}" \
    --clump "${SUMSTATS}" \
    --clump-kb 10000 \
    --clump-r2 0.01 \
    --out "${OUT_PREFIX}"

echo "Done: ${OUT_PREFIX}.clumped"

# ── Get allele frequencies for the clumped SNPs ───────────────────────────────
SNPLIST="${SCREENED_DIR}/UKB_${MODALITY}_allphases_${roi}_snplist.txt"
FREQ_PREFIX="${SCREENED_DIR}/UKB_${MODALITY}_allphases_${roi}_freq"

if [[ -f "${SNPLIST}" ]]; then
    plink \
        --bfile "${LD_REF_BFILE}" \
        --extract "${SNPLIST}" \
        --freq \
        --out "${FREQ_PREFIX}"
    echo "Frequencies written to: ${FREQ_PREFIX}.frq"
fi
