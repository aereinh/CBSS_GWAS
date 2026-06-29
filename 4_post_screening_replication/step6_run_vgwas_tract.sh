#!/bin/bash
# =============================================================================
# step6_run_vgwas_tract.sh
#
# Run targeted voxel GWAS in the validation cohort for a single fiber cluster,
# restricting to candidate SNPs near discovery lead loci for the tract that
# cluster belongs to (from step 5).
#
# Called by step6_run.slurm with --tract, --cluster, --n-ldrs, --modality.
# =============================================================================

# Parse named arguments
while [[ $# -gt 0 ]]; do
    case "$1" in
        --tract)    TRACT="$2";    shift 2 ;;
        --cluster)  CLUSTER="$2";  shift 2 ;;
        --n-ldrs)   N_LDRS="$2";   shift 2 ;;
        --modality) MODALITY="$2"; shift 2 ;;
        *)          echo "Unknown argument: $1"; exit 1 ;;
    esac
done

# ── User-configurable paths ────────────────────────────────────────────────
HEIG_DIR="/path/to/heig"
MAIN_DIR="/path/to/output/output_${MODALITY}"
FPCA_DIR="${MAIN_DIR}/fpca"
LDR_DIR_VAL="${MAIN_DIR}/ldrs_val"
SUMSTATS_DIR_VAL="${MAIN_DIR}/sumstats_val"
REP_OUT="${MAIN_DIR}/replication/tract"
SIG_THRESH=1e-3    # Relaxed threshold for validation (adjust as needed)
# ──────────────────────────────────────────────────────────────────────────

mkdir -p "${REP_OUT}"

CANDIDATE_SNPS="${MAIN_DIR}/post_screen/disc/tract/UKB_${MODALITY}_allphases_Tract_${TRACT}_candidate_sig_snps.txt"

if [ ! -f "${CANDIDATE_SNPS}" ]; then
    echo "No candidate SNP file for tract ${TRACT}. Skipping cluster ${CLUSTER}."
    exit 0
fi

python "${HEIG_DIR}/heig.py" \
    --voxel-gwas \
    --out          "${REP_OUT}/UKB_${MODALITY}_phase6_${CLUSTER}_vgwas" \
    --n-ldrs       "${N_LDRS}" \
    --ldr-sumstats "${SUMSTATS_DIR_VAL}/UKB_${MODALITY}_phase6_${CLUSTER}" \
    --bases        "${FPCA_DIR}/UKB_${MODALITY}_allphases_${CLUSTER}_bases_top"* \
    --ldr-cov      "${LDR_DIR_VAL}/UKB_${MODALITY}_phase6_${CLUSTER}_ldr_cov_top"* \
    --extract      "${CANDIDATE_SNPS}" \
    --threads      4 \
    --sig-thresh   "${SIG_THRESH}"
