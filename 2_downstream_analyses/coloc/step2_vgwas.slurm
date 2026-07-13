#!/bin/bash
# =============================================================================
# step2_vgwas.slurm
#
# Array job: run targeted single-voxel GWAS for each locus using HEIG.
# Reads the per-locus best-voxel CSV produced by step 1 and runs HEIG
# --voxel-gwas restricted to a ±250 kb window around the locus centre.
#
# Requires:
#   - HEIG installation with heig.py
#   - selLDRs.txt: two-column CSV (cluster_id, n_ldrs) for LDR counts
# =============================================================================
#SBATCH --job-name=coloc_vgwas
#SBATCH --output=logs/coloc_vgwas_%a.out
#SBATCH --error=logs/coloc_vgwas_%a.err
#SBATCH --partition=general
#SBATCH --time=00:30:00
#SBATCH --mem=12G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-315%50      # adjust to number of loci

# ── User-configurable paths ───────────────────────────────────────────────────
MODALITY="dti_FA"
HEIG_DIR="/path/to/heig"                           # HEIG installation directory
SEL_DIR="/path/to/coloc/sel_voxels"               # Step 1 output directory
LDR_SUMSTATS_DIR="/path/to/sumstats"              # LDR sumstats (prefix per cluster)
BASES_DIR="/path/to/fpca"                         # fPCA bases directory
LDR_COV_DIR="/path/to/ldr"                        # LDR covariance directory
SEL_LDRS_FILE="/path/to/selLDRs.txt"             # Two-column CSV: cluster_id,n_ldrs
OUT_DIR="/path/to/coloc/voxelgwas"
WINDOW=250000
# ─────────────────────────────────────────────────────────────────────────────

module purge
module load anaconda
conda activate heig

LOCUS_ID=${SLURM_ARRAY_TASK_ID}
mkdir -p "${OUT_DIR}" logs

# Find per-locus selected-voxel CSV from step 1
CSV_GLOB="${SEL_DIR}/locus_${LOCUS_ID}_*.csv"
CSV_FP="$(ls -1 ${CSV_GLOB} 2>/dev/null | head -n 1 || true)"
if [[ -z "${CSV_FP}" ]]; then
    echo "No selected-voxel file found for locus ${LOCUS_ID}. Skipping."
    exit 0
fi

# Parse fields from CSV
read_csv_field() {
    local fp="$1" col="$2"
    awk -F',' -v col="$col" '
        NR==1 { for(i=1;i<=NF;i++) { gsub(/^[ \t"]+|[ \t"]+$/, "", $i); if($i==col) c=i }; next }
        NR==2 { if(!c) exit 2; v=$c; gsub(/^[ \t"]+|[ \t"]+$/, "", v); print v; exit }
    ' "$fp"
}

CLUSTER="$(read_csv_field "${CSV_FP}" "cluster")"
VOX="$(read_csv_field    "${CSV_FP}" "INDEX")"
CHR="$(read_csv_field    "${CSV_FP}" "locus_chr")"
POS="$(read_csv_field    "${CSV_FP}" "locus_bp")"
STATUS="$(read_csv_field "${CSV_FP}" "status")"

if [[ "${STATUS}" != "OK" ]]; then
    echo "Locus ${LOCUS_ID} status=${STATUS}. Skipping."
    exit 0
fi

if [[ -z "${CLUSTER}" || -z "${VOX}" || -z "${CHR}" || -z "${POS}" ]]; then
    echo "Failed to parse CSV fields for locus ${LOCUS_ID}. Skipping."
    exit 1
fi

CHR="${CHR#chr}"
START=$(( POS - WINDOW )); (( START < 1 )) && START=1
END=$(( POS + WINDOW ))
CHR_INTERVAL="${CHR}:${START},${CHR}:${END}"

# Look up number of LDRs for this cluster
N_LDRS=$(awk -F',' -v c="${CLUSTER}" 'NR>1 && $1==c { print $2; exit }' "${SEL_LDRS_FILE}")
if [[ -z "${N_LDRS}" ]]; then
    echo "Could not find N_LDRS for cluster ${CLUSTER}. Skipping."
    exit 1
fi

OUT_PREFIX="${OUT_DIR}/UKB_${MODALITY}_LOC${LOCUS_ID}_CLST${CLUSTER}_VOX${VOX}"

echo "Locus ${LOCUS_ID} | cluster=${CLUSTER} | vox=${VOX} | chr${CHR}:${START}-${END} | N_LDRS=${N_LDRS}"

python "${HEIG_DIR}/heig.py" \
    --voxel-gwas \
    --out          "${OUT_PREFIX}" \
    --n-ldrs       "${N_LDRS}" \
    --ldr-sumstats "${LDR_SUMSTATS_DIR}/UKB_${MODALITY}_allphases_${CLUSTER}" \
    --bases        "${BASES_DIR}/UKB_${MODALITY}_allphases_${CLUSTER}_bases_top"* \
    --ldr-cov      "${LDR_COV_DIR}/UKB_${MODALITY}_allphases_${CLUSTER}_ldr_cov_top"* \
    --chr-interval "${CHR_INTERVAL}" \
    --voxels       "${VOX}"

echo "Done: locus ${LOCUS_ID}"
