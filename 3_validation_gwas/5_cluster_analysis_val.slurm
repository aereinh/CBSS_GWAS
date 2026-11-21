#!/bin/bash
# =============================================================================
# 5_cluster_analysis_val.slurm
#
# Run cluster analysis on validation (Phase 7) GWAS results to generate a
# null cluster-size distribution for replication screening.
#
# Runs N_BATCHES independent bootstrap batches per cluster to build a robust
# null.  Each batch is run sequentially within the job to avoid excessive
# job submissions; results are written to the final output directory only
# after a successful run (enabling safe reruns on failure).
#
# Array index = fiber cluster ID (1-indexed).
# =============================================================================
#SBATCH --job-name=cluster_val
#SBATCH --output=logs/cluster_val_%a.out
#SBATCH --error=logs/cluster_val_%a.err
#SBATCH --mem=60G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --time=20:00:00
#SBATCH --array=1-430%30

# ── User-configurable paths ────────────────────────────────────────────────
MODALITY="dti_FA"
HEIG_DIR="/path/to/heig"
FPCA_DIR="/path/to/output/fpca"           # Discovery fPCA bases
LDR_DIR_VAL="/path/to/output/ldrs_val"
REL_DIR_VAL="/path/to/output/rel_val"
CLUSTER_OUT_VAL="/path/to/output/cluster_val"
COVAR_FILE="/path/to/pheno_covar/covariates.txt"
CAT_COVAR="sex,assess_center"
GENO_MT="/path/to/genotype/ukb_cal_allchr_cleaned.mt"
SPARK_CONF="/path/to/spark_config.json"
TMP_ROOT="/path/to/tmp"
NLDR_FILE_VAL="/path/to/output/selLDRs_val.txt"

SIG_THRESH=1e-3        # Screening threshold for validation
N_BATCHES=5            # Number of bootstrap batches to run
N_BOOTSTRAP=10         # Bootstraps per batch
# ──────────────────────────────────────────────────────────────────────────

module purge
module load anaconda
conda activate heig

mkdir -p "${CLUSTER_OUT_VAL}"

IDX=${SLURM_ARRAY_TASK_ID}
N_LDRS=$(awk -F',' -v i=${IDX} 'NR == i+1 {print $2}' "${NLDR_FILE_VAL}")

# Per-job scratch space
JOB_TMP="${TMP_ROOT}/cluster_val_${SLURM_JOB_ID}_${IDX}"
mkdir -p "${JOB_TMP}"
trap "rm -rf '${JOB_TMP}'" EXIT

for BATCH in $(seq 1 ${N_BATCHES}); do
    FINAL_TXT="${CLUSTER_OUT_VAL}/UKB_${MODALITY}_phase7_${IDX}_${BATCH}.txt"

    # Skip if already complete
    if [[ -s "${FINAL_TXT}" ]]; then
        echo "[SKIP] Batch ${BATCH} already exists: ${FINAL_TXT}"
        continue
    fi

    BATCH_TMP="${JOB_TMP}/batch_${BATCH}"
    mkdir -p "${BATCH_TMP}/spark"

    export TMPDIR="${BATCH_TMP}"
    export SPARK_LOCAL_DIRS="${BATCH_TMP}/spark"
    export _JAVA_OPTIONS="-Djava.io.tmpdir=${BATCH_TMP}"

    RUN_PREFIX="${BATCH_TMP}/UKB_${MODALITY}_phase7_${IDX}_${BATCH}"

    echo "Running batch ${BATCH} for cluster ${IDX} (N_LDRS=${N_LDRS})"

    python "${HEIG_DIR}/heig.py" \
        --cluster \
        --out          "${RUN_PREFIX}" \
        --n-ldrs       "${N_LDRS}" \
        --bases        "${FPCA_DIR}/UKB_${MODALITY}_phase16_${IDX}_bases_top"* \
        --ldr-cov      "${LDR_DIR_VAL}/UKB_${MODALITY}_phase7_${IDX}_ldr_cov_top"* \
        --covar        "${COVAR_FILE}" \
        --cat-covar-list "${CAT_COVAR}" \
        --ldrs         "${LDR_DIR_VAL}/UKB_${MODALITY}_phase7_${IDX}_ldr_top"* \
        --geno-mt      "${GENO_MT}" \
        --grch37 \
        --spark-conf   "${SPARK_CONF}" \
        --loco-preds   "${REL_DIR_VAL}/UKB_${MODALITY}_phase7_${IDX}_ldr_loco_preds.h5" \
        --sig-thresh   "${SIG_THRESH}" \
        --n-bootstrap  "${N_BOOTSTRAP}"

    if [[ -s "${RUN_PREFIX}.txt" ]]; then
        cp -f "${RUN_PREFIX}.txt" "${FINAL_TXT}"
        echo "[OK] Wrote: ${FINAL_TXT}"
    else
        echo "[ERROR] Missing output for cluster ${IDX} batch ${BATCH}" >&2
        exit 1
    fi

    rm -rf "${BATCH_TMP}"
done
