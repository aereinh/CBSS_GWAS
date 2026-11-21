#!/bin/bash
# =============================================================================
# 5d_cluster_analysis.slurm
#
# Run spatially-informed cluster analysis on voxel GWAS results using HEIG.
# Generates a null distribution of cluster sizes via bootstrapping, which is
# used in the post-screening step to determine genome-wide significant clusters.
#
# Note: this step is computationally intensive (Spark-backed).
# Set BATCH to run multiple independent bootstrap batches for robustness;
# results from all batches are combined in the post-screening step.
#
# Array index = fiber cluster ID (1-indexed).
# =============================================================================
#SBATCH --job-name=cluster
#SBATCH --output=logs/cluster_%a.out
#SBATCH --error=logs/cluster_%a.err
#SBATCH --mem=30G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --time=8:00:00
#SBATCH --array=1-430%30

# ── User-configurable paths ────────────────────────────────────────────────
MODALITY="dti_FA"
HEIG_DIR="/path/to/heig"
FPCA_DIR="/path/to/output/fpca"
LDR_DIR="/path/to/output/ldrs"
REL_DIR="/path/to/output/rel"
CLUSTER_OUT="/path/to/output/cluster"
COVAR_FILE="/path/to/pheno_covar/covariates.txt"
CAT_COVAR="sex,assess_center"
GENO_MT="/path/to/genotype/ukb_cal_allchr_cleaned.mt"
SPARK_CONF="/path/to/spark_config.json"
TMP_ROOT="/path/to/tmp"
NLDR_FILE="/path/to/output/selLDRs.txt"

# Bootstrap settings
N_BOOTSTRAP=10      # Bootstraps per batch
BATCH=1             # Batch number (run multiple batches for larger null)
SIG_THRESH=1e-5     # Cluster-level screening threshold
# ──────────────────────────────────────────────────────────────────────────

module purge
module load anaconda
conda activate heig

mkdir -p "${CLUSTER_OUT}"

export TMPDIR="${TMP_ROOT}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
export SPARK_LOCAL_DIRS="${TMPDIR}/spark"
mkdir -p "${SPARK_LOCAL_DIRS}"

IDX=${SLURM_ARRAY_TASK_ID}
N_LDRS=$(awk -F',' -v i=${IDX} 'NR == i+1 {print $2}' "${NLDR_FILE}")

python "${HEIG_DIR}/heig.py" \
    --cluster \
    --out          "${CLUSTER_OUT}/UKB_${MODALITY}_phase16_${IDX}_${BATCH}" \
    --n-ldrs       "${N_LDRS}" \
    --bases        "${FPCA_DIR}/UKB_${MODALITY}_phase16_${IDX}_bases_top"* \
    --ldr-cov      "${LDR_DIR}/UKB_${MODALITY}_phase16_${IDX}_ldr_cov_top"* \
    --covar        "${COVAR_FILE}" \
    --cat-covar-list "${CAT_COVAR}" \
    --ldrs         "${LDR_DIR}/UKB_${MODALITY}_phase16_${IDX}_ldr_top"* \
    --geno-mt      "${GENO_MT}" \
    --grch37 \
    --spark-conf   "${SPARK_CONF}" \
    --loco-preds   "${REL_DIR}/UKB_${MODALITY}_phase16_${IDX}_ldr_loco_preds.h5" \
    --sig-thresh   "${SIG_THRESH}" \
    --n-bootstrap  "${N_BOOTSTRAP}"

rm -rf "${CLUSTER_OUT}/UKB_${MODALITY}_phase16_${IDX}_${BATCH}_tmp/"
rm -rf "${TMPDIR}"
