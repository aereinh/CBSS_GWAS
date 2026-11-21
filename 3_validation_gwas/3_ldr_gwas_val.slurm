#!/bin/bash
# =============================================================================
# 3_ldr_gwas_val.slurm
#
# Run LDR-level GWAS on the validation (Phase 7) sample.
# Mirrors 5a_ldr_gwas.slurm but uses Phase 7 LDRs and LOCO predictions.
#
# Array index = fiber cluster ID (1-indexed).
# =============================================================================
#SBATCH --job-name=gwas_val
#SBATCH --output=logs/gwas_val_%a.out
#SBATCH --error=logs/gwas_val_%a.err
#SBATCH --mem=45G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --time=15:00:00
#SBATCH --array=1-430%20

# ── User-configurable paths ────────────────────────────────────────────────
MODALITY="dti_FA"
HEIG_DIR="/path/to/heig"
LDR_DIR_VAL="/path/to/output/ldrs_val"
REL_DIR_VAL="/path/to/output/rel_val"
GWAS_OUT_VAL="/path/to/output/gwas_val"
COVAR_FILE="/path/to/pheno_covar/covariates.txt"
CAT_COVAR="sex,assess_center"
GENO_MT="/path/to/genotype/ukb_imputed_allchr_qc.mt"
SPARK_CONF="/path/to/spark_config.json"
TMP_ROOT="/path/to/tmp"
NLDR_FILE_VAL="/path/to/output/selLDRs_val.txt"
# ──────────────────────────────────────────────────────────────────────────

module purge
module load anaconda
conda activate heig

export JAVA_HOME="${CONDA_PREFIX}"
export PATH="${CONDA_PREFIX}/bin:${PATH}"

mkdir -p "${GWAS_OUT_VAL}"

export TMPDIR="${TMP_ROOT}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
export SPARK_LOCAL_DIRS="${TMPDIR}/spark"
mkdir -p "${SPARK_LOCAL_DIRS}"

IDX=${SLURM_ARRAY_TASK_ID}
N_LDRS=$(awk -F',' -v i=${IDX} 'NR == i+1 {print $2}' "${NLDR_FILE_VAL}")

python "${HEIG_DIR}/heig.py" \
    --gwas \
    --out          "${GWAS_OUT_VAL}/UKB_${MODALITY}_phase7_${IDX}" \
    --n-ldrs       "${N_LDRS}" \
    --covar        "${COVAR_FILE}" \
    --cat-covar-list "${CAT_COVAR}" \
    --ldrs         "${LDR_DIR_VAL}/UKB_${MODALITY}_phase7_${IDX}_ldr_top"* \
    --geno-mt      "${GENO_MT}" \
    --grch37 \
    --spark-conf   "${SPARK_CONF}" \
    --loco-preds   "${REL_DIR_VAL}/UKB_${MODALITY}_phase7_${IDX}_ldr_loco_preds.h5" \
    --threads      4 \
    --not-save-genotype-data

rm -rf "${GWAS_OUT_VAL}/UKB_${MODALITY}_phase7_${IDX}_tmp/"
rm -rf "${GWAS_OUT_VAL}/UKB_${MODALITY}_phase7_${IDX}_spark/"
rm -rf "${TMPDIR}"
