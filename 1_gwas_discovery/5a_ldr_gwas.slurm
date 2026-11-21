#!/bin/bash
# =============================================================================
# 5a_ldr_gwas.slurm
#
# Run LDR-level genome-wide association study (GWAS) using HEIG.
# Outputs per-cluster GWAS summary statistics (.parquet) to be preprocessed
# in step 5b before voxel-level analysis.
#
# Uses LOCO predictions from step 4 to correct for relatedness/stratification.
# The full imputed genotype matrix is used here (vs. the smaller SNP set for
# relatedness).
#
# Array index = fiber cluster ID (1-indexed).
# =============================================================================
#SBATCH --job-name=ldr_gwas
#SBATCH --output=logs/gwas_%a.out
#SBATCH --error=logs/gwas_%a.err
#SBATCH --mem=45G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --time=15:00:00
#SBATCH --array=1-430%20

# ── User-configurable paths ────────────────────────────────────────────────
MODALITY="dti_FA"
HEIG_DIR="/path/to/heig"
LDR_DIR="/path/to/output/ldrs"
REL_DIR="/path/to/output/rel"
GWAS_OUT="/path/to/output/gwas"
COVAR_FILE="/path/to/pheno_covar/covariates.txt"
CAT_COVAR="sex,assess_center"
GENO_MT="/path/to/genotype/ukb_imputed_allchr_qc.mt"   # Full imputed genotype matrix
SPARK_CONF="/path/to/spark_config.json"
TMP_ROOT="/path/to/tmp"
NLDR_FILE="/path/to/output/selLDRs.txt"
# ──────────────────────────────────────────────────────────────────────────

module purge
module load anaconda
conda activate heig

export JAVA_HOME="${CONDA_PREFIX}"
export PATH="${CONDA_PREFIX}/bin:${PATH}"

mkdir -p "${GWAS_OUT}"

export TMPDIR="${TMP_ROOT}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
export SPARK_LOCAL_DIRS="${TMPDIR}/spark"
mkdir -p "${SPARK_LOCAL_DIRS}"

IDX=${SLURM_ARRAY_TASK_ID}
N_LDRS=$(awk -F',' -v i=${IDX} 'NR == i+1 {print $2}' "${NLDR_FILE}")

python "${HEIG_DIR}/heig.py" \
    --gwas \
    --out          "${GWAS_OUT}/UKB_${MODALITY}_phase16_${IDX}" \
    --n-ldrs       "${N_LDRS}" \
    --covar        "${COVAR_FILE}" \
    --cat-covar-list "${CAT_COVAR}" \
    --ldrs         "${LDR_DIR}/UKB_${MODALITY}_phase16_${IDX}_ldr_top"* \
    --geno-mt      "${GENO_MT}" \
    --grch37 \
    --spark-conf   "${SPARK_CONF}" \
    --loco-preds   "${REL_DIR}/UKB_${MODALITY}_phase16_${IDX}_ldr_loco_preds.h5" \
    --threads      4 \
    --not-save-genotype-data

rm -rf "${GWAS_OUT}/UKB_${MODALITY}_phase16_${IDX}_tmp/"
rm -rf "${GWAS_OUT}/UKB_${MODALITY}_phase16_${IDX}_spark/"
rm -rf "${TMPDIR}"
