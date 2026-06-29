#!/bin/bash
# =============================================================================
# 2_relatedness_val.slurm
#
# Compute LOCO predictions for the validation (Phase 6) sample.
# Mirrors 4_relatedness.slurm but uses validation LDRs and the smaller
# genotype reference panel.
#
# Array index = fiber cluster ID (1-indexed).
# =============================================================================
#SBATCH --job-name=rel_val
#SBATCH --output=logs/rel_val_%a.out
#SBATCH --error=logs/rel_val_%a.err
#SBATCH --mem=45G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --time=10:00:00
#SBATCH --array=1-430%20

# ── User-configurable paths ────────────────────────────────────────────────
MODALITY="dti_FA"
HEIG_DIR="/path/to/heig"
LDR_DIR_VAL="/path/to/output/ldrs_val"
REL_OUT_VAL="/path/to/output/rel_val"
COVAR_FILE="/path/to/pheno_covar/covariates.txt"
CAT_COVAR="sex,assess_center"
GENO_MT="/path/to/genotype/geno_1ksnps_qc"
SPARK_CONF="/path/to/spark_config.json"
TMP_ROOT="/path/to/tmp"
NLDR_FILE_VAL="/path/to/output/selLDRs_val.txt"
# ──────────────────────────────────────────────────────────────────────────

module purge
module load anaconda
conda activate heig

export JAVA_HOME="${CONDA_PREFIX}"
export PATH="${CONDA_PREFIX}/bin:${PATH}"

mkdir -p "${REL_OUT_VAL}"

export TMPDIR="${TMP_ROOT}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
export SPARK_LOCAL_DIRS="${TMPDIR}/spark"
mkdir -p "${SPARK_LOCAL_DIRS}"

IDX=${SLURM_ARRAY_TASK_ID}
N_LDRS=$(awk -F',' -v i=${IDX} 'NR == i+1 {print $2}' "${NLDR_FILE_VAL}")

python "${HEIG_DIR}/heig.py" \
    --relatedness \
    --out          "${REL_OUT_VAL}/UKB_${MODALITY}_phase6_${IDX}" \
    --n-ldrs       "${N_LDRS}" \
    --covar        "${COVAR_FILE}" \
    --cat-covar-list "${CAT_COVAR}" \
    --ldrs         "${LDR_DIR_VAL}/UKB_${MODALITY}_phase6_${IDX}_ldr_top"* \
    --geno-mt      "${GENO_MT}" \
    --grch37 \
    --spark-conf   "${SPARK_CONF}" \
    --threads      4

rm -rf "${REL_OUT_VAL}/UKB_${MODALITY}_phase6_${IDX}_tmp/"
rm -rf "${REL_OUT_VAL}/UKB_${MODALITY}_phase6_${IDX}_spark/"
rm -rf "${TMPDIR}"
