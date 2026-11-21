#!/bin/bash
# =============================================================================
# 4_relatedness.slurm
#
# Compute leave-one-chromosome-out (LOCO) genetic predictions for each fiber
# cluster using HEIG.  These predictions account for population stratification
# and cryptic relatedness and are used as inputs to the GWAS step.
#
# Requires a genotype matrix table (.mt) in Hail/Spark format.
# See HEIG documentation for genotype file preparation.
#
# Array index = fiber cluster ID (1-indexed).
# =============================================================================
#SBATCH --job-name=relatedness
#SBATCH --output=logs/rel_%a.out
#SBATCH --error=logs/rel_%a.err
#SBATCH --mem=45G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --time=10:00:00
#SBATCH --array=1-430%20

# ── User-configurable paths ────────────────────────────────────────────────
MODALITY="dti_FA"
HEIG_DIR="/path/to/heig"
LDR_DIR="/path/to/output/ldrs"
REL_OUT="/path/to/output/rel"
COVAR_FILE="/path/to/pheno_covar/covariates.txt"
CAT_COVAR="sex,assess_center"
GENO_MT="/path/to/genotype/geno_1ksnps_qc"   # Hail matrix table for relatedness
SPARK_CONF="/path/to/spark_config.json"
TMP_ROOT="/path/to/tmp"
NLDR_FILE="/path/to/output/selLDRs.txt"      # Output of step 3b
# ──────────────────────────────────────────────────────────────────────────

module purge
module load anaconda
conda activate heig

export JAVA_HOME="${CONDA_PREFIX}"
export PATH="${CONDA_PREFIX}/bin:${PATH}"

mkdir -p "${REL_OUT}"

# Per-job temp dirs (Spark needs scratch space)
export TMPDIR="${TMP_ROOT}/${SLURM_JOB_ID}_${SLURM_ARRAY_TASK_ID}"
export SPARK_LOCAL_DIRS="${TMPDIR}/spark"
mkdir -p "${SPARK_LOCAL_DIRS}"

IDX=${SLURM_ARRAY_TASK_ID}
N_LDRS=$(awk -F',' -v i=${IDX} 'NR == i+1 {print $2}' "${NLDR_FILE}")

python "${HEIG_DIR}/heig.py" \
    --relatedness \
    --out          "${REL_OUT}/UKB_${MODALITY}_phase16_${IDX}" \
    --n-ldrs       "${N_LDRS}" \
    --covar        "${COVAR_FILE}" \
    --cat-covar-list "${CAT_COVAR}" \
    --ldrs         "${LDR_DIR}/UKB_${MODALITY}_phase16_${IDX}_ldr_top"* \
    --geno-mt      "${GENO_MT}" \
    --grch37 \
    --spark-conf   "${SPARK_CONF}" \
    --threads      4

# Clean up Spark and temp dirs
rm -rf "${REL_OUT}/UKB_${MODALITY}_phase16_${IDX}_tmp/"
rm -rf "${REL_OUT}/UKB_${MODALITY}_phase16_${IDX}_spark/"
rm -rf "${TMPDIR}"
