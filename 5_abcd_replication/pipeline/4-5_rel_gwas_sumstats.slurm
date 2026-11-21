#!/bin/bash
#SBATCH --job-name=abcd_gwas
#SBATCH --output=/path/to/UKB_PSC/logs/abcd/%x_%A_%a.out
#SBATCH --error=/path/to/UKB_PSC/logs/abcd/%x_%A_%a.err
#SBATCH --mem=35G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --time=12:00:00
#SBATCH --array=1-430%60

# Steps 4-5: relatedness (LOCO) -> LDR GWAS -> HEIG sumstats, all with COVAR (config.sh).
# Run with:  bash submit.sh <MOD> gwas
# Safe to resubmit: finished clusters are skipped, an existing LOCO file is reused,
# and an existing GWAS parquet is reused (only the sumstats step reruns).

: "${MOD:?Set MOD, e.g. bash submit.sh FA gwas}"
source "${PIPE_DIR:?submit with submit.sh}/config.sh"
set_modality "$MOD" || exit 1
load_heig_env
export JAVA_HOME="$CONDA_PREFIX"
export PATH="$CONDA_PREFIX/bin:$PATH"
set -eo pipefail

IDX=${SLURM_ARRAY_TASK_ID}
TAG=ABCD_${MOD}_${IDX}
LOCO=${REL_DIR}/${TAG}_ldr_loco_preds.h5
SUMSTATS_OUT=${SUMSTATS_DIR}/${TAG}

# sumstats outputs = any <TAG>.* file in SUMSTATS_DIR other than the log
sumstats_done() {
  local f
  for f in "${SUMSTATS_OUT}".*; do
    [[ -e "$f" && "$f" != *.log ]] && return 0
  done
  return 1
}

if sumstats_done && [[ -z "${OVERWRITE:-}" ]]; then
  echo "[SKIP] sumstats already exist for ${TAG} in ${SUMSTATS_DIR} (set OVERWRITE=1 to redo)"
  exit 0
fi

N_LDRS=$(n_ldrs_val "${IDX}")
LDR_FILES=( "${LDR_DIR}/${TAG}"_ldr_top*.txt )
LDRS=${LDR_FILES[0]}
[[ -e "${LDRS}" ]] || die "no LDR file ${LDR_DIR}/${TAG}_ldr_top*.txt (run the ldr step first)"

echo "MOD=${MOD} IDX=${IDX} N_LDRS=${N_LDRS}"
echo "LDRS=${LDRS}"
echo "COVAR=${COVAR}  CAT=${CAT_COVARS}"
echo "JAVA=$(which java)"; java -version 2>&1 | head -n 1

mkdir -p "${REL_DIR}" "${GWAS_DIR}" "${SUMSTATS_DIR}"
TMPBASE=${TMP_ROOT}/${SLURM_JOB_ID}_${IDX}
trap 'rm -rf "${TMPBASE}" 2>/dev/null || true' EXIT

# Remove HEIG/Spark tmp dirs without failing the job. On NFS (/work), files the
# exiting JVM still holds open become .nfsXXXX and can't be deleted yet
# ("Device or resource busy"); retry a few times, then warn and move on.
cleanup_tmp() {
  local d try
  for d in "$@"; do
    [[ -e "$d" ]] || continue
    for try in 1 2 3; do
      rm -rf "$d" 2>/dev/null && break
      sleep 20
    done
    if [[ -e "$d" ]]; then echo "[WARN] could not fully remove $d (NFS busy files); delete it later" >&2; fi
  done
  return 0
}

# ---- 4. Relatedness (skipped if the LOCO file already exists)
if [[ -s "${LOCO}" ]]; then
  echo "[SKIP] relatedness: ${LOCO} exists"
else
  export TMPDIR=${TMPBASE}/rel
  export SPARK_LOCAL_DIRS=${TMPDIR}/spark
  mkdir -p "${SPARK_LOCAL_DIRS}"

  python heig.py \
    --relatedness \
    --out "${REL_DIR}/${TAG}" \
    --n-ldrs "${N_LDRS}" \
    --covar "${COVAR}" \
    --cat-covar-list "${CAT_COVARS}" \
    --ldrs "${LDRS}" \
    --geno-mt "${GENO_MT}" \
    --grch37 \
    --spark-conf "${SPARK_CONF}" \
    --threads 4 \
    --maf-min 0.01 \
    --maf-max 0.5

  cleanup_tmp "${REL_DIR}/${TAG}_tmp" "${REL_DIR}/${TAG}_spark"
  [[ -s "${LOCO}" ]] || die "heig --relatedness did not write ${LOCO}"
fi

# HEIG truncates the long NDAR_INV IDs; restore them from the LDR file (idempotent)
python "${PIPE_DIR}/fix_h5_ids.py" --h5 "${LOCO}" --ids-txt "${LDRS}"

# ---- 5.1 LDR GWAS (skipped if the parquet already exists, e.g. after a failed tmp cleanup)
GWAS_PQ=${GWAS_DIR}/${TAG}.parquet
if [[ -d "${GWAS_PQ}" && -z "${OVERWRITE:-}" ]]; then
  echo "[SKIP] GWAS: ${GWAS_PQ} exists (set OVERWRITE=1 to redo)"
  cleanup_tmp "${GWAS_DIR}/${TAG}_tmp" "${GWAS_DIR}/${TAG}_spark"   # leftovers from the earlier run
else
export TMPDIR=${TMPBASE}/gwas
export SPARK_LOCAL_DIRS=${TMPDIR}/spark
mkdir -p "${SPARK_LOCAL_DIRS}"

python heig.py \
  --gwas \
  --out "${GWAS_DIR}/${TAG}" \
  --n-ldrs "${N_LDRS}" \
  --covar "${COVAR}" \
  --cat-covar-list "${CAT_COVARS}" \
  --ldrs "${LDRS}" \
  --geno-mt "${GENO_MT}" \
  --grch37 \
  --spark-conf "${SPARK_CONF}" \
  --loco-preds "${LOCO}" \
  --threads 4 \
  --not-save-genotype-data \
  --maf-min 0.01 \
  --maf-max 0.5

cleanup_tmp "${GWAS_DIR}/${TAG}_tmp" "${GWAS_DIR}/${TAG}_spark"
fi
[[ -d "${GWAS_DIR}/${TAG}.parquet" ]] || die "heig --gwas did not write ${GWAS_DIR}/${TAG}.parquet"

# ---- 5.2 Sumstats
python heig.py \
  --sumstats \
  --out "${SUMSTATS_OUT}" \
  --ldr-gwas-heig "${GWAS_DIR}/${TAG}.parquet"

if sumstats_done; then
  echo "sumstats written for ${TAG}; deleting ${GWAS_DIR}/${TAG}.parquet"
  cleanup_tmp "${GWAS_DIR}/${TAG}.parquet"
else
  die "sumstats failed for ${TAG}; keeping ${GWAS_DIR}/${TAG}.parquet"
fi
