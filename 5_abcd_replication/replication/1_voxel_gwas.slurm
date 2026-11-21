#!/bin/bash
#SBATCH --job-name=abcd_rep_vgwas
#SBATCH --time=00:30:00
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1-430%80

# Step 1: reconstruct ABCD location-level associations for the discovery candidate
# SNPs of the cluster's tract and keep those with P < SIG_THRESH.
#   -> REP_DIR/ABCD_<MOD>_<cluster>_vgwas_<TAG>.txt
source "${REP_SCRIPTS:?run through submit_rep.sh}/rep_config.sh"
set -eo pipefail

CLUSTER_ID=${SLURM_ARRAY_TASK_ID}
OUT=${REP_DIR}/ABCD_${MOD}_${CLUSTER_ID}_vgwas_${TAG}
DONE=${OUT}.done

if [[ -z "${OVERWRITE:-}" && -e "${DONE}" ]]; then echo "[SKIP] cluster ${CLUSTER_ID} already done"; exit 0; fi

TRACT=$(awk -F',' -v id="$CLUSTER_ID" 'NR>1{gsub(/"/,""); if(($2+0)==id){print $NF; exit}}' "${CLUSTER_TRACT_CSV}")
if [[ -z "$TRACT" || "$TRACT" == "NA" ]]; then
  echo "No assigned tract for cluster ${CLUSTER_ID}: skipping."; mkdir -p "${REP_DIR}"; echo "no tract" > "${DONE}"; exit 0
fi
N_LDRS=$(n_ldrs_val "${CLUSTER_ID}")
EXTRACT=${DISC_TRACT_DIR}/${DISC_PREFIX}_${TRACT}_candidate_sig_snps.txt
[[ -s "${EXTRACT}" ]] || die "missing discovery candidate SNP list ${EXTRACT}"
[[ -s "${SUMSTATS_DIR}/ABCD_${MOD}_${CLUSTER_ID}.snpinfo" ]] || die "missing ABCD LDR sumstats for cluster ${CLUSTER_ID} in ${SUMSTATS_DIR}"
compgen -G "${LDR_DIR}/ABCD_${MOD}_${CLUSTER_ID}_ldr_cov_top*" > /dev/null || die "missing ${LDR_DIR}/ABCD_${MOD}_${CLUSTER_ID}_ldr_cov_top*"
echo "Cluster ${CLUSTER_ID} -> tract ${TRACT}; n_ldrs=${N_LDRS}; sumstats=${SUMSTATS_DIR}"

module purge
module load gcc/11.2.0
module load anaconda
conda activate "${CONDA_ENV}"
export LD_LIBRARY_PATH="$CONDA_PREFIX/lib:$LD_LIBRARY_PATH"
export TMP_DIR=${TMP_ROOT}
mkdir -p "${REP_DIR}"

cd "$(dirname "${VGWAS_HEIG}")"
python "${VGWAS_HEIG}" \
  --voxel-gwas \
  --out "${OUT}" \
  --n-ldrs "${N_LDRS}" \
  --ldr-sumstats "${SUMSTATS_DIR}/ABCD_${MOD}_${CLUSTER_ID}" \
  --bases ${UKB_BASES_PREFIX}_${CLUSTER_ID}_bases_top* \
  --ldr-cov "${LDR_DIR}/ABCD_${MOD}_${CLUSTER_ID}"_ldr_cov_top* \
  --extract "${EXTRACT}" \
  --threads 4 \
  --sig-thresh "${SIG_THRESH}"

# HEIG writes no table when nothing passes the threshold; the .done file marks the task as finished either way
if [[ -s "${OUT}.txt" ]]; then echo "rows: $(($(wc -l < "${OUT}.txt") - 1))" > "${DONE}"; else echo "no association below ${SIG_THRESH}" > "${DONE}"; fi
cat "${DONE}"
