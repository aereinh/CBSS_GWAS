#!/bin/bash
# =============================================================================
# step6_run.slurm
#
# For each fiber cluster, look up its assigned tract and run targeted
# validation voxel GWAS (step6_run_vgwas_tract.sh).
# Clusters not assigned to a tract are skipped.
#
# Array index = fiber cluster ID (1-indexed).
# =============================================================================
#SBATCH --job-name=vgwas_rep
#SBATCH --output=logs/step6_%a.out
#SBATCH --error=logs/step6_%a.err
#SBATCH --time=00:30:00
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=4
#SBATCH --mem=16G
#SBATCH --array=1-430%80

# ── User-configurable ─────────────────────────────────────────────────────
MODALITY="dti_FA"
NLDR_FILE_VAL="/path/to/output/selLDRs_val.txt"
CLUSTER_TO_TRACT="/path/to/cluster_tract_summary.csv"   # columns: cluster, tract_assign
SCRIPTS_DIR="/path/to/scripts/4_post_screening_replication"
# ──────────────────────────────────────────────────────────────────────────

module purge
module load gcc/11.2.0
module load anaconda
conda activate heig
export LD_LIBRARY_PATH="${CONDA_PREFIX}/lib:${LD_LIBRARY_PATH}"

CLUSTER_ID=${SLURM_ARRAY_TASK_ID}
N_LDRS=$(awk -F',' -v i=${CLUSTER_ID} 'NR == i+1 {print $2}' "${NLDR_FILE_VAL}")

# Resolve tract name for this cluster
TRACT=$(awk -F',' -v id="${CLUSTER_ID}" '
    NR>1 {
        gsub(/"/, "", $0)
        if (($2+0) == id) { print $NF; exit }
    }
' "${CLUSTER_TO_TRACT}")

if [ -z "${TRACT}" ] || [ "${TRACT}" = "NA" ]; then
    echo "No assigned tract for cluster ${CLUSTER_ID} — skipping."
    exit 0
fi

echo "Cluster ${CLUSTER_ID} → Tract: ${TRACT}"

bash "${SCRIPTS_DIR}/step6_run_vgwas_tract.sh" \
    --tract    "${TRACT}" \
    --cluster  "${CLUSTER_ID}" \
    --n-ldrs   "${N_LDRS}" \
    --modality "${MODALITY}"
