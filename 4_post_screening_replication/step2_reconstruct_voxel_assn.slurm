#!/bin/bash
# =============================================================================
# step2_reconstruct_voxel_assn.slurm
#
# For each fiber cluster, reconstruct full voxel-level association statistics
# for the subset of significant SNPs identified in step 1, using a relaxed
# per-voxel threshold (1e-5).  This produces the data needed for cluster-size
# screening in step 3.
#
# Array index = fiber cluster ID (1-indexed).
# =============================================================================
#SBATCH --job-name=reconstruct
#SBATCH --output=logs/reconstruct_%a.out
#SBATCH --error=logs/reconstruct_%a.err
#SBATCH --time=00:30:00
#SBATCH --mem=25G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --array=1-430%50

# ── User-configurable paths ────────────────────────────────────────────────
MODALITY="dti_FA"
HEIG_DIR="/path/to/heig"
SUMSTATS_DIR="/path/to/output/sumstats"
FPCA_DIR="/path/to/output/fpca"
LDR_DIR="/path/to/output/ldrs"
MAIN_DIR="/path/to/output"
NLDR_FILE="/path/to/output/selLDRs.txt"
# ──────────────────────────────────────────────────────────────────────────

module purge
module load anaconda
conda activate heig

IDX=${SLURM_ARRAY_TASK_ID}
N_LDRS=$(awk -F',' -v i=${IDX} 'NR == i+1 {print $2}' "${NLDR_FILE}")

SNP_LIST="${MAIN_DIR}/output_${MODALITY}/post_screen/disc/fibercluster/UKB_${MODALITY}_allphases_${IDX}_sig_unique_snps.txt"
OUT_FILE="${MAIN_DIR}/output_${MODALITY}/post_screen/disc/fibercluster/UKB_${MODALITY}_allphases_${IDX}_vgwas_1e5"

if [ ! -f "${SNP_LIST}" ]; then
    echo "No significant SNP list for cluster ${IDX}. Skipping."
    exit 0
fi

python "${HEIG_DIR}/heig.py" \
    --voxel-gwas \
    --out          "${OUT_FILE}" \
    --n-ldrs       "${N_LDRS}" \
    --ldr-sumstats "${SUMSTATS_DIR}/UKB_${MODALITY}_allphases_${IDX}" \
    --bases        "${FPCA_DIR}/UKB_${MODALITY}_allphases_${IDX}_bases_top"* \
    --ldr-cov      "${LDR_DIR}/UKB_${MODALITY}_allphases_${IDX}_ldr_cov_top"* \
    --extract      "${SNP_LIST}" \
    --threads      4 \
    --sig-thresh   1e-05
