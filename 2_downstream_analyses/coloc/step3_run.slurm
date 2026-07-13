#!/bin/bash
# =============================================================================
# step3_run.slurm
#
# Array job: merge imaging voxel GWAS sumstats with trait sumstats for coloc.
# Requires a pre-built filelist of step 2 output files (one per line).
#
# Build filelist before submitting:
#   ls /path/to/coloc/voxelgwas/UKB_dti_FA_LOC*_CLST*_VOX*.txt \
#     > /path/to/coloc/tmp/voxelgwas_files.txt
# Then submit:
#   sbatch step3_run.slurm
# =============================================================================
#SBATCH --job-name=coloc_merge
#SBATCH --output=logs/coloc_merge_%a.out
#SBATCH --error=logs/coloc_merge_%a.err
#SBATCH --partition=general
#SBATCH --time=02:00:00
#SBATCH --mem=32G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-315%80      # adjust upper bound to number of loci

# ── User-configurable paths ───────────────────────────────────────────────────
FILELIST="/path/to/coloc/tmp/voxelgwas_files.txt"   # one voxelgwas .txt file per line
TRAIT_DIR="/path/to/coloc/trait_sumstats"            # *_coloc.txt files
OUT_DIR="/path/to/coloc/coloc_sumstats"
SCRIPTS_DIR="/path/to/scripts/2_downstream_analyses/coloc"
# ─────────────────────────────────────────────────────────────────────────────

module purge
module load anaconda
conda activate heig

mkdir -p "${OUT_DIR}" logs

NFILES=$(wc -l < "${FILELIST}")
if (( SLURM_ARRAY_TASK_ID > NFILES )); then
    echo "Task ${SLURM_ARRAY_TASK_ID} > NFILES ${NFILES}. Nothing to do."
    exit 0
fi

SUMSTATS_FP=$(sed -n "${SLURM_ARRAY_TASK_ID}p" "${FILELIST}")
echo "Processing: ${SUMSTATS_FP}"

python "${SCRIPTS_DIR}/step3_merge_sumstats.py" \
    --vox-sumstats "${SUMSTATS_FP}" \
    --trait-dir    "${TRAIT_DIR}" \
    --out-dir      "${OUT_DIR}"

echo "Done: $(date)"
