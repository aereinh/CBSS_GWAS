#!/bin/bash
# =============================================================================
# step1_run.slurm
#
# Array job: find the best voxel per locus for colocalization.
# Array range should cover the number of loci in your loci-info CSV.
# =============================================================================
#SBATCH --job-name=coloc_sel_vox
#SBATCH --output=logs/coloc_sel_vox_%a.out
#SBATCH --error=logs/coloc_sel_vox_%a.err
#SBATCH --partition=general
#SBATCH --time=01:00:00
#SBATCH --mem=8G
#SBATCH --cpus-per-task=1
#SBATCH --array=1-315%50      # adjust upper bound to number of loci

# ── User-configurable paths ───────────────────────────────────────────────────
MODALITY="dti_FA"
VOX_DIR="/path/to/output/post_screen/disc/fibercluster"
LOCI_CSV="/path/to/loci/ALL_TRACTS_locus_summary.csv"
OUT_DIR="/path/to/coloc/sel_voxels"
SCRIPTS_DIR="/path/to/scripts/2_downstream_analyses/coloc"
# ─────────────────────────────────────────────────────────────────────────────

module purge
module load anaconda
conda activate heig

mkdir -p "${OUT_DIR}" logs

python "${SCRIPTS_DIR}/step1_select_best_voxel.py" \
    --modality  "${MODALITY}" \
    --locus-id  "${SLURM_ARRAY_TASK_ID}" \
    --vox-dir   "${VOX_DIR}" \
    --loci-info "${LOCI_CSV}" \
    --out-dir   "${OUT_DIR}"
