#!/bin/bash
# =============================================================================
# step2_run.slurm
#
# Spatial closure with 10,000 rotation draws, one array task per modality
# (FA, MD = dti_MD, AD = dti_L1).  Run step3_merge_and_fdr.py afterwards.
# =============================================================================
#SBATCH --job-name=closure
#SBATCH --partition=general
#SBATCH --time=12:00:00
#SBATCH --mem=30G
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --array=0-2
#SBATCH --output=logs/closure_%A_%a.out
#SBATCH --error=logs/closure_%A_%a.err

# ── User-configurable paths ────────────────────────────────────────────────
SCRIPTS_DIR="/path/to/scripts/8_druggability_spatial_closure"
COLOC_TABLE="/path/to/coloc/coloc_hits_strong_all_pairs.tsv"   # colocalized locus x trait pairs
OUT_DIR="/path/to/output/closure"
# ──────────────────────────────────────────────────────────────────────────

module load anaconda
conda activate heig
mkdir -p "${OUT_DIR}" logs

MODS=(FA dti_MD dti_L1)
M=${MODS[$SLURM_ARRAY_TASK_ID]}

python "${SCRIPTS_DIR}/step2_spatial_closure.py" \
    --coloc-table "${COLOC_TABLE}" \
    --modalities "$M" \
    --n-rot 10000 \
    --seed 1 \
    --out "${OUT_DIR}/closure_${M}.csv"
