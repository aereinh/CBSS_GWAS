#!/bin/bash
# =============================================================================
# aggregate_gc_tracts.slurm
#
# Tract-level and network-pair-level aggregation of the location-level genetic
# correlations: 3 diffusion measures x 14 traits x (17 tracts | 66 network pairs).
# Paths are set in PipelineConfig at the top of aggregate_gc_tracts.py.
# =============================================================================
#SBATCH --job-name=gc_tract_agg
#SBATCH --output=logs/gc_tract_agg_%j.out
#SBATCH --error=logs/gc_tract_agg_%j.err
#SBATCH --time=08:00:00
#SBATCH --mem=30G
#SBATCH --cpus-per-task=4
#SBATCH --ntasks=1

# ── User-configurable paths ────────────────────────────────────────────────
SCRIPTS_DIR="/path/to/scripts/2_downstream_analyses/herigc_postprocess"
OUT_DIR="/path/to/UKB_PSC/result_data/agg"
# ──────────────────────────────────────────────────────────────────────────

set -euo pipefail
module load anaconda
conda activate heig
mkdir -p "${OUT_DIR}" logs

# Tracts (cluster_tract_summary.csv column tract_assign)
python -u "${SCRIPTS_DIR}/aggregate_gc_tracts.py" \
    --out "${OUT_DIR}/tract_assign_gc_results.csv"

# Functional network pairs (column network_pair)
python -u "${SCRIPTS_DIR}/aggregate_gc_tracts.py" \
    --group-col network_pair \
    --out "${OUT_DIR}/network_pair_gc_results.csv"
