#!/bin/bash
#SBATCH --job-name=abcd_rep_tract
#SBATCH --time=01:00:00
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=25G
#SBATCH --array=0-16

# Steps 2-3, one tract per task:
#   2: locus-level ("genetic") replication within +-250 kb of each discovery lead
#   3: spatial overlap of the discovery footprint
source "${REP_SCRIPTS:?run through submit_rep.sh}/rep_config.sh"
set -eo pipefail
TRACT=${TRACTS[$SLURM_ARRAY_TASK_ID]}
echo "Tract: ${TRACT}"
module purge; module load anaconda; conda activate "${CONDA_ENV}"

python "${REP_SCRIPTS}/2_replication_one_tract.py"     --tract "${TRACT}"
python "${REP_SCRIPTS}/3_spatial_overlap_one_tract.py" --tract "${TRACT}"
