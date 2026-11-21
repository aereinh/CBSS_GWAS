#!/bin/bash
#SBATCH --job-name=abcd_rep_agg
#SBATCH --time=00:30:00
#SBATCH --partition=general
#SBATCH --ntasks=1
#SBATCH --cpus-per-task=1
#SBATCH --mem=8G

# Step 4: merge tracts into loci and write every number the manuscript reports for ABCD.
source "${REP_SCRIPTS:?run through submit_rep.sh}/rep_config.sh"
set -eo pipefail
module purge; module load anaconda; conda activate "${CONDA_ENV}"
python "${REP_SCRIPTS}/4_aggregate_loci.py"
