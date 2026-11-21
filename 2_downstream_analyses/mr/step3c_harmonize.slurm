#!/bin/bash
#SBATCH --job-name=mr_harmonize
#SBATCH --mem=48G
#SBATCH --cpus-per-task=2
#SBATCH --partition=general
#SBATCH --time=01:00:00

# Step 3c: one array task per trait (index into TRAITS in mr_config.sh).
# Submitted by submit_mr.sh prep
source "${MR_DIR:?submit with submit_mr.sh}/mr_config.sh"
set -eo pipefail
TRAIT=${TRAITS[$((SLURM_ARRAY_TASK_ID - 1))]}
echo "trait: ${TRAIT}"
load_r
Rscript "${MR_DIR}/step3c_harmonize_trait.R" "${TRAIT}"
