#!/bin/bash
#SBATCH --job-name=mr_instruments
#SBATCH --mem=16G
#SBATCH --cpus-per-task=1
#SBATCH --partition=general
#SBATCH --time=02:00:00

# Step 3b: union of all clumped instruments (no reference panel needed: allele frequencies come from
# the trait files; for Alzheimer's from the Nagel neuroticism file, see eaf_from in trait_config.R).
# Submitted by submit_mr.sh prep
source "${MR_DIR:?submit with submit_mr.sh}/mr_config.sh"
set -eo pipefail
mkdir -p "${HARM_DIR}"

load_r
Rscript "${MR_DIR}/step3b_build_instruments.R"

[[ -s "${HARM_DIR}/instruments.tsv" ]] || { echo "ERROR: no instruments written" >&2; exit 1; }
