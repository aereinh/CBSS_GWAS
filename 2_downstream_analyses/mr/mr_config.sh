#!/bin/bash
# =====================================================================
# Two-sample MR with local (pre-downloaded) outcome summary statistics.
# Shared settings: every .slurm / .sh file here sources this file, and the
# R scripts read the same values from the environment.
# =====================================================================
export WORK=/path/to/UKB_PSC

# ---- inputs
export RAW_DIR=${WORK}/trait_sumstats/raw            # put the downloaded trait files HERE (unzipped or .gz)
export SCREEN_DIR=${WORK}/ldr-sumstats/screened      # step 1-2 outputs (screened LDR sumstats + clump files)
export EXP_PREFIX=UKB_dti_FA_phase16_                # prefix of the step 1-2 files (UKB_<MODALITY>_phase16_)
# (no reference panel is needed for the outcomes: allele frequencies come from the trait files,
#  see eaf_from in trait_config.R)
export N_FIBERS=430

# ---- outputs
export HARM_DIR=${WORK}/trait_sumstats/harmonized    # instrument list + one small harmonized file per trait
export OUT_ROOT=${WORK}/mr_v2                         # MR results: ${OUT_ROOT}/<trait>/...
export LOG_DIR=${WORK}/logs/mr_v2

# ---- traits (same 14 as the genetic-correlation panel; names must match trait_config.R)
TRAITS=(ADHD bip mdd CP EA insomnia ASD SCZ ALZ2 intelligence depressed_affect depression neuroticism worry)

# ---- where these scripts live
export MR_DIR="${MR_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

load_r()     { module purge; module load r; }
load_plink() { module load plink/1.90b6.21; }
