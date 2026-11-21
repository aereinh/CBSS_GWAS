#!/bin/bash
# =====================================================================
# ABCD replication pipeline: shared settings.
# Every .slurm file sources this file. Change paths HERE, not in the .slurm files.
# =====================================================================

# ---- where this folder lives
PIPE_DIR="${PIPE_DIR:-$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)}"

# ---- software
HEIG_DIR=/path/to/heig
CONDA_ENV=heig
SPARK_CONF=/path/to/spark_config.json

# ---- roots
WORK=/path/to/UKB_PSC
OVERFLOW=/path/to/UKB_PSC_GWAS
LOG_ROOT=${WORK}/logs/abcd
TMP_ROOT=/path/to/tmp

# ---- inputs shared by all modalities
N_CLUSTERS=430
COORD_DIR=/path/to/UKB_PSC_GWAS/shared/coords
GENO_MT=${WORK}/output/genotype/genotype_abcd_qc.mt
CLUST_CSV=/path/to/CBSS_atlas/OtherModality/template/Fiber_atlas_100points_FullName.csv

# ---- covariates (used by make-ldr, relatedness and GWAS)
#   File columns: FID IID sex age age2 agesex age2sex PC1..PC40 site   (made by 00_prep_covar_site.py)
COVAR=${WORK}/covar/AgeGender_ABCD_withPCs_site.txt
CAT_COVARS=sex,site
COV_TAG=site      # tag of the covariate set; results go to output_ABCD_<MOD>_site

# ---------------------------------------------------------------------
# Per-modality settings.  MOD is the short name used for ABCD outputs.
#   RAW_DIR             ABCD *_Aligned.mat files + ABCD.txt (only needed for the format step)
#   UKB_BASES_PREFIX    UKB phase 1-6 fPCA bases:   <prefix>_<IDX>_bases_*
#   UKB_FPCA_LOG_PREFIX UKB phase 1-6 fPCA log:     <prefix>_<IDX>.log (starting number of LDRs)
#   OUTLIER_THRESH      |value| > threshold set to NA before read-image (empty = no cleaning)
#   REUSE_SEL_LDRS      if set, the select step copies this file instead of recomputing
# ---------------------------------------------------------------------
set_modality() {
  MOD="$1"
  local OM=/path/to/CBSS_atlas/OtherModality
  case "$MOD" in
    FA)
      UKB_MOD=FA
      RAW_DIR=${OM}/step2_FA_ABCD
      UKB_BASES_PREFIX=/path/to/UKB_PSC_GWAS/shared/LDR/FA/phase1to6/fpca/UKB_FA_phase16
      UKB_FPCA_LOG_PREFIX=${OVERFLOW}/FA/fpca/UKB_FA_phase16
      OUTLIER_THRESH=""
      # optional: reuse a previously selected number of LDRs per cluster
      REUSE_SEL_LDRS=${WORK}/output/output_ABCD_FA/selLDRs_val.txt
      ;;
    MD)
      UKB_MOD=dti_MD
      RAW_DIR=${OM}/OtherModality/dti_MD/step2_dti_MD_ABCD
      UKB_BASES_PREFIX=${WORK}/output/output_dti_MD/fpca/UKB_dti_MD_phase16
      UKB_FPCA_LOG_PREFIX=${WORK}/output/output_dti_MD/fpca/UKB_dti_MD_phase16
      OUTLIER_THRESH=100     # same threshold as the UKB dti_* runs
      REUSE_SEL_LDRS=""
      ;;
    AD)
      UKB_MOD=dti_L1
      RAW_DIR=${OM}/OtherModality/dti_L1/step2_dti_L1_ABCD
      UKB_BASES_PREFIX=${WORK}/output/output_dti_L1/fpca/UKB_dti_L1_phase16
      UKB_FPCA_LOG_PREFIX=${WORK}/output/output_dti_L1/fpca/UKB_dti_L1_phase16
      OUTLIER_THRESH=100
      REUSE_SEL_LDRS=""
      ;;
    # To add RD or MO, copy the MD block and swap dti_MD -> dti_RD / dti_MO.
    *)
      echo "ERROR: unknown modality '$MOD' (expected FA, MD or AD)" >&2
      return 1
      ;;
  esac

  # ABCD paths (derived; no need to edit)
  IMG_TXT_DIR=${WORK}/Modalities/ABCD_${MOD}                 # fibersClust_<IDX>.txt
  IMG_DIR=${WORK}/output/output_ABCD_${MOD}/images           # images don't depend on covariates, so shared across COV_TAGs
  RUN_NAME=ABCD_${MOD}${COV_TAG:+_${COV_TAG}}
  OUT_DIR=${WORK}/output/output_${RUN_NAME}
  LDR_DIR=${OUT_DIR}/ldrs
  REL_DIR=${OUT_DIR}/rel
  GWAS_DIR=${OUT_DIR}/gwas
  SEL_LDRS=${OUT_DIR}/selLDRs_val.txt
  SUMSTATS_DIR=${OVERFLOW}/${RUN_NAME}/sumstats
}

# ---- helpers
load_heig_env() {
  module purge
  module load anaconda
  conda activate "${CONDA_ENV}"
  cd "${HEIG_DIR}"
}

# selected number of LDRs for cluster $1 (row i+1, column 2 of selLDRs_val.txt)
n_ldrs_val() {
  local n
  n=$(awk -F',' -v i="$1" 'NR == i+1 {print $2}' "${SEL_LDRS}")
  if [[ -z "$n" || "$n" == "NA" ]]; then
    echo "ERROR: no selected number of LDRs for cluster $1 in ${SEL_LDRS}" >&2
    return 1
  fi
  echo "$n"
}

die() { echo "ERROR: $*" >&2; exit 1; }
