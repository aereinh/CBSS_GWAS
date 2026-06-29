#!/bin/bash
# =============================================================================
# genetic_correlation.slurm
#
# Estimate genetic correlations between white-matter voxel phenotypes and
# external GWAS traits using HEIG (--heri-gc with --y2-sumstats).
#
# Traits included: ADHD, Bipolar Disorder, MDD, Cognitive Performance (CP),
# Educational Attainment (EA), Insomnia, ASD, Schizophrenia (SCZ),
# Alzheimer's Disease, Intelligence, Depressed Affect, Depression,
# Neuroticism, Worry.
#
# Add --overlap flag for traits with UK Biobank sample overlap.
#
# Array index = fiber cluster ID (1-indexed).
# =============================================================================
#SBATCH --job-name=gen_corr
#SBATCH --output=logs/gc_%a.out
#SBATCH --error=logs/gc_%a.err
#SBATCH --mem=25G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --time=3:00:00
#SBATCH --array=1-430%80

# ── User-configurable paths ────────────────────────────────────────────────
MODALITY="dti_FA"
HEIG_DIR="/path/to/heig"
SUMSTATS_DIR="/path/to/output/sumstats"
FPCA_DIR="/path/to/output/fpca"
LDR_DIR="/path/to/output/ldrs"
GC_OUT="/path/to/output/genetic_correlation"
NLDR_FILE="/path/to/output/selLDRs.txt"

# LD matrices
LD_INV="/path/to/ldmatrix/ukb_ld_inv_chr{1:22}"
LD="/path/to/ldmatrix/ukb_ld_chr{1:22}"

# Directory containing external GWAS summary statistics
EXT_SUMSTATS_DIR="/path/to/external_sumstats"
# ──────────────────────────────────────────────────────────────────────────

module purge
module load anaconda
conda activate heig

mkdir -p "${GC_OUT}"

IDX=${SLURM_ARRAY_TASK_ID}
N_LDRS=$(awk -F',' -v i=${IDX} 'NR == i+1 {print $2}' "${NLDR_FILE}")

# Helper function to run genetic correlation for one external trait
run_gc() {
    local TRAIT_NAME="$1"
    local SUMSTATS_FILE="$2"
    local OVERLAP_FLAG="$3"   # "--overlap" or ""

    local OUT_DIR="${GC_OUT}/${TRAIT_NAME}"
    mkdir -p "${OUT_DIR}"

    python "${HEIG_DIR}/heig.py" \
        --heri-gc \
        --out          "${OUT_DIR}/UKB_${MODALITY}_allphases_${IDX}_${TRAIT_NAME}" \
        --n-ldrs       "${N_LDRS}" \
        --ldr-sumstats "${SUMSTATS_DIR}/UKB_${MODALITY}_allphases_${IDX}" \
        --bases        "${FPCA_DIR}/UKB_${MODALITY}_allphases_${IDX}_bases_top"* \
        --ldr-cov      "${LDR_DIR}/UKB_${MODALITY}_allphases_${IDX}_ldr_cov_top"* \
        --threads      4 \
        --ld-inv       "${LD_INV}" \
        --ld           "${LD}" \
        --y2-sumstats  "${EXT_SUMSTATS_DIR}/${SUMSTATS_FILE}" \
        ${OVERLAP_FLAG}
}

# ── Psychiatric / cognitive traits ───────────────────────────────────────────
# Add --overlap for traits with UK Biobank sample overlap
run_gc "ADHD"             "ADHD2022_iPSYCH_deCODE_PGC.meta"                        "--overlap"
run_gc "bipolar"          "daner_bip_pgc3_nm_noukbiobank"                           ""
run_gc "MDD"              "daner_pgc_mdd_meta_w2_no23andMe_rmUKBB"                  ""
run_gc "CP"               "GWAS_CP_all"                                              "--overlap"
run_gc "EA"               "GWAS_EA_excl23andMe"                                      "--overlap"
run_gc "insomnia"         "insomnia_ukb2b_EUR_sumstats_with_chrX_mac_100"            "--overlap"
run_gc "ASD"              "iPSYCH-PGC_ASD_Nov2017"                                  ""
run_gc "SCZ"              "PGC3_SCZ_wave3.european.autosome.public.v3.vcf"           ""
run_gc "ALZ"              "PGCALZ2sumstatsExcluding23andMe"                          "--overlap"
run_gc "intelligence"     "SavageJansen_2018_intelligence_metaanalysis"              "--overlap"
run_gc "depressed_affect" "sumstats_depressed_affect_ctg_format"                    "--overlap"
run_gc "depression"       "sumstats_depression_ctg_format"                           "--overlap"
run_gc "neuroticism"      "sumstats_neuroticism_ctg_format"                          "--overlap"
run_gc "worry"            "sumstats_worry_ctg_format"                                "--overlap"
