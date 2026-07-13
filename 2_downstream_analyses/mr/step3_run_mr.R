#!/usr/bin/env Rscript
# =============================================================================
# step3_run_mr.R
#
# Run two-sample Mendelian Randomization using TwoSampleMR, with white matter
# fiber cluster LDR sumstats as the exposure and psychiatric/cognitive traits
# from OpenGWAS as outcomes.
#
# Usage:
#   Rscript step3_run_mr.R <cluster_id>
#
# Requirements:
#   - R packages: TwoSampleMR, data.table, dplyr
#   - An OpenGWAS JWT token (obtain at https://api.opengwas.io/api/#/auth/get_token)
#     Set as environment variable OPENGWAS_JWT, or edit the jwt variable below.
#
# Outputs:
#   One TSV file per cluster: UKB_{MODALITY}_allphases_{cluster_id}_MR_results.tsv
# =============================================================================

suppressPackageStartupMessages({
  library(TwoSampleMR)
  library(data.table)
  library(dplyr)
})

# ── User-configurable paths ───────────────────────────────────────────────────
MODALITY     <- "dti_FA"
SCREENED_DIR <- "/path/to/output/mr/screened"     # step 1 output directory
OUT_DIR      <- "/path/to/output/mr/results"      # MR results output directory

# OpenGWAS JWT — obtain from https://api.opengwas.io/api/#/auth/get_token
# Recommended: set the OPENGWAS_JWT environment variable rather than hard-coding.
jwt <- Sys.getenv("OPENGWAS_JWT", unset = "")
if (nchar(jwt) == 0) {
  stop("OPENGWAS_JWT environment variable not set. See https://api.opengwas.io/api/#/auth/get_token")
}

# 14 psychiatric/cognitive outcome traits from OpenGWAS
OUTCOMES <- c(
  "ieu-a-1183",          # ADHD
  "ieu-b-5111",          # Bipolar disorder
  "ieu-a-1187",          # MDD
  "finn-b-G6_CP",        # Cognitive performance
  "ebi-a-GCST90029012",  # Educational attainment
  "ukb-b-3957",          # Insomnia
  "ieu-a-1185",          # ASD
  "ieu-b-5099",          # Schizophrenia
  "ebi-a-GCST90027158",  # Alzheimer's disease (2022)
  "ebi-a-GCST006250",    # Intelligence
  "ebi-a-GCST006475",    # Depressed affect
  "ukb-b-12064",         # Depression
  "ukb-b-4630",          # Neuroticism
  "ukb-b-13653"          # Worry
)

MR_METHODS <- c(
  "mr_ivw_fe", "mr_ivw_mre", "mr_ivw", "mr_ivw_radial",
  "mr_egger_regression",
  "mr_simple_median", "mr_weighted_median", "mr_penalised_weighted_median",
  "mr_simple_mode", "mr_simple_mode_nome", "mr_weighted_mode", "mr_weighted_mode_nome",
  "mr_grip"
)
# ─────────────────────────────────────────────────────────────────────────────

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 1) {
  stop("Usage: Rscript step3_run_mr.R <cluster_id>")
}
cluster_id <- args[1]

dir.create(OUT_DIR, recursive = TRUE, showWarnings = FALSE)
out_file <- file.path(OUT_DIR, paste0("UKB_", MODALITY, "_allphases_", cluster_id, "_MR_results.tsv"))

# ── Step 1: Load exposure (DTI fiber cluster) sumstats ───────────────────────
sumstats_file <- file.path(SCREENED_DIR, paste0("UKB_", MODALITY, "_allphases_", cluster_id, ".txt"))
clump_file    <- file.path(SCREENED_DIR, paste0("UKB_", MODALITY, "_allphases_", cluster_id, "_clump_r2_0.01.clumped"))

if (!file.exists(sumstats_file)) {
  cat("No screened sumstats found:", sumstats_file, "— skipping.\n")
  quit(save = "no", status = 0)
}
if (!file.exists(clump_file)) {
  cat("No clump file found:", clump_file, "— skipping.\n")
  quit(save = "no", status = 0)
}

dti_sumstats <- fread(sumstats_file)
dti_clump    <- fread(clump_file)
dti_sumstats <- dti_sumstats %>% filter(SNP %in% unique(dti_clump$SNP))
dti_sumstats$SE <- dti_sumstats$BETA / dti_sumstats$Z
dti_sumstats <- as.data.frame(dti_sumstats)

if (nrow(dti_sumstats) == 0) {
  cat("No significant clumped SNPs for cluster", cluster_id, "— skipping.\n")
  quit(save = "no", status = 0)
}

dti_exp <- format_data(
  dti_sumstats,
  type             = "exposure",
  snp_col          = "SNP",
  beta_col         = "BETA",
  se_col           = "SE",
  effect_allele_col = "A1",
  other_allele_col = "A2",
  pval_col         = "P",
  samplesize_col   = "N",
  z_col            = "Z",
  chr_col          = "CHR",
  pos_col          = "POS",
  phenotype_col    = "INDEX"
)

# ── Step 2: Extract outcome data from OpenGWAS ───────────────────────────────
cat("Querying OpenGWAS for", length(OUTCOMES), "outcomes ...\n")
outcome_data <- extract_outcome_data(
  snps          = unique(dti_sumstats$SNP),
  outcomes      = OUTCOMES,
  opengwas_jwt  = jwt
)

if (is.null(outcome_data)) {
  cat("No matching SNPs found in outcome sumstats. Exiting.\n")
  quit(save = "no", status = 0)
}

# ── Step 3: Harmonise ────────────────────────────────────────────────────────
dat <- harmonise_data(dti_exp, outcome_data)

# ── Step 4: Run MR ───────────────────────────────────────────────────────────
res <- mr(dat, method_list = MR_METHODS)

# ── Save ─────────────────────────────────────────────────────────────────────
fwrite(res, file = out_file, sep = "\t")
cat("MR results written to:", out_file, "\n")
