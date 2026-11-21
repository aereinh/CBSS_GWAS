# =====================================================================
# One entry per trait: where the raw file is and what its columns are called.
# Column names below were checked against the headers of the downloaded files
# (all 14 traits). If a file is replaced, rerun `Rscript step3a_peek_headers.R` and recheck.
# step3c refuses traits with verified = FALSE.
#
# Fields
#   file         start of the file name in RAW_DIR (any extension: .txt, .tsv, .gz); must match one file
#   label        name written to the results
#   snp          rsID column (NULL if the file has none)
#   chr, pos     chromosome and base-pair position columns (GRCh37/hg19, like the LDR GWAS)
#   ea, oa       effect allele and other (non-effect) allele columns
#   effect       effect column; effect_type says what it is: "beta", "or" (odds ratio) or "z"
#   se           standard-error column (of beta / log OR); NULL if absent
#   p            P-value column
#   n            per-SNP sample-size column, or NULL
#   ncase, ncontrol   per-SNP case / control count columns, or NULL
#   n_const      total sample size to use when the file has no N column
#   eaf          effect-allele frequency column (a regular expression is allowed, e.g. "^FRQ_U_"), or NULL
#   eaf_from     only when eaf = NULL and effect_type = "z": name of another trait in this list whose file
#                supplies the allele frequency (matched by chr:pos + alleles, flipped to this trait's effect allele)
#   match_by     "pos" (chr:pos + alleles; default) or "rsid"
#   exclude      regions to remove from the outcome, "chr:start-end" (GRCh37); NULL for none
#   verified     TRUE once the column names have been checked AND which allele is the effect allele
# =====================================================================
APOE_REGION <- "19:44400000-46500000"   # GRCh37; removed from the Alzheimer's outcome

TRAIT_CONFIG <- list(

  ADHD = list(   # Demontis 2023 (daner format: A1 = effect allele, OR)
    file = "ADHD2022_iPSYCH_deCODE_PGC.meta", label = "ADHD",
    snp = "SNP", chr = "CHR", pos = "BP", ea = "A1", oa = "A2",
    effect = "OR", effect_type = "or", se = "SE", p = "P",
    n = NULL, ncase = "Nca", ncontrol = "Nco", n_const = 225534, eaf = "^FRQ_U_",
    match_by = "pos", exclude = NULL, verified = TRUE),

  bip = list(    # Mullins 2021, no UKB (daner format)
    file = "daner_bip_pgc3_nm_noukbiobank", label = "Bipolar disorder",
    snp = "SNP", chr = "CHR", pos = "BP", ea = "A1", oa = "A2",
    effect = "OR", effect_type = "or", se = "SE", p = "P",
    n = NULL, ncase = "Nca", ncontrol = "Nco", n_const = 413466, eaf = "^FRQ_U_",
    match_by = "pos", exclude = NULL, verified = TRUE),

  mdd = list(    # Wray 2018, no 23andMe, no UKB (daner format)
    file = "daner_pgc_mdd_meta_w2_no23andMe_rmUKBB", label = "Major depressive disorder",
    snp = "SNP", chr = "CHR", pos = "BP", ea = "A1", oa = "A2",
    effect = "OR", effect_type = "or", se = "SE", p = "P",
    n = NULL, ncase = "Nca", ncontrol = "Nco", n_const = 480359, eaf = "^FRQ_U_",
    match_by = "pos", exclude = NULL, verified = TRUE),

  CP = list(     # Lee 2018 cognitive performance (SSGAC: A1 = effect allele)
    file = "GWAS_CP_all", label = "Cognitive performance",
    snp = "MarkerName", chr = "CHR", pos = "POS", ea = "A1", oa = "A2",
    effect = "Beta", effect_type = "beta", se = "SE", p = "Pval",
    n = NULL, ncase = NULL, ncontrol = NULL, n_const = 257828, eaf = "EAF",
    match_by = "pos", exclude = NULL, verified = TRUE),

  EA = list(     # Lee 2018 educational attainment, excluding 23andMe
    file = "GWAS_EA_excl23andMe", label = "Educational attainment",
    snp = "MarkerName", chr = "CHR", pos = "POS", ea = "A1", oa = "A2",
    effect = "Beta", effect_type = "beta", se = "SE", p = "Pval",
    n = NULL, ncase = NULL, ncontrol = NULL, n_const = 766345, eaf = "EAF",
    match_by = "pos", exclude = NULL, verified = TRUE),

  insomnia = list(   # Watanabe 2022 (PLINK output: A1 = tested/effect allele; SNP is chr:pos:a1_a2; MAF is minor-allele freq, not EAF)
    file = "insomnia_ukb2b_EUR_sumstats_20190311_with_chrX_mac_100", label = "Insomnia",
    snp = "RSID_UKB", chr = "CHR", pos = "BP", ea = "A1", oa = "A2",
    effect = "OR", effect_type = "or", se = "SE", p = "P",
    n = "NMISS", ncase = NULL, ncontrol = NULL, n_const = 386988, eaf = NULL,
    match_by = "pos", exclude = NULL, verified = TRUE),

  ASD = list(    # Grove 2019 (no frequency or N columns)
    file = "iPSYCH-PGC_ASD_Nov2017", label = "Autism spectrum disorder",
    snp = "SNP", chr = "CHR", pos = "BP", ea = "A1", oa = "A2",
    effect = "OR", effect_type = "or", se = "SE", p = "P",
    n = NULL, ncase = NULL, ncontrol = NULL, n_const = 46350, eaf = NULL,
    match_by = "pos", exclude = NULL, verified = TRUE),

  SCZ = list(    # Trubetskoy 2022 (lines starting with ## are skipped automatically)
    file = "PGC3_SCZ_wave3.european.autosome.public.v3", label = "Schizophrenia",
    snp = "ID", chr = "CHROM", pos = "POS", ea = "A1", oa = "A2",
    effect = "BETA", effect_type = "beta", se = "SE", p = "PVAL",
    n = NULL, ncase = "NCAS", ncontrol = "NCON", n_const = 161405, eaf = "FCON",
    match_by = "pos", exclude = NULL, verified = TRUE),

  ALZ2 = list(   # Wightman 2021 excl. 23andMe. Checked against the real file: no rsID, no beta/SE, no EAF.
    file = "PGCALZ2sumstatsExcluding23andMe", label = "Alzheimer's disease",
    snp = NULL, chr = "chr", pos = "PosGRCh37", ea = "testedAllele", oa = "otherAllele",
    effect = "z", effect_type = "z", se = NULL, p = "p",
    n = "N", ncase = NULL, ncontrol = NULL, n_const = NULL, eaf = NULL,
    eaf_from = "neuroticism",   # no EAF in the file: take it from Nagel 2018 (EAF_UKB, UK Biobank Europeans), aligned by position + alleles
    match_by = "pos", exclude = APOE_REGION, verified = TRUE),

  intelligence = list(   # Savage 2018
    file = "SavageJansen_2018_intelligence_metaanalysis", label = "Intelligence",
    snp = "SNP", chr = "CHR", pos = "POS", ea = "A1", oa = "A2",
    effect = "stdBeta", effect_type = "beta", se = "SE", p = "P",
    n = "N_analyzed", ncase = NULL, ncontrol = NULL, n_const = 269867, eaf = "EAF_HRC",
    match_by = "pos", exclude = NULL, verified = TRUE),

  # Nagel 2018 (CTG format; A1 = effect allele, EAF refers to A1)
  depressed_affect = list(
    file = "sumstats_depressed_affect_ctg_format", label = "Depressed affect",
    snp = "RSID", chr = "CHR", pos = "POS", ea = "A1", oa = "A2",
    effect = "BETA", effect_type = "beta", se = "SE", p = "P",
    n = "N", ncase = NULL, ncontrol = NULL, n_const = 357957, eaf = "^EAF$",
    match_by = "pos", exclude = NULL, verified = TRUE),

  depression = list(
    file = "sumstats_depression_ctg_format", label = "Depression",
    snp = "RSID", chr = "CHR", pos = "POS", ea = "A1", oa = "A2",
    effect = "Z", effect_type = "z", se = NULL, p = "P",
    n = "N", ncase = NULL, ncontrol = NULL, n_const = 381455, eaf = "^EAF_UKB$",
    match_by = "pos", exclude = NULL, verified = TRUE),

  neuroticism = list(
    file = "sumstats_neuroticism_ctg_format", label = "Neuroticism",
    snp = "RSID", chr = "CHR", pos = "POS", ea = "A1", oa = "A2",
    effect = "Z", effect_type = "z", se = NULL, p = "P",
    n = "N", ncase = NULL, ncontrol = NULL, n_const = 390278, eaf = "^EAF_UKB$",
    match_by = "pos", exclude = NULL, verified = TRUE),

  worry = list(
    file = "sumstats_worry_ctg_format", label = "Worry",
    snp = "RSID", chr = "CHR", pos = "POS", ea = "A1", oa = "A2",
    effect = "BETA", effect_type = "beta", se = "SE", p = "P",
    n = "N", ncase = NULL, ncontrol = NULL, n_const = 348219, eaf = "^EAF$",
    match_by = "pos", exclude = NULL, verified = TRUE)
)

# ---- helpers shared by the step-3 scripts ----------------------------------
env <- function(name, default) { v <- Sys.getenv(name, ""); if (nzchar(v)) v else default }
WORK       <- env("WORK", "/path/to/UKB_PSC")
RAW_DIR    <- env("RAW_DIR",    file.path(WORK, "trait_sumstats/raw"))
HARM_DIR   <- env("HARM_DIR",   file.path(WORK, "trait_sumstats/harmonized"))
SCREEN_DIR <- env("SCREEN_DIR", file.path(WORK, "ldr-sumstats/screened"))
OUT_ROOT   <- env("OUT_ROOT",   file.path(WORK, "mr_v2"))
EXP_PREFIX <- env("EXP_PREFIX", "UKB_dti_FA_phase16_")
N_FIBERS   <- as.integer(env("N_FIBERS", "430"))

# the single raw file whose name starts with cfg$file (archives are ignored)
find_raw_file <- function(cfg) {
  f <- list.files(RAW_DIR, full.names = TRUE)
  f <- f[startsWith(basename(f), cfg$file) & !grepl("\\.(zip|tar|tgz|md5|pdf|readme)$", f, ignore.case = TRUE)]
  if (length(f) > 1 && any(grepl("\\.gz$", f))) {      # same file kept both zipped and unzipped: use the unzipped one
    g <- f[!grepl("\\.gz$", f)]; if (length(g) == 1) f <- g
  }
  if (length(f) != 1) stop(sprintf("Expected exactly one file starting with '%s' in %s, found %d: %s",
                                   cfg$file, RAW_DIR, length(f), paste(basename(f), collapse = ", ")))
  f
}
# shell command that streams the file without '##' comment lines
read_cmd <- function(f) sprintf("%s %s | grep -v '^##'", if (grepl("\\.gz$", f)) "zcat" else "cat", shQuote(f))
