#!/usr/bin/env Rscript
# =============================================================================
# step4_coloc_analysis.R
#
# Run colocalization (coloc.abf) between imaging voxel GWAS and external trait
# sumstats for all locus × trait files produced by step 3.
#
# Requires:
#   - trait_metadata.json: named list where each key is a trait name and the
#     value is either "quant" (quantitative) or a numeric case prevalence in
#     (0, 1) for case-control traits.  Trait names must match the suffixes of
#     files in LOCUS_DIR (e.g., *_ADHD.txt implies key "ADHD").
#
# Output:
#   coloc_results.txt with four labelled sections:
#     [RESULTS]  — successful colocalization results
#     [SKIPPED]  — files that could not be processed (with reason)
#     [WARNINGS_TOP20]  — most frequent coloc warnings
#     [WARNINGS_ALL]    — full warnings log
#
# Usage:
#   Rscript step4_coloc_analysis.R
# =============================================================================

suppressPackageStartupMessages({
  library(coloc)
  library(jsonlite)
})

# ── User-configurable paths ───────────────────────────────────────────────────
LOCUS_DIR     <- "/path/to/coloc/coloc_sumstats"    # step 3 output directory
METADATA_PATH <- "/path/to/coloc/trait_metadata.json"
OUT_FILE      <- "/path/to/coloc/coloc_results.txt"
# ─────────────────────────────────────────────────────────────────────────────


# ── Helpers ──────────────────────────────────────────────────────────────────

infer_trait <- function(file_base, trait_keys_sorted) {
  for (trait in trait_keys_sorted) {
    if (endsWith(file_base, paste0("_", trait, ".txt"))) return(trait)
  }
  NA_character_
}

first_col <- function(df, candidates) {
  for (nm in candidates) if (nm %in% names(df)) return(nm)
  NA_character_
}

as_num <- function(x) suppressWarnings(as.numeric(gsub(",", "", as.character(x), fixed = TRUE)))


# ── Load trait metadata ───────────────────────────────────────────────────────

meta <- fromJSON(txt = paste(readLines(METADATA_PATH, warn = FALSE), collapse = "\n"),
                 simplifyVector = TRUE)
trait_keys        <- names(meta)
trait_keys_sorted <- trait_keys[order(nchar(trait_keys), decreasing = TRUE)]  # longest first for suffix matching

# ── Process locus files ───────────────────────────────────────────────────────

locus_files <- list.files(LOCUS_DIR, pattern = "\\.txt$", full.names = TRUE)
results     <- list()
skipped     <- list()
warn_log    <- list()

for (file in locus_files) {
  file_base <- basename(file)
  trait     <- infer_trait(file_base, trait_keys_sorted)

  skip <- function(reason) {
    skipped[[length(skipped) + 1]] <<- data.frame(file = file_base, trait = trait,
                                                   reason = reason, stringsAsFactors = FALSE)
  }

  if (is.na(trait)) { skip("trait not found via suffix match"); next }

  meta_val <- meta[[trait]]
  if (identical(meta_val, "NA") || is.na(meta_val)) { skip("trait prevalence is NA"); next }

  if (identical(tolower(as.character(meta_val)), "quant")) {
    trait_type <- "quant"; trait_s <- NA_real_
  } else {
    val <- suppressWarnings(as.numeric(meta_val))
    if (!is.na(val) && val > 0 && val < 1) {
      trait_type <- "cc"; trait_s <- val
    } else {
      skip(paste0("unsupported metadata value: ", meta_val)); next
    }
  }

  dat <- tryCatch(
    read.table(file, header = TRUE, sep = "\t", stringsAsFactors = FALSE, check.names = FALSE),
    error = function(e) e
  )
  if (inherits(dat, "error")) { skip(paste0("read failed: ", dat$message)); next }

  req <- c("SNP", "BETA_trait", "SE_trait", "MAF", "BETA_img", "SE_img", "N_y")
  # Support both naming conventions (_trait/_img or _x/_y)
  for (sfx in list(c("_trait", "_img"), c("_x", "_y"))) {
    cands <- c("SNP", paste0("BETA", sfx[1]), paste0("SE", sfx[1]),
               "MAF", paste0("BETA", sfx[2]), paste0("SE", sfx[2]), "N_y")
    if (all(cands %in% names(dat))) { req <- cands; break }
  }
  missing <- setdiff(req, names(dat))
  if (length(missing) > 0) { skip(paste0("missing columns: ", paste(missing, collapse = ","))); next }

  beta_t   <- as_num(dat[[req[2]]])
  se_t     <- abs(as_num(dat[[req[3]]]))
  maf      <- as_num(dat[["MAF"]]); maf[maf > 0.5] <- 1 - maf[maf > 0.5]
  beta_img <- as_num(dat[[req[5]]])
  se_img   <- abs(as_num(dat[[req[6]]]))
  n_img    <- as_num(dat[["N_y"]])
  n_t_col  <- first_col(dat, c("N_x", "N"))
  n_t      <- if (!is.na(n_t_col)) as_num(dat[[n_t_col]]) else rep(NA_real_, nrow(dat))

  keep <- complete.cases(dat[["SNP"]], beta_t, se_t^2, beta_img, se_img^2, maf, n_img)
  if (trait_type == "quant") keep <- keep & !is.na(n_t)
  if (sum(keep) < 100) { skip(paste0("too few SNPs: ", sum(keep))); next }

  d1 <- list(snp = dat[["SNP"]][keep], beta = beta_img[keep], varbeta = se_img[keep]^2,
              MAF = maf[keep], N = n_img[keep], type = "quant")
  d2 <- list(snp = dat[["SNP"]][keep], beta = beta_t[keep],   varbeta = se_t[keep]^2,
              MAF = maf[keep], type = trait_type)
  if (trait_type == "cc") {
    d2$s <- trait_s
    if (!is.na(n_t_col)) d2$N <- n_t[keep]
  } else {
    d2$N <- n_t[keep]
  }

  res_coloc <- tryCatch(
    withCallingHandlers(
      coloc.abf(dataset1 = d1, dataset2 = d2),
      warning = function(w) {
        warn_log[[length(warn_log) + 1]] <<- data.frame(
          file = file_base, trait = trait, warning = conditionMessage(w), stringsAsFactors = FALSE)
        invokeRestart("muffleWarning")
      }
    ),
    error = function(e) e
  )
  if (inherits(res_coloc, "error")) { skip(paste0("coloc.abf error: ", res_coloc$message)); next }

  s   <- res_coloc$summary
  r   <- res_coloc$results
  top <- r[which.max(r$SNP.PP.H4), ]

  results[[length(results) + 1]] <- data.frame(
    file = file_base, trait = trait, trait_type = trait_type,
    s = ifelse(trait_type == "cc", trait_s, NA_real_),
    nsnps      = as.numeric(s[["nsnps"]]),
    PP.H0.abf  = as.numeric(s[["PP.H0.abf"]]),
    PP.H1.abf  = as.numeric(s[["PP.H1.abf"]]),
    PP.H2.abf  = as.numeric(s[["PP.H2.abf"]]),
    PP.H3.abf  = as.numeric(s[["PP.H3.abf"]]),
    PP.H4.abf  = as.numeric(s[["PP.H4.abf"]]),
    top_snp_h4    = top$snp,
    top_snp_h4_pp = as.numeric(top$SNP.PP.H4),
    stringsAsFactors = FALSE
  )
}

# ── Write output ──────────────────────────────────────────────────────────────

results_df <- if (length(results) > 0) do.call(rbind, results) else data.frame()
skipped_df <- if (length(skipped) > 0) do.call(rbind, skipped) else data.frame()
warn_df    <- if (length(warn_log) > 0) do.call(rbind, warn_log) else data.frame()

dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)
con <- file(OUT_FILE, open = "wt")
on.exit(close(con), add = TRUE)

writeLines(c(
  "# coloc analysis results",
  paste0("# locus_dir: ", LOCUS_DIR),
  paste0("# n_input_files: ", length(locus_files)),
  paste0("# n_success: ", nrow(results_df)),
  paste0("# n_skipped: ", nrow(skipped_df)),
  ""
), con)

writeLines("[RESULTS]", con)
if (nrow(results_df) > 0) write.table(results_df, con, sep = "\t", row.names = FALSE, quote = FALSE) else writeLines("None", con)
writeLines("", con)

writeLines("[SKIPPED]", con)
if (nrow(skipped_df) > 0) write.table(skipped_df, con, sep = "\t", row.names = FALSE, quote = FALSE) else writeLines("None", con)
writeLines("", con)

writeLines("[WARNINGS_TOP20]", con)
if (nrow(warn_df) > 0) {
  top20 <- aggregate(list(count = warn_df$warning),
                     by = list(warning = warn_df$warning), FUN = length)
  top20 <- head(top20[order(-top20$count), ], 20)
  write.table(top20, con, sep = "\t", row.names = FALSE, quote = FALSE)
} else writeLines("None", con)
writeLines("", con)

writeLines("[WARNINGS_ALL]", con)
if (nrow(warn_df) > 0) write.table(warn_df, con, sep = "\t", row.names = FALSE, quote = FALSE) else writeLines("None", con)

message("Coloc results written to: ", OUT_FILE)
