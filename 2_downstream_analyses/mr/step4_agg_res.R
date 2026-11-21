#!/usr/bin/env Rscript
# Step 4: merge the per-fiber MR results.
#   Rscript step4_agg_res.R [trait ...]      (default: every trait in trait_config.R)
#   -> OUT_ROOT/<trait>/<trait>_mr_all.csv   one file per trait
#   -> OUT_ROOT/agg/mr_all_results.csv       all traits
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(sub("^--file=", "", args[grep("^--file=", args)]))
source(file.path(here, "trait_config.R"))
want <- commandArgs(trailingOnly = TRUE); if (!length(want)) want <- names(TRAIT_CONFIG)

all <- list()
for (tr in want) {
  files <- list.files(file.path(OUT_ROOT, tr), pattern = "_MR_results\\.tsv$", full.names = TRUE)
  if (!length(files)) { cat(sprintf("%-18s no result files\n", tr)); next }
  d <- rbindlist(lapply(files, fread), fill = TRUE)
  d <- d[, .(id.exposure, id.outcome, outcome, exposure, method, nsnp, b, se, pval, fiber_id, trait)]
  fwrite(d, file.path(OUT_ROOT, tr, paste0(tr, "_mr_all.csv")))
  cat(sprintf("%-18s %4d fiber clusters, %8d rows\n", tr, uniqueN(d$fiber_id), nrow(d)))
  all[[tr]] <- d
}
if (length(all)) {
  dir.create(file.path(OUT_ROOT, "agg"), showWarnings = FALSE)
  out <- file.path(OUT_ROOT, "agg", "mr_all_results.csv")
  fwrite(rbindlist(all), out)
  cat("Merged data saved to:", out, "\n")
}
