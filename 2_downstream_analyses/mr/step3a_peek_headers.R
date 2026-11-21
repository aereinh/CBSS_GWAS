#!/usr/bin/env Rscript
# Step 3a: print the header and first two rows of every raw trait file, next to the
# column names currently set in trait_config.R, so the config can be filled in.
#   Rscript step3a_peek_headers.R            # all traits
#   Rscript step3a_peek_headers.R CP worry   # some traits
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(sub("^--file=", "", args[grep("^--file=", args)]))
source(file.path(here, "trait_config.R"))
want <- commandArgs(trailingOnly = TRUE); if (!length(want)) want <- names(TRAIT_CONFIG)

for (tr in want) {
  cfg <- TRAIT_CONFIG[[tr]]
  cat("\n==================== ", tr, " (verified = ", cfg$verified, ")\n", sep = "")
  f <- tryCatch(find_raw_file(cfg), error = function(e) { cat("  [MISS] ", conditionMessage(e), "\n"); NULL })
  if (is.null(f)) next
  cat("  file: ", f, "\n", sep = "")
  lines <- system(paste(read_cmd(f), "| head -n 3"), intern = TRUE)
  cat(paste0("  | ", lines), sep = "\n")
  hdr <- strsplit(trimws(lines[1]), "[ \t,]+")[[1]]
  for (k in c("snp", "chr", "pos", "ea", "oa", "effect", "se", "p", "n", "ncase", "ncontrol", "eaf")) {
    v <- cfg[[k]]; if (is.null(v)) next
    ok <- if (k == "eaf") any(grepl(v, hdr)) else v %in% hdr
    cat(sprintf("  [%s] %-8s = %s\n", if (ok) " OK " else "MISS", k, v))
  }
}
