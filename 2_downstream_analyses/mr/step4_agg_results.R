#!/usr/bin/env Rscript
# =============================================================================
# step4_agg_results.R
#
# Aggregate MR results across all fiber clusters into a single CSV.
# Reads all *_MR_results.tsv files in MR_DIR and writes one combined file.
#
# Usage:
#   Rscript step4_agg_results.R
# =============================================================================

# ── User-configurable paths ───────────────────────────────────────────────────
MODALITY   <- "dti_FA"
MR_DIR     <- "/path/to/output/mr/results"   # Directory containing per-cluster *_MR_results.tsv
OUT_FILE   <- "/path/to/output/mr/agg/mr_all_results.csv"
N_CLUSTERS <- 430
# ─────────────────────────────────────────────────────────────────────────────

dir.create(dirname(OUT_FILE), recursive = TRUE, showWarnings = FALSE)

merged_list <- list()

for (cluster_id in seq_len(N_CLUSTERS)) {
  pattern <- paste0("UKB_", MODALITY, "_allphases_", cluster_id, "_MR_results.tsv")
  files   <- list.files(MR_DIR, pattern = glob2rx(pattern), full.names = TRUE)

  if (length(files) == 0) next

  for (f in files) {
    df <- tryCatch(read.csv(f, sep = "\t"), error = function(e) NULL)
    if (is.null(df) || nrow(df) == 0) next
    df$cluster_id <- cluster_id
    merged_list[[length(merged_list) + 1]] <- df
  }
}

if (length(merged_list) == 0) {
  stop("No MR result files found in: ", MR_DIR)
}

merged <- do.call(rbind, merged_list)
write.csv(merged, OUT_FILE, row.names = FALSE)
cat("Merged", nrow(merged), "rows from", length(merged_list), "files.\n")
cat("Written to:", OUT_FILE, "\n")
