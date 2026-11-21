#!/usr/bin/env Rscript
# =============================================================================
# 3b_select_num_ldrs.R
#
# Select the optimal number of LDRs for each fiber cluster by jointly
# considering:
#   - Image variance preserved (target: 0.80–0.90)
#   - Mean reconstruction correlation between LDR-reconstructed and raw images
#     (target: 0.85–0.95)
#
# Reads fPCA and LDR log files, extracts reported metrics, and writes a CSV
# with the selected number of LDRs for downstream steps.
#
# Usage:
#   Rscript 3b_select_num_ldrs.R
#   (Edit the Run section below to set MODALITY and directory paths.)
# =============================================================================

library(dplyr)

# ── Helpers ──────────────────────────────────────────────────────────────────

# Extract numeric cluster index from a log filename
get_cluster_idx <- function(files, prefix, suffix = ".log") {
  as.numeric(gsub(suffix, "", gsub(prefix, "", files, fixed = TRUE), fixed = TRUE))
}

# Pull the 2-row metric array immediately following key_phrase in a log
get_arr <- function(logs, key_phrase, indices, i) {
  log_i   <- logs[[which(indices == i)]]
  ldr_str <- log_i[which(log_i == key_phrase) + 1:2]
  ldr_arr <- do.call(rbind, lapply(strsplit(ldr_str, "\\s+"), as.numeric))
  ldr_arr <- ldr_arr[, !is.nan(ldr_arr[1, ]) & !is.nan(ldr_arr[2, ])]
  ldr_arr
}


# ── Core selection function ───────────────────────────────────────────────────

select_n_ldrs <- function(fpca_logs, ldr_logs, fpca_idx, ldr_idx,
                          var_thresh = c(0.80, 0.90),
                          rc_thresh  = c(0.85, 0.95)) {

  KEY_FPCA <- "The number of LDRs for preserving varying proportions of image variance:"
  KEY_LDR  <- "Mean correlation between reconstructed images and raw images using varying numbers of LDRs:"

  sel_ldr  <- rep(NA_real_, max(fpca_idx))
  min_ldr  <- rep(NA_real_, max(fpca_idx))
  flags    <- rep("None",   max(fpca_idx))

  for (i in sort(fpca_idx)) {
    if (!any(fpca_idx == i) || !any(ldr_idx == i)) next

    var_arr <- get_arr(fpca_logs, KEY_FPCA, fpca_idx, i)
    rc_arr  <- tryCatch(
      get_arr(ldr_logs, KEY_LDR, ldr_idx, i),
      error = function(e) var_arr   # fall back to variance array if unavailable
    )
    if (is.null(dim(rc_arr))) rc_arr <- var_arr

    # Candidate LDRs satisfying both constraints
    cands_var <- var_arr[2, var_arr[1, ] >= var_thresh[1] & var_arr[1, ] <= var_thresh[2]]
    cands_rc  <- rc_arr[1, rc_arr[2, ]   >= rc_thresh[1]  & rc_arr[2, ]  <= rc_thresh[2]]

    if (length(cands_rc) == 0) {
      if (all(rc_arr[2, ] < rc_thresh[1])) {
        # All RCs too low — use upper end of variance range
        cand_range <- max(cands_var)
        flags[i]   <- "No LDRs >=0.85 RC; returning LDR for ~0.95 variance"
      } else if (all(rc_arr[2, ] > rc_thresh[2])) {
        # All RCs too high — use lower end of variance range
        cand_range <- min(cands_var)
        flags[i]   <- "No LDRs <=0.95 RC; returning LDR for ~0.80 variance"
      } else {
        # Find nearest bracketing LDRs
        lo <- max(rc_arr[1, rc_arr[2, ] < rc_thresh[1]])
        hi <- min(rc_arr[1, rc_arr[2, ] > rc_thresh[2]])
        cand_range <- c(
          max(min(cands_var), lo),
          min(max(cands_var), hi)
        )
        flags[i] <- "No LDRs in 0.85-0.95 RC; using nearest range"
      }
    } else {
      cand_range <- c(
        max(min(cands_var), min(cands_rc)),
        min(max(cands_var), max(cands_rc))
      )
    }

    min_ldr[i] <- cand_range[1]
    sel_ldr[i] <- max(min_ldr[i], round(mean(cand_range)))
  }

  data.frame(
    cluster = sort(fpca_idx),
    min_ldr = min_ldr[sort(fpca_idx)],
    sel_ldr = sel_ldr[sort(fpca_idx)],
    flag    = flags[sort(fpca_idx)]
  )
}


# ── Run ──────────────────────────────────────────────────────────────────────
# Edit these variables for your setup

MODALITY <- "dti_FA"
FPCA_DIR <- paste0("/path/to/output/fpca/")
LDR_DIR  <- paste0("/path/to/output/ldrs/")
OUT_FILE <- paste0("/path/to/output/selLDRs.txt")

# Log file pattern (adjust prefix if phase tag differs)
LOG_PREFIX_FPCA <- paste0("UKB_", MODALITY, "_phase16_")
LOG_PREFIX_LDR  <- paste0("UKB_", MODALITY, "_phase16_")

fpca_log_files <- list.files(FPCA_DIR, pattern = "\\.log$", full.names = FALSE)
ldr_log_files  <- list.files(LDR_DIR,  pattern = "\\.log$", full.names = FALSE)

fpca_idx <- get_cluster_idx(fpca_log_files, prefix = LOG_PREFIX_FPCA)
ldr_idx  <- get_cluster_idx(ldr_log_files,  prefix = LOG_PREFIX_LDR)

fpca_logs <- lapply(fpca_log_files, function(x) readLines(file.path(FPCA_DIR, x)))
ldr_logs  <- lapply(ldr_log_files,  function(x) readLines(file.path(LDR_DIR,  x)))

result <- select_n_ldrs(fpca_logs, ldr_logs, fpca_idx, ldr_idx)

# Write CSV (row i+1 corresponds to cluster i for awk-based lookup in SLURM scripts)
write.csv(result$min_ldr, file = OUT_FILE, row.names = FALSE)
message("Wrote: ", OUT_FILE)
print(result)
