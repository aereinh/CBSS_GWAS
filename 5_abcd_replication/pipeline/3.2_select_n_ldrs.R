# Step 3.2: choose the number of LDRs per fiber cluster for ABCD.
# Same selection rule as 3.2_select_num_ldrs_from_logs.ipynb (UKB fPCA variance + ABCD
# reconstruction correlation); writes the minLDR column, as before.
#
# Usage (normally via 3.2_select_n_ldrs.slurm):
#   Rscript 3.2_select_n_ldrs.R <ukb_fpca_dir> <ukb_fpca_prefix> <abcd_ldr_dir> <abcd_ldr_prefix> <out_file>

args <- commandArgs(trailingOnly = TRUE)
if (length(args) != 5) stop("need 5 arguments: fpca_dir fpca_prefix ldr_dir ldr_prefix out_file")
fpca_dir    <- args[1]
fpca_prefix <- args[2]   # e.g. UKB_FA_phase16
ldr_dir     <- args[3]
ldr_prefix  <- args[4]   # e.g. ABCD_FA
out_file    <- args[5]

get_fib_cl_ind <- function(file, prefix_key, suffix_key = ".log") {
  as.numeric(gsub(suffix_key, "", gsub(prefix_key, "", file), fixed = TRUE))
}

get_arr <- function(logs, key_phrase, ind, i) {
  log_i <- logs[[which(ind == i)]]
  ldr_str <- log_i[which(log_i == key_phrase) + 1:2]
  ldr_arr <- do.call(rbind, lapply(strsplit(ldr_str, "\\s+"), as.numeric))
  ldr_arr <- ldr_arr[, !is.nan(ldr_arr[1, ]) & !is.nan(ldr_arr[2, ])]
  ldr_arr
}

key_phrase_fpca <- "The number of LDRs for preserving varying proportions of image variance:"
key_phrase_ldr  <- "Mean correlation between reconstructed images and raw images using varying numbers of LDRs:"
var_thresh <- c(0.8, 0.9)
rc_thresh  <- c(0.85, 0.95)

# Selection rule copied unchanged from the original notebook (and the UKB R script),
# including the all(ldr_cands_rc) comparisons, so ABCD and UKB use identical logic.
get_nLDRs <- function(fpca_logs, ldr_logs, fpca_ind, ldr_ind) {
  minLDR <- rep(NA, length(fpca_ind))
  selLDR <- rep(NA, length(fpca_ind))
  flags  <- rep("None", length(fpca_ind))
  for (i in sort(fpca_ind)) {
    if (any(fpca_ind == i) & any(ldr_ind == i)) {
      ldr_var_arr <- get_arr(fpca_logs, key_phrase_fpca, fpca_ind, i)
      ldr_rc_arr  <- get_arr(ldr_logs, key_phrase_ldr, ldr_ind, i)
      if (is.null(dim(ldr_rc_arr))) ldr_rc_arr <- ldr_var_arr

      ldr_cands_var <- ldr_var_arr[2, which(ldr_var_arr[1, ] >= var_thresh[1] &
                                              ldr_var_arr[1, ] <= var_thresh[2])]
      ldr_cands_rc  <- ldr_rc_arr[1, which(ldr_rc_arr[2, ] >= rc_thresh[1] &
                                             ldr_rc_arr[2, ] <= rc_thresh[2])]
      if (length(ldr_cands_rc) == 0) {
        if (all(ldr_rc_arr[2, ] < rc_thresh[1])) {
          ldr_cand_range <- max(ldr_cands_var)
          flags[i] <- "No LDRs with >0.85 reconstruction correlation -- returning LDR with ~0.95 variance preserved"
        } else if (all(ldr_rc_arr[2, ] > rc_thresh[2])) {
          ldr_cand_range <- min(ldr_cands_var)
          flags[i] <- "No LDRs with <0.95 reconstruction correlation -- returning LDR with ~0.80 variance preserved"
        } else {
          min_ldr_rc <- max(ldr_rc_arr[1, which(ldr_rc_arr[2, ] < rc_thresh[1])])
          max_ldr_rc <- min(ldr_rc_arr[1, which(ldr_rc_arr[2, ] > rc_thresh[2])])
          ldr_cands_rc <- c(min_ldr_rc, max_ldr_rc)
          if (all(ldr_cands_rc) < min(ldr_cands_var)) {
            ldr_cand_range <- min(ldr_cands_var)
          } else if (all(ldr_cands_rc) > max(ldr_cands_var)) {
            ldr_cand_range <- max(ldr_cands_var)
          } else {
            ldr_cand_range <- c(max(min(ldr_cands_var), min(ldr_cands_rc)),
                                min(max(ldr_cands_var), max(ldr_cands_rc)))
          }
          flags[i] <- "No LDRs found with 0.85-0.95 reconstruction correlation -- using nearest range"
        }
      } else {
        ldr_cand_range <- c(max(min(ldr_cands_var), min(ldr_cands_rc)),
                            min(max(ldr_cands_var), max(ldr_cands_rc)))
      }
      minLDR[i] <- ldr_cand_range[1]
      selLDR[i] <- max(minLDR[i], round(mean(ldr_cand_range)))
    }
  }
  data.frame(ind = sort(fpca_ind), minLDR = minLDR, selLDR = selLDR, flag = flags)
}

# only <prefix>_<number>.log, so stray logs in the folder are ignored
fpca_log_files <- list.files(fpca_dir, pattern = paste0("^", fpca_prefix, "_[0-9]+\\.log$"))
ldr_log_files  <- list.files(ldr_dir,  pattern = paste0("^", ldr_prefix,  "_[0-9]+\\.log$"))
cat("UKB fPCA logs:", length(fpca_log_files), " ABCD LDR logs:", length(ldr_log_files), "\n")

fpca_ind <- get_fib_cl_ind(fpca_log_files, paste0(fpca_prefix, "_"))
ldr_ind  <- get_fib_cl_ind(ldr_log_files,  paste0(ldr_prefix, "_"))

fpca_logs <- lapply(fpca_log_files, function(x) readLines(file.path(fpca_dir, x)))
ldr_logs  <- lapply(ldr_log_files,  function(x) readLines(file.path(ldr_dir, x)))

nldr_dat <- get_nLDRs(fpca_logs, ldr_logs, fpca_ind, ldr_ind)

# same format as before: row i+1, column 2 = number of LDRs for cluster i
write.csv(nldr_dat$minLDR, file = out_file)
write.csv(nldr_dat, file = sub("\\.txt$", "_detail.csv", out_file), row.names = FALSE)

n_na <- sum(is.na(nldr_dat$minLDR))
cat("Wrote", out_file, "\n")
cat("Clusters with a selected number:", sum(!is.na(nldr_dat$minLDR)), " missing:", n_na, "\n")
if (n_na > 0) cat("Missing clusters:", nldr_dat$ind[is.na(nldr_dat$minLDR)], "\n")
print(table(nldr_dat$flag))
