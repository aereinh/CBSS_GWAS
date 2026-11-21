#!/usr/bin/env Rscript
# =============================================================================
# step11_aggregate_loci.R
#
# Merge the per-tract replication results (steps 8-9) into loci and summarize
# replication per locus.
#   1. Read the per-tract locus-level replication files (step 8) and the
#      per-tract spatial-overlap files (step 9).
#   2. Join them at tract x lead-variant level.
#   3. Merge lead variants of different tracts that lie within 500 kb on the
#      same chromosome into one locus.
#   4. Per locus: lead variant (smallest P), tracts, whether the locus is
#      replicated in any tract, and the overlap fraction
#      (summed overlap locations / summed discovery locations).
#
# Usage:
#   Rscript step11_aggregate_loci.R <tract_dir> <modality> <tag> [rep_study]
#   e.g. Rscript step11_aggregate_loci.R /path/to/post_screen/tract dti_FA 1e4 phase7
#
# Outputs (in <tract_dir>):
#   ALL_TRACTS_tract_locus_replication_merged_<modality>_<tag>.csv
#   ALL_TRACTS_locus_summary_500kb_merged_<modality>_<tag>.csv
# =============================================================================

suppressMessages({
  library(dplyr)
  library(readr)
  library(stringr)
  library(purrr)
  library(tidyr)
})

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 3) stop("Usage: Rscript step11_aggregate_loci.R <tract_dir> <modality> <tag> [rep_study]")
tract_dir <- args[1]
MODALITY  <- args[2]            # e.g. dti_FA, dti_MD, dti_L1
TAG       <- args[3]            # threshold tag of the replication files, e.g. 1e4
REP_STUDY <- if (length(args) >= 4) args[4] else "phase7"

# Window for "same locus" across tracts
merge_bp <- 500000

# -----------------------------
# Helpers
# -----------------------------

# Extract tract name from filename like:
# UKB_dti_FA_phase16_Tract_ACR_clump_replication_by_phase7_1e4.txt
get_tract <- function(path) {
  bn <- basename(path)
  # tract is between "Tract_" and "_clump"
  str_match(bn, "Tract_([^_]+)_clump")[,2]
}

# Greedy merge of loci within merge_bp on each chromosome.
# Returns same data with a new integer locus_id.
assign_locus_id <- function(df, merge_bp = 500000) {
  df <- df %>% arrange(CHR, BP)
  
  locus_id <- integer(nrow(df))
  current_id <- 0
  
  for (chr in unique(df$CHR)) {
    idx <- which(df$CHR == chr)
    if (length(idx) == 0) next
    
    # Within chr, sorted by BP
    idx <- idx[order(df$BP[idx])]
    
    current_id <- current_id + 1
    locus_id[idx[1]] <- current_id
    last_bp <- df$BP[idx[1]]
    
    if (length(idx) > 1) {
      for (j in idx[-1]) {
        if ((df$BP[j] - last_bp) <= merge_bp) {
          locus_id[j] <- current_id
        } else {
          current_id <- current_id + 1
          locus_id[j] <- current_id
        }
        last_bp <- df$BP[j]
      }
    }
  }
  
  df$locus_id <- locus_id
  df
}

# Safe reader with required columns check
read_tsv_safe <- function(path) {
  tryCatch(
    read_tsv(path, show_col_types = FALSE, progress = FALSE),
    error = function(e) {
      message("Failed reading: ", path, "  |  ", conditionMessage(e))
      NULL
    }
  )
}

# -----------------------------
# 1) Find files
# -----------------------------

# helper to build pattern
make_pattern <- function(base_pattern, modality, tag) {
  paste0(".*", modality, ".*", ".*", base_pattern, tag, "\\.txt$")
}

gen_files <- list.files(
  tract_dir,
  pattern = make_pattern(paste0("clump_replication_by_",REP_STUDY,"_"), MODALITY, TAG),
  full.names = TRUE
)

spat_files <- list.files(
  tract_dir,
  pattern = make_pattern("clump_replication_spatial_overlap_", MODALITY, TAG),
  full.names = TRUE
)

if (length(gen_files) == 0) stop("No genetic replication files found in: ", tract_dir)
if (length(spat_files) == 0) message("Warning: no spatial overlap files found in: ", tract_dir)


# -----------------------------
# 2) Read genetic replication results
# -----------------------------
gen_df <- map_dfr(gen_files, function(f) {
  df <- read_tsv_safe(f)
  if (is.null(df)) return(NULL)
  
  tract <- get_tract(f)
  
  # Expect columns like CHR, SNP, BP, P, replication (plus maybe F)
  # Normalize names a bit if needed:
  df <- df %>%
    rename_with(~toupper(.x)) %>%
    mutate(TRACT = tract)
  
  # Position column may be BP or POS
  if (!("BP" %in% names(df)) && ("POS" %in% names(df))) df <- df %>% rename(BP = POS)
  
  required <- c("TRACT", "CHR", "BP", "SNP", "P", "REPLICATION")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    message("Skipping (missing cols: ", paste(missing, collapse=", "), "): ", f)
    return(NULL)
  }
  
  df %>%
    transmute(
      TRACT,
      CHR = as.integer(CHR),
      BP  = as.integer(BP),
      SNP = as.character(SNP),
      P   = as.numeric(P),
      replication_gen = as.character(REPLICATION)
    )
})

if (nrow(gen_df) == 0) stop("No usable genetic replication rows after parsing.")


# -----------------------------
# 3) Read spatial overlap results
# -----------------------------
spat_df <- map_dfr(spat_files, function(f) {
  df <- read_tsv_safe(f)
  if (is.null(df)) return(NULL)
  
  tract <- get_tract(f)
  
  # expected columns: chr, snp, pos, p, disc_voxels, rep_voxels, overlap_voxels, overlap_percent
  df <- df %>%
    rename_with(~tolower(.x)) %>%
    mutate(tract = tract)
  
  # normalize
  if (!("pos" %in% names(df)) && ("bp" %in% names(df))) df <- df %>% rename(pos = bp)
  
  required <- c("tract", "chr", "pos", "snp", "p", "disc_voxels", "rep_voxels", "overlap_voxels", "overlap_percent")
  missing <- setdiff(required, names(df))
  if (length(missing) > 0) {
    message("Skipping spatial (missing cols: ", paste(missing, collapse=", "), "): ", f)
    return(NULL)
  }
  
  df %>%
    transmute(
      TRACT = tract,
      CHR = as.integer(chr),
      BP  = as.integer(pos),
      SNP = as.character(snp),
      P   = as.numeric(p),
      disc_voxels = as.integer(disc_voxels),
      rep_voxels  = as.integer(rep_voxels),
      overlap_voxels = as.integer(overlap_voxels),
      overlap_percent = as.numeric(overlap_percent)
    )
})

# -----------------------------
# 4) Merge genetic + spatial at tract-locus level
# -----------------------------
# Primary merge keys: TRACT + CHR + BP (and maybe SNP)
# BP is the lead-variant position produced by the tract-level clumping step.
tract_locus <- gen_df %>%
  left_join(
    spat_df,
    by = c("TRACT", "CHR", "BP", "SNP", "P") # excluded "SNP", "P"
  )
cat("Unmatched loci (no spatial row):", sum(is.na(tract_locus$disc_voxels)), "\n")

# If join is too strict (SNP or P differs across outputs), fallback join on TRACT/CHR/BP
# Keep the strict join if it worked; otherwise attempt a looser one.
if (any(is.na(tract_locus$overlap_voxels)) && nrow(spat_df) > 0) {
  message("Note: some spatial rows did not join on (TRACT,CHR,BP,SNP,P). Trying fallback join on (TRACT,CHR,BP)...")
  tract_locus <- gen_df %>%
    left_join(
      spat_df %>%
        dplyr::select(TRACT, CHR, BP, disc_voxels, rep_voxels, overlap_voxels, overlap_percent),
      by = c("TRACT", "CHR", "BP")
    )
}

# Replace NAs for spatial metrics with 0 (means: no spatial overlap file row)
tract_locus <- tract_locus %>%
  mutate(
    disc_voxels = coalesce(disc_voxels, 0L),
    rep_voxels = coalesce(rep_voxels, 0L),
    overlap_voxels = coalesce(overlap_voxels, 0L),
    overlap_percent = coalesce(overlap_percent, 0)
  )

# -----------------------------
# 5) Define "global locus" IDs across tracts within 500kb
# -----------------------------
# IMPORTANT: do this across ALL tracts, by chr and BP.
tract_locus2 <- assign_locus_id(tract_locus, merge_bp = merge_bp)

# Add an easy label for reporting
tract_locus2 <- tract_locus2 %>%
  group_by(locus_id) %>%
  mutate(
    locus_chr = first(CHR),
    locus_bp_min = min(BP),
    locus_bp_max = max(BP),
    locus_bp_mid = as.integer(round(mean(BP)))
  ) %>%
  ungroup()

# -----------------------------
# 6) Locus-level summary across tracts
# -----------------------------
locus_summary <- tract_locus2 %>%
  group_by(locus_id, locus_chr, locus_bp_min, locus_bp_max, locus_bp_mid) %>%
  summarise(
    .idx = which.min(replace_na(P, Inf)),
    snp = SNP[.idx],
    chr = CHR[.idx],
    bp  = BP[.idx],
    p   = min(P, na.rm = TRUE),
    n_tracts = n_distinct(TRACT),
    tracts = paste(sort(unique(TRACT)), collapse = ","),
    gen_rep_any = any(replication_gen == "Yes"),
    n_tracts_gen_rep = n_distinct(TRACT[replication_gen == "Yes"], na.rm=T),
    disc_voxels_sum = sum(disc_voxels, na.rm = TRUE),
    rep_voxels_sum  = sum(rep_voxels, na.rm = TRUE),
    overlap_voxels_sum = sum(overlap_voxels, na.rm = TRUE),
    .groups = "drop"
  ) %>%
  mutate(
    overlap_percent = ifelse(disc_voxels_sum > 0,
                             overlap_voxels_sum / disc_voxels_sum,
                             NA_real_)
  ) %>%
  dplyr::select(-c(.idx, locus_chr, locus_bp_min, locus_bp_max, locus_bp_mid)) %>%
  arrange(chr, bp)

# -----------------------------
# 7) Save outputs
# -----------------------------
if (REP_STUDY != "ABCD") {
  out1 <- file.path(tract_dir, paste0("ALL_TRACTS_tract_locus_replication_merged_",MODALITY,"_",TAG,".csv"))
  out2 <- file.path(tract_dir, paste0("ALL_TRACTS_locus_summary_500kb_merged_",MODALITY,"_",TAG,".csv"))
} else {
  out1 <- file.path(tract_dir, paste0("ALL_TRACTS_tract_locus_replication_merged_",MODALITY,"_",TAG,"_ABCD.csv"))
  out2 <- file.path(tract_dir, paste0("ALL_TRACTS_locus_summary_500kb_merged_",MODALITY,"_",TAG,"_ABCD.csv"))
}
write.csv(tract_locus2, out1)
write.csv(locus_summary, out2)

message("Wrote tract-level merged table: ", out1)
