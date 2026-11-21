#!/usr/bin/env Rscript
# Step 3b: list every instrument used by any fiber cluster (the LD-clumped SNPs of step 2),
# with its position and alleles from the screened LDR sumstats.
#   -> HARM_DIR/instruments.tsv   (SNP CHR POS A1 A2 n_fibers)
#   -> HARM_DIR/instruments_snplist.txt   (for plink --freq, run by the .slurm file)
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(sub("^--file=", "", args[grep("^--file=", args)]))
source(file.path(here, "trait_config.R"))
dir.create(HARM_DIR, recursive = TRUE, showWarnings = FALSE)

res <- vector("list", N_FIBERS); n_ok <- 0
for (fid in seq_len(N_FIBERS)) {
  sfile <- file.path(SCREEN_DIR, paste0(EXP_PREFIX, fid, ".txt"))
  cfile <- file.path(SCREEN_DIR, paste0(EXP_PREFIX, fid, "_clump_r2_0.01.clumped"))
  if (!file.exists(sfile) || !file.exists(cfile)) next
  clumped <- unique(fread(cfile, select = "SNP")$SNP)
  d <- unique(fread(sfile, select = c("SNP", "CHR", "POS", "A1", "A2")))[SNP %in% clumped]
  d[, fiber_id := fid]
  res[[fid]] <- d; n_ok <- n_ok + 1
  if (fid %% 50 == 0) cat("read", fid, "\n")
}
all <- rbindlist(res)
if (!nrow(all)) stop("No clumped instruments found in ", SCREEN_DIR)
ins <- all[, .(n_fibers = uniqueN(fiber_id)), by = .(SNP, CHR, POS, A1, A2)]
ins[, `:=`(CHR = as.integer(CHR), POS = as.integer(POS), A1 = toupper(A1), A2 = toupper(A2))]
dup <- ins[, .N, by = SNP][N > 1]
if (nrow(dup)) warning(nrow(dup), " SNP ids appear with more than one position/allele pair; keeping all")
fwrite(ins, file.path(HARM_DIR, "instruments.tsv"), sep = "\t")
fwrite(unique(ins[, .(SNP)]), file.path(HARM_DIR, "instruments_snplist.txt"), col.names = FALSE)
cat(sprintf("%d fiber clusters with instruments; %d unique instruments -> %s\n",
            n_ok, uniqueN(ins$SNP), file.path(HARM_DIR, "instruments.tsv")))
