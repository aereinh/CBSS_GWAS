#!/usr/bin/env Rscript
# Step 3d: MR for one fiber cluster (all its LDR exposures) and one trait, using the
# harmonized local trait file.
#   Rscript step3d_run_mr.R <fiber_id> <trait>
#   -> OUT_ROOT/<trait>/<EXP_PREFIX><fiber_id>_<trait>_MR_results.tsv
# Twelve estimators, as described in the manuscript.
suppressMessages({ library(TwoSampleMR); library(data.table); library(dplyr) })
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(sub("^--file=", "", args[grep("^--file=", args)]))
source(file.path(here, "trait_config.R"))
a <- commandArgs(trailingOnly = TRUE)
if (length(a) != 2) stop("Usage: Rscript step3d_run_mr.R <fiber_id> <trait>")
fiber_id <- a[1]; tr <- a[2]
if (!tr %in% names(TRAIT_CONFIG)) stop("Unknown trait ", tr)

out_dir <- file.path(OUT_ROOT, tr); dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)
output_path <- file.path(out_dir, paste0(EXP_PREFIX, fiber_id, "_", tr, "_MR_results.tsv"))
done_path   <- sub("_MR_results.tsv$", "_MR.done", output_path)   # marks tasks that finished with nothing to report
bye <- function(msg) { cat(msg, "\n"); writeLines(msg, done_path); quit(save = "no", status = 0) }

# ---- exposure: clumped instruments of this fiber cluster, one exposure per LDR (INDEX)
sfile <- file.path(SCREEN_DIR, paste0(EXP_PREFIX, fiber_id, ".txt"))
cfile <- file.path(SCREEN_DIR, paste0(EXP_PREFIX, fiber_id, "_clump_r2_0.01.clumped"))
if (!file.exists(sfile) || !file.exists(cfile)) bye("No screened/clumped LDR sumstats for this fiber cluster (no genome-wide significant SNP).")
dti_sumstats <- fread(sfile)
dti_loci <- fread(cfile)
dti_sumstats <- dti_sumstats %>% filter(SNP %in% unique(dti_loci$SNP))
dti_sumstats$SE <- dti_sumstats$BETA / dti_sumstats$Z
dti_sumstats <- as.data.frame(dti_sumstats)

dti_exp_dat <- format_data(
  dti_sumstats, type = "exposure",
  snp_col = "SNP", beta_col = "BETA", se_col = "SE",
  effect_allele_col = "A1", other_allele_col = "A2",
  pval_col = "P", samplesize_col = "N", z_col = "Z",
  chr_col = "CHR", pos_col = "POS", phenotype_col = "INDEX"
)

# ---- outcome: harmonized local trait file (already restricted to instruments; APOE removed for ALZ2)
hfile <- file.path(HARM_DIR, paste0(tr, ".tsv"))
if (!file.exists(hfile)) stop("Missing ", hfile, " (run step3c for ", tr, ")")
trait_dat <- as.data.frame(fread(hfile))
trait_dat <- trait_dat[trait_dat$SNP %in% unique(dti_sumstats$SNP), ]
if (nrow(trait_dat) == 0) bye("No instrument of this fiber cluster is present in the trait sumstats.")

disease_sumstats <- format_data(
  trait_dat, type = "outcome",
  snp_col = "SNP", beta_col = "beta", se_col = "se",
  effect_allele_col = "effect_allele", other_allele_col = "other_allele",
  eaf_col = "eaf", pval_col = "pval", samplesize_col = "samplesize",
  chr_col = "CHR", pos_col = "POS", phenotype_col = "trait"
)
disease_sumstats$id.outcome <- tr
disease_sumstats$outcome <- tr

# ---- harmonise and run MR
dat <- harmonise_data(dti_exp_dat, disease_sumstats)
if (is.null(dat) || nrow(dat) == 0 || !any(dat$mr_keep)) bye("No SNP left after harmonisation.")

mr_methods <- c("mr_ivw_fe", "mr_ivw", "mr_ivw_radial",   # 12 estimators 
                "mr_egger_regression",
                "mr_simple_median", "mr_weighted_median", "mr_penalised_weighted_median",
                "mr_simple_mode", "mr_simple_mode_nome", "mr_weighted_mode", "mr_weighted_mode_nome",
                "mr_grip")
# Fit each exposure x estimator separately, so that one estimator failing on one exposure (e.g. radial IVW's
# bootstrap optimiser on a near-perfect fit: "'xmin' not less than 'xmax'") does not stop the whole cluster.
# Failed fits are left out of the results and listed in <...>_MR_failed_fits.tsv.
res_list <- list(); failed <- list()
for (ex in unique(dat$id.exposure)) {
  d1 <- dat[dat$id.exposure == ex, ]
  for (m in mr_methods) {
    r <- tryCatch(suppressWarnings(mr(d1, method_list = m)),
                  error = function(e) { failed[[length(failed) + 1]] <<- data.frame(id.exposure = ex, exposure = d1$exposure[1],
                                          method = m, error = conditionMessage(e)); NULL })
    if (!is.null(r) && nrow(r) > 0) res_list[[length(res_list) + 1]] <- r
  }
}
fail_path <- sub("_MR_results.tsv$", "_MR_failed_fits.tsv", output_path)
if (length(failed)) {
  fwrite(do.call(rbind, failed), fail_path, sep = "\t")
  cat(length(failed), "exposure x estimator fit(s) failed and were skipped; see", fail_path, "\n")
} else if (file.exists(fail_path)) invisible(file.remove(fail_path))
if (!length(res_list)) bye("mr() returned no results (too few SNPs per exposure).")
res <- do.call(rbind, res_list)

res$fiber_id <- as.integer(fiber_id)
res$trait <- tr
fwrite(res, file = output_path, sep = "\t")
if (file.exists(done_path)) invisible(file.remove(done_path))
cat("MR results saved to:", output_path, "\n")
