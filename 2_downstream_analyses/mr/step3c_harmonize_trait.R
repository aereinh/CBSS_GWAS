#!/usr/bin/env Rscript
# Step 3c: reduce one raw trait file to the instrument SNPs and write it in one standard format.
#   Rscript step3c_harmonize_trait.R <trait>
#   -> HARM_DIR/<trait>.tsv   (SNP CHR POS effect_allele other_allele beta se pval eaf samplesize ncase ncontrol trait)
#   -> HARM_DIR/<trait>.log   (what was matched, dropped and excluded)
# SNP is the rsID used in the LDR GWAS; rows are matched to the instruments by chr:pos + alleles
# (either strand), so the trait file needs no rsID. Allele ORIENTATION is left to TwoSampleMR.
suppressMessages(library(data.table))
args <- commandArgs(trailingOnly = FALSE)
here <- dirname(sub("^--file=", "", args[grep("^--file=", args)]))
source(file.path(here, "trait_config.R"))
tr <- commandArgs(trailingOnly = TRUE)[1]
if (is.na(tr) || !tr %in% names(TRAIT_CONFIG)) stop("Usage: step3c_harmonize_trait.R <trait>; traits: ", paste(names(TRAIT_CONFIG), collapse = ", "))
cfg <- TRAIT_CONFIG[[tr]]
if (!isTRUE(cfg$verified)) stop(tr, ": verified = FALSE in trait_config.R. Run step3a_peek_headers.R, check the column names and the effect allele, then set verified = TRUE.")

logf <- file.path(HARM_DIR, paste0(tr, ".log")); if (file.exists(logf)) invisible(file.remove(logf))
say <- function(...) { m <- sprintf(...); cat(m, "\n"); cat(m, "\n", file = logf, append = TRUE) }

ins <- fread(file.path(HARM_DIR, "instruments.tsv"))
f <- find_raw_file(cfg)
say("Trait %s (%s): %s", tr, cfg$label, f)

# ---- read only the columns needed
hdr <- names(fread(cmd = paste(read_cmd(f), "| head -n 5")))
col <- function(k, required = FALSE) {
  v <- cfg[[k]]; if (is.null(v)) { if (required) stop(tr, ": config field '", k, "' is required"); return(NULL) }
  hit <- if (k == "eaf") grep(v, hdr, value = TRUE) else intersect(v, hdr)
  if (length(hit) != 1) stop(sprintf("%s: column for '%s' (%s) matched %d columns.\nHeader: %s", tr, k, v, length(hit), paste(hdr, collapse = " ")))
  hit
}
by_pos <- !identical(cfg$match_by, "rsid")
cols <- list(snp = col("snp", !by_pos), chr = col("chr", by_pos), pos = col("pos", by_pos), ea = col("ea", TRUE), oa = col("oa", TRUE),
             effect = col("effect", TRUE), se = col("se"), p = col("p"), n = col("n"),
             ncase = col("ncase"), ncontrol = col("ncontrol"), eaf = col("eaf"))
cols <- cols[!sapply(cols, is.null)]
d <- fread(cmd = read_cmd(f), select = unname(unlist(cols)), showProgress = FALSE)
setnames(d, unname(unlist(cols)), names(cols))
say("Read %d rows", nrow(d))
d[, `:=`(ea = toupper(ea), oa = toupper(oa))]
if ("chr" %in% names(d)) {
  d[, chr := sub("^chr", "", as.character(chr), ignore.case = TRUE)]
  d <- d[chr %in% as.character(1:22)][, `:=`(chr = as.integer(chr), pos = as.integer(pos))]
}

# ---- explicit region exclusion (APOE for Alzheimer's)
if (!is.null(cfg$exclude)) {
  if (!"chr" %in% names(d)) stop(tr, ": 'exclude' needs chr and pos columns")
  for (reg in cfg$exclude) {
    rg <- as.numeric(strsplit(reg, "[:-]")[[1]])
    hit <- d$chr == rg[1] & d$pos >= rg[2] & d$pos <= rg[3]
    n_ins <- nrow(ins[CHR == rg[1] & POS >= rg[2] & POS <= rg[3]])
    say("Excluded region %s (GRCh37): removed %d variants from the trait file; %d instrument(s) lie in it and are dropped for this trait", reg, sum(hit), n_ins)
    d <- d[!hit]
    stopifnot(!any(d$chr == rg[1] & d$pos >= rg[2] & d$pos <= rg[3]))
  }
}

# ---- match to the instruments
if (by_pos) {
  m <- merge(ins, d, by.x = c("CHR", "POS"), by.y = c("chr", "pos"))
} else {
  m <- merge(ins, d, by.x = "SNP", by.y = "snp")
}
comp <- function(x) chartr("ACGT", "TGCA", x)
same   <- (m$ea == m$A1 & m$oa == m$A2) | (m$ea == m$A2 & m$oa == m$A1)
flip   <- (comp(m$ea) == m$A1 & comp(m$oa) == m$A2) | (comp(m$ea) == m$A2 & comp(m$oa) == m$A1)
say("Instruments: %d; found in trait file: %d; alleles agree: %d (%d only on the other strand); allele mismatch dropped: %d",
    uniqueN(ins$SNP), uniqueN(m$SNP), uniqueN(m$SNP[same | flip]), uniqueN(m$SNP[flip & !same]), uniqueN(m$SNP[!(same | flip)]))
m <- m[same | flip]

# ---- sample size
m[, samplesize := NA_real_]
if ("n" %in% names(m)) m[, samplesize := as.numeric(n)]
if (all(c("ncase", "ncontrol") %in% names(m))) m[is.na(samplesize), samplesize := as.numeric(ncase) + as.numeric(ncontrol)]
if (!is.null(cfg$n_const)) m[is.na(samplesize), samplesize := cfg$n_const]
for (k in c("ncase", "ncontrol")) if (!k %in% names(m)) m[, (k) := NA_real_]

# ---- effect-allele frequency: the file's own column, else another trait's file (cfg$eaf_from)
if (!"eaf" %in% names(m)) {
  m[, eaf := NA_real_]
  src <- cfg$eaf_from
  if (!is.null(src)) {
    scfg <- TRAIT_CONFIG[[src]]
    if (is.null(scfg) || is.null(scfg$eaf)) stop(tr, ": eaf_from = '", src, "' is not a trait with an eaf column")
    sf <- find_raw_file(scfg)
    shdr <- names(fread(cmd = paste(read_cmd(sf), "| head -n 5")))
    seaf <- grep(scfg$eaf, shdr, value = TRUE)
    if (length(seaf) != 1) stop(src, ": eaf column '", scfg$eaf, "' matched ", length(seaf), " columns")
    s <- fread(cmd = read_cmd(sf), select = c(scfg$chr, scfg$pos, scfg$ea, scfg$oa, seaf), showProgress = FALSE)
    setnames(s, c("chr", "pos", "s_ea", "s_oa", "s_eaf"))
    s[, chr := suppressWarnings(as.integer(sub("^chr", "", as.character(chr), ignore.case = TRUE)))]
    s <- s[!is.na(chr)][, `:=`(pos = as.integer(pos), s_ea = toupper(s_ea), s_oa = toupper(s_oa), s_eaf = as.numeric(s_eaf))]
    s <- unique(s[paste(chr, pos) %chin% paste(m$CHR, m$POS)], by = c("chr", "pos", "s_ea", "s_oa"))
    m[, row_id := .I]
    j <- merge(m[, .(row_id, CHR, POS, ea, oa)], s, by.x = c("CHR", "POS"), by.y = c("chr", "pos"))
    j[, f := fifelse(s_ea == ea & s_oa == oa, s_eaf,
             fifelse(s_ea == oa & s_oa == ea, 1 - s_eaf,
             fifelse(comp(s_ea) == ea & comp(s_oa) == oa, s_eaf,
             fifelse(comp(s_ea) == oa & comp(s_oa) == ea, 1 - s_eaf, NA_real_))))]
    j[s_ea == comp(s_oa), f := NA_real_]          # A/T and C/G: strand cannot be told from alleles alone
    j <- j[!is.na(f)][, .(f = f[1]), by = row_id]
    m[j, on = "row_id", eaf := i.f]
    m[, row_id := NULL]
    say("EAF taken from %s (%s, column %s) for %d of %d instrument rows; %d without a usable frequency (absent there, allele mismatch or A/T, C/G)",
        src, basename(sf), seaf, sum(!is.na(m$eaf)), nrow(m), sum(is.na(m$eaf)))
  } else if (cfg$effect_type == "z") stop(tr, ": effect_type = 'z' needs allele frequencies: set eaf or eaf_from in trait_config.R")
}
m[, eaf := as.numeric(eaf)]

# ---- beta and standard error
m[, effect := as.numeric(effect)]
if ("p" %in% names(m)) m[, pval := as.numeric(p)] else m[, pval := NA_real_]
if (cfg$effect_type == "z") {
  m[, se := 1 / sqrt(2 * eaf * (1 - eaf) * (samplesize + effect^2))]
  m[, beta := effect * se]
  m[is.na(pval), pval := 2 * pnorm(-abs(effect))]
  say("beta/se derived from z, N and EAF (standardized scale): beta = z / sqrt(2 p (1-p) (N + z^2))")
} else {
  m[, beta := if (cfg$effect_type == "or") log(effect) else effect]
  if ("se" %in% names(m)) m[, se := as.numeric(se)] else m[, se := abs(beta / qnorm(pval / 2))]
  m[is.na(pval), pval := 2 * pnorm(-abs(beta / se))]
}
bad <- !is.finite(m$beta) | !is.finite(m$se) | m$se <= 0
if (any(bad)) say("Dropped %d rows without a usable beta/se (e.g. missing EAF or N)", sum(bad))
m <- m[!bad]

# ---- one row per instrument (keep the smallest P if a position is listed twice)
setorder(m, SNP, pval); m <- unique(m, by = "SNP")
out <- m[, .(SNP, CHR, POS, effect_allele = ea, other_allele = oa, beta, se, pval, eaf, samplesize, ncase, ncontrol, trait = tr)]
fwrite(out, file.path(HARM_DIR, paste0(tr, ".tsv")), sep = "\t")
say("Wrote %d instruments -> %s", nrow(out), file.path(HARM_DIR, paste0(tr, ".tsv")))
