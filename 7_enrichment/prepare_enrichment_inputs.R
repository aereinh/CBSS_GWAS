#!/usr/bin/env Rscript
# =============================================================================
# prepare_enrichment_inputs.R
#
# Prepare the gene, protein and metabolite lists that are submitted to the
# enrichment and network tools (g:Profiler, STRING, MetaboAnalyst; see
# README.md in this folder for versions and settings).
#
# Inputs
#   <results_dir>/sel_features_indeffnum_corrected_sig05.csv
#       multi-omics features with significant location-level associations
#       (6_multiomics/3_screen_top_features.py); columns VARIABLE, type, n_signif_assoc
#   <enrich_dir>/metabolite_names_full.csv         UKB field ID -> metabolite name (column ID)
#   <enrich_dir>/genom/tract/<TRACT>_genes.csv     FUMA / g:Profiler gene table per tract (column converted_alias)
#   <enrich_dir>/genom/snp_to_gene_replicated.csv  genes mapped to the replicated loci (column gene_ids)
#   <loci_file>   spatially replicated loci per tract (columns snp, tract)
#
# Outputs (in <out_dir>)
#   lead_snps_<TRACT>.txt       lead variants of the replicated loci, per tract
#   genes_<TRACT>.txt           mapped genes per tract
#   genes_replicated_loci.txt   genes mapped to all replicated loci
#   proteins_selected.txt       FA-associated plasma proteins
#   metabolites_selected.csv    FA-associated metabolites with names
#   genes_genomic_proteomic.txt union of GWAS-mapped genes and protein-coding genes
#
# Usage:
#   Rscript prepare_enrichment_inputs.R <results_dir> <enrich_dir> <loci_file> <out_dir>
# =============================================================================
suppressMessages(library(dplyr))

args <- commandArgs(trailingOnly = TRUE)
if (length(args) < 4) stop("Usage: Rscript prepare_enrichment_inputs.R <results_dir> <enrich_dir> <loci_file> <out_dir>")
results_dir <- args[1]; enrich_dir <- args[2]; loci_file <- args[3]; out_dir <- args[4]
dir.create(out_dir, recursive = TRUE, showWarnings = FALSE)

n_locations <- 609000   # 6,090 atlas fibers x 100 points

# ---- multi-omics features ---------------------------------------------------
feat_info <- read.csv(file.path(results_dir, "sel_features_indeffnum_corrected_sig05.csv"))
metab <- feat_info %>% filter(type == "metabolomics")
prot  <- feat_info %>% filter(type == "proteomics")

# metabolites associated with at least 1% of locations, with names
metab_names_dat <- read.csv(file.path(enrich_dir, "metabolite_names_full.csv"))
metab_sel <- metab %>% filter(n_signif_assoc >= 0.01 * n_locations)
metab_sel_dat <- metab_sel %>%
  mutate(ID = as.numeric(substr(VARIABLE, 1, 5))) %>%
  left_join(metab_names_dat, by = "ID")
write.csv(metab_sel_dat, file.path(out_dir, "metabolites_selected.csv"), row.names = FALSE)

# proteins associated with at least 10 locations
prot_sel <- prot %>% filter(n_signif_assoc >= 10) %>% select(VARIABLE)
writeLines(prot_sel$VARIABLE, file.path(out_dir, "proteins_selected.txt"))

# ---- lead variants and mapped genes per tract --------------------------------
loci_info_tract <- read.csv(loci_file)
loci_list <- split(loci_info_tract$snp, loci_info_tract$tract)

extract_genes <- function(gene_dat) {
  genelist <- gene_dat$converted_alias
  genelist <- genelist[genelist != "" & genelist != "None"]
  genelist <- unlist(strsplit(genelist, split = ","))
  unique(genelist)
}

tracts <- c('ACR','ALIC','BCC','CGH','CP','EC','FX','GCC','IFO','PCR','PLIC','PTR','RLIC','SCC','SCR','SLF','SS')
for (tr in tracts) {
  writeLines(as.character(loci_list[[tr]]), file.path(out_dir, paste0("lead_snps_", tr, ".txt")))
  gene_dat <- read.csv(file.path(enrich_dir, "genom", "tract", paste0(tr, "_genes.csv")))
  writeLines(extract_genes(gene_dat), file.path(out_dir, paste0("genes_", tr, ".txt")))
}

# ---- genes at all replicated loci, and union with the protein-coding genes ---
gene_dat <- read.csv(file.path(enrich_dir, "genom", "snp_to_gene_replicated.csv"))
genelist <- gene_dat$gene_ids
genelist <- genelist[genelist != ""]
genelist <- unique(unlist(strsplit(genelist, split = ",")))
writeLines(genelist, file.path(out_dir, "genes_replicated_loci.txt"))

# Ensembl IDs of the genes encoding the proteins with the widest FA associations
# (at least 0.1% of locations)
prot_genelist <- c("ENSG00000108688", "ENSG00000175866")
writeLines(unique(c(genelist, prot_genelist)), file.path(out_dir, "genes_genomic_proteomic.txt"))
