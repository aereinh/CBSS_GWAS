# Enrichment and network analyses

The enrichment and network analyses were run with public web tools. This folder holds the script that prepares their inputs and records the settings used.

`prepare_enrichment_inputs.R` writes the lead variants and mapped genes per tract, the genes at all replicated loci, the FA-associated proteins and the FA-associated metabolites.

| Analysis | Tool | Input | Settings |
|---|---|---|---|
| SNP-to-gene mapping | FUMA v1.5.7 (SNP2GENE) | Summary statistics at the representative location of each locus | cis-eQTL mapping, GTEx v8, FDR < 0.05 |
| Gene-set enrichment | g:Profiler | Genes per tract; genes at all replicated loci; protein-coding genes of the associated proteins | GO (BP, MF, CC), KEGG, Reactome; g:SCS correction, adjusted P < 0.05 |
| Protein–protein interaction network | STRING v12.0 | FA-associated plasma proteins | Confidence score > 0.4 |
| Gene–metabolite–disease network | MetaboAnalyst (Network Analysis) | Genes and FA-associated metabolites | HMDB gene–metabolite relations; no pathway pre-filtering |
