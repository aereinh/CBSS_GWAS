# Mendelian randomization

Two-sample MR of white-matter FA on 14 outcomes. Exposures are the low-dimensional representations (LDRs) of each fiber cluster, each LDR a separate exposure. Outcomes are the 14 GWAS used for genetic correlation, read from local summary-statistics files.

All paths are set in `mr_config.sh`; the outcome files and their column names are described in `trait_config.R`.

## Steps

| Step | Script | What it does |
|---|---|---|
| 1 | `step1_screen_ldr_sumstats.py` | Genome-wide significant variants of the LDR GWAS of one cluster |
| 2 | `step2_ld_clump.sh` | LD clumping (r² < 0.01 within 10 Mb) to independent instruments |
| 3a | `step3a_peek_headers.R` | Prints the header of each outcome file next to the column names in `trait_config.R` |
| 3b | `step3b_build_instruments.R` | Union of the instruments of all clusters (SNP, CHR, POS, A1, A2) |
| 3c | `step3c_harmonize_trait.R <trait>` | Reads one outcome file and keeps the instruments, matched by GRCh37 chr:pos and alleles on either strand; converts odds ratios to log odds ratios, and z-scores to effects and standard errors |
| 3d | `step3d_run_mr.R <cluster> <trait>` | TwoSampleMR with 12 estimators; each exposure × estimator is fitted separately, and fits that cannot be computed are listed in `*_MR_failed_fits.tsv` |
| 4 | `step4_agg_res.R` | Per-trait and merged result tables |
| 5 | `step5_tract_level_inference.py` | Tract-level inference (rule below) |

## Run

```bash
sbatch steps1-2_instruments.slurm      # steps 1-2, one task per fiber cluster

module load r
Rscript step3a_peek_headers.R          # check the outcome files against trait_config.R
bash submit_mr.sh prep                 # 3b, then 3c for the 14 outcomes
bash check_mr.sh                       # progress
bash submit_mr.sh run                  # 3d: one array per outcome
bash submit_mr.sh agg                  # 4
python step5_tract_level_inference.py \
    --mr-results /path/to/UKB_PSC/mr_v2/agg/mr_all_results.csv \
    --cluster-tract /path/to/cluster_tract_summary.csv \
    --out-dir /path/to/UKB_PSC/mr_v2/agg
```

One outcome at a time: `bash submit_mr.sh harmonize CP`, then `bash submit_mr.sh run CP`. Selected clusters: `bash submit_mr.sh run ALZ2 --array 17,203`.

## Outcome handling

- Outcome files are assumed to be GRCh37, like the LDR GWAS. The number of instruments found in each file is written to `<trait>.log`.
- z-score-only files (Alzheimer's disease, depression, neuroticism): `beta = z / sqrt(2p(1-p)(N + z^2))`, with `p` the effect-allele frequency; these effects are on the standardized scale. The Alzheimer's file has no allele frequency, so it is taken from another file of the panel (`eaf_from` in `trait_config.R`), aligned by position and alleles.
- Alzheimer's disease: the APOE region (chr19:44,400,000–46,500,000, GRCh37) is removed from the outcome before matching (`exclude` in `trait_config.R`).
- Without exposure allele frequencies, TwoSampleMR drops palindromic (A/T, C/G) instruments.

## Tract-level rule (step 5)

1. Keep cluster-level results with at least 5 instruments.
2. For each tract × outcome × estimator, combine the cluster-level P-values with the Cauchy combination test and compute the −log10(P)-weighted fraction of fits whose sign agrees with the fixed-effect IVW estimate of the same exposure.
3. Benjamini–Hochberg FDR across the tract × outcome grid, separately per estimator.
4. A pair is reported when fixed-effect IVW has q < 0.05 and at least 3 of the 12 estimators have q < 0.05 with weighted sign concordance 1.
