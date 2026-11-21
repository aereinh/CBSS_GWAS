# ABCD replication

GWAS of the CBSS fiber-cluster representations in the ABCD Study and replication of the UK Biobank discovery loci. In the paper this analysis is reported for FA, with imaging site as a covariate.

## genotype/

- `merge_bed_files.slurm` — merge the per-chromosome PLINK files.
- `make_mt.slurm` — Hail matrix table of the ABCD genotypes for HEIG.

## pipeline/

All paths and covariates are in `config.sh`. Every step is submitted in the same way:

    bash submit.sh <MOD> <STEP> [--array SPEC] [--after JOBID]

| STEP | Script | What it does |
|---|---|---|
| format | `0_format_images.slurm` | `<k>_Aligned.mat` → `fibersClust_<IDX>.txt` |
| read | `1_read_images.slurm` | text → `images.h5`; restores the full participant IDs (`fix_h5_ids.py`) |
| ldr | `3.1_make_ldrs.slurm` | LDRs from the UK Biobank Phases 1–6 fPCA bases, adjusted for covariates |
| select | `3.2_select_n_ldrs.slurm` | number of LDRs per cluster (`selLDRs_val.txt`) |
| gwas | `4-5_rel_gwas_sumstats.slurm` | relatedness (LOCO) → GWAS → summary statistics |

Covariates: sex, age, age², age × sex, age² × sex, 40 genetic principal components and site (categorical: sex, site). `00_prep_covar_site.py` adds the site of the earliest visit to the covariate file.

    bash check_inputs.sh FA                      # inputs and finished clusters per step
    J1=$(bash submit.sh FA format)
    J2=$(bash submit.sh FA read   --after $J1)
    J3=$(bash submit.sh FA ldr    --after $J2)
    J4=$(bash submit.sh FA select --after $J3)
    J5=$(bash submit.sh FA gwas   --after $J4)

`run_all.sh <MOD> [first_step]` chains all steps. Each step skips clusters whose output exists (`OVERWRITE=1` forces a redo); `--array 17,203` resubmits selected clusters.

## replication/

Settings are in `rep_config.sh`, which reads `../pipeline/config.sh`.

1. `1_voxel_gwas.slurm` — for each fiber cluster with a tract assignment, reconstruct the ABCD location-level associations of the discovery candidate SNPs of that tract and keep P < 0.01. No cluster-size screen is applied in ABCD.
2. `2_replication_one_tract.py` — a discovery locus is replicated when an ABCD association lies within ±250 kb of its lead variant.
3. `3_spatial_overlap_one_tract.py` — share of the discovery footprint of the locus that is also in the ABCD footprint.
4. `4_aggregate_loci.py` — tract-level leads within 500 kb are merged into loci; overlap fraction = summed overlap / summed discovery locations; a locus is spatially replicated when the fraction exceeds 0.5. Also writes the tract tables, a Fisher test of tract over-representation among the reproduced loci, and the PhenoGram input for the ideogram.

    bash check_rep.sh
    bash submit_rep.sh all        # 1 → 2–3 → 4

`4_aggregate_loci.py` compares ABCD with UK Biobank Phase 7 when the Phase 7 per-locus table (`ALL_TRACTS_locus_summary_500kb_merged_FA_1e4.csv`, from `4_post_screening_replication/step11_aggregate_loci.R`) and the top colocalized trait per locus (`coloc_top_trait_per_locus.tsv`) are placed in `replication/ref/`.
