# A spatially resolved genomic and molecular atlas of human white matter: analysis code

Analysis code for the white-matter atlas built on the connectome-based spatial statistics (CBSS) representation of diffusion MRI in the UK Biobank, with validation in an independent UK Biobank release and in the ABCD Study. The repository covers image preparation from the CBSS outputs, functional PCA and low-dimensional representations (LDRs), location-level GWAS and spatial screening, replication, heritability and genetic correlation, colocalization, Mendelian randomization, multi-omics association, enrichment, druggability and spatial closure.

The imaging-genetics steps use the [HEIG](https://github.com/Zhiwen-Owen-Jiang/heig) framework and are written for SLURM-based HPC clusters.

---

## Overview

White-matter fibers of the CBSS atlas (Li et al. 2026; 6,090 fibers, 100 points each) are grouped into 430 fiber clusters. For each cluster and diffusion measure (FA, MD, AD), HEIG computes LDRs of the along-fiber profiles and runs a genome-wide scan of each LDR. Location-level association statistics are then reconstructed from the LDR summary statistics, screened with a spatial cluster-size null distribution, and replicated in held-out samples.

Samples and naming used throughout the code:

| Sample | Description | File tag |
|---|---|---|
| Discovery | UK Biobank imaging Phases 1–6 | `phase16` |
| Validation | UK Biobank imaging Phase 7 | `phase7` |
| ABCD | Adolescent Brain Cognitive Development Study (FA) | `ABCD` |

Diffusion measures are named `dti_FA`, `dti_MD` and `dti_L1` (axial diffusivity, AD) in file names.

---

## Data availability

Individual-level UK Biobank and ABCD data cannot be redistributed. They are available to approved researchers through the UK Biobank access procedure and the NBDC Data Sharing Platform. The steps that need individual-level data are folders 0, 1, 3, 4 (reconstruction steps), 5 and 6. The downstream analyses (heritability, genetic correlation, Mendelian randomization, colocalization, enrichment, druggability and spatial closure) can be reproduced from the released summary statistics.

The CBSS atlas and the projection of individual images onto it (tractography, fiber clustering, alignment) come from the connectome-based spatial statistics framework of Li et al. (2026; see Citation). The CBSS code is not part of this repository and has not been released yet. This repository starts from the CBSS outputs: the along-fiber measurements of each participant and the atlas tables (`Fiber_atlas_100points_FullName.csv`, the atlas tractogram).

---

## Software

| Tool | Version used | Purpose |
|------|---------|---------|
| [HEIG](https://github.com/Zhiwen-Owen-Jiang/heig) | 1.3.0 | fPCA, LDRs, GWAS, location-level reconstruction, heritability, genetic correlation |
| Tabular association module (Zhiwen Jiang and Hongtu Zhu) | 0.0.1 | LDR and location-level association with tabular (omics) features; not yet publicly released |
| [PLINK](https://www.cog-genomics.org/plink/) | 1.90b6 (2.00 for merging ABCD genotypes) | LD clumping, genotype handling |
| Python | 3.12.4 | Preprocessing and post-processing (numpy, pandas, scipy, dipy, requests, matplotlib) |
| R | 4.5.0 | LDR selection, MR, colocalization, locus aggregation (data.table, dplyr, tidyr, readr, purrr, stringr, jsonlite) |
| TwoSampleMR | 0.7.0 | Mendelian randomization |
| coloc | 5.2.3 | Colocalization (`coloc.abf`) |
| FUMA | 1.5.7 (web) | SNP-to-gene mapping |
| g:Profiler, STRING (v12.0), MetaboAnalyst | web | Enrichment and network analyses |
| Open Targets Platform | 26.06 | Druggability annotation |

HEIG additionally requires Java (via Conda) and Apache Spark for the GWAS and cluster-analysis steps. The tabular association module (`main.py` in the multi-omics scripts) was developed by Zhiwen (Owen) Jiang and Hongtu Zhu as a companion to HEIG. It is not yet publicly available and is not included here; the multi-omics scripts show how it was called.

---

## Repository structure

```
.
├── 0_preprocessing/                 # From CBSS outputs to HEIG inputs
│   ├── format_images.py / .slurm    # .mat along-fiber measurements -> one text file per cluster and phase (Phases 1–7)
│   ├── qc.py                        # Scan raw .mat files for outlier values
│   ├── remove_outliers.py           # Replace outliers in the text files with NA
│   ├── get_coords_clusters.py       # 3D coordinates per fiber cluster
│   └── get_coords_tracts.py         # Coordinates aggregated to tract level
│
├── 1_gwas_discovery/                # Discovery GWAS (Phases 1–6)
│   ├── 1_read_images.slurm          # Read images per phase; merge Phases 1–6; read Phase 7
│   ├── 2_fpca.slurm                 # Functional PCA
│   ├── 3a_make_ldrs.slurm           # Low-dimensional representations
│   ├── 3b_select_num_ldrs.R         # Number of LDRs per cluster
│   ├── 4_relatedness.slurm          # LOCO predictions
│   ├── 5a_ldr_gwas.slurm            # GWAS of each LDR
│   ├── 5b_preprocess_sumstats.slurm # LDR summary statistics in HEIG format
│   ├── 5c_voxelgwas_screening.slurm # Location-level GWAS (all locations)
│   └── 5d_cluster_analysis.slurm    # Bootstrap cluster-size null distribution
│
├── 2_downstream_analyses/
│   ├── heritability.slurm           # SNP heritability per location
│   ├── genetic_correlation.slurm    # Genetic correlation with 14 traits
│   ├── herigc_postprocess/          # Per-cluster output -> fibers x points matrices; tract and network-pair aggregation of r_g
│   ├── mr/                          # Mendelian randomization (see mr/README.md)
│   └── coloc/                       # Colocalization with 19 traits
│
├── 3_validation_gwas/               # Validation GWAS (Phase 7), discovery fPCA bases
│
├── 4_post_screening_replication/    # Spatial screening, locus definition, Phase 7 replication
│   ├── step1 … step10               # see "Running the pipeline"
│   └── step11_aggregate_loci.R      # Merge tract-level leads into loci; per-locus replication summary
│
├── 5_abcd_replication/              # ABCD GWAS (site-adjusted) and replication (see its README.md)
│   ├── genotype/                    # Genotype merging and Hail matrix table
│   ├── pipeline/                    # Images -> LDRs -> GWAS -> summary statistics
│   └── replication/                 # Location-level reconstruction, locus replication, spatial overlap
│
├── 6_multiomics/                    # Proteomics, metabolomics, clinical biomarkers (see its README.md)
├── 7_enrichment/                    # Inputs for g:Profiler, STRING, MetaboAnalyst (see its README.md)
├── 8_druggability_spatial_closure/  # Open Targets annotation; footprint vs r_g-map closure test
│
└── utilities/
    ├── agg_vgwas_singlesnp.py       # Whole-brain Z-score matrix for a single SNP
    └── extract_snplist_from_gwas.py # Extract a SNP list from GWAS results
```

---

## Configuration

Scripts in folders 0–4 have a `User-configurable paths` block near the top. Edit before running:

- `HEIG_DIR` — HEIG installation directory
- `MODALITY` — `dti_FA`, `dti_MD` or `dti_L1`
- all `*_DIR`, `*_OUT` and `*_FILE` variables

Python scripts take paths as command-line arguments. The MR, ABCD and spatial-closure folders keep all paths in one configuration file each (`mr_config.sh`, `pipeline/config.sh`, `replication/rep_config.sh`) or at the top of the script. All paths in this repository are placeholders (`/path/to/...`).

---

## Running the pipeline

### 0. Preprocessing

```bash
# Along-fiber measurements (.mat) -> fibersClust_<cluster>.txt, for Phases 1–7
sbatch 0_preprocessing/format_images.slurm

# QC scan of the raw .mat files for outlier values
python 0_preprocessing/qc.py \
    --phase Phase1 \
    --modality dti_FA \
    --work-root /path/to/OtherModality/ \
    --cluster-csv /path/to/Fiber_atlas_100points_FullName.csv \
    --out-prefix /path/to/output/Phase1_dti_FA_badscan

# Coordinates per cluster (once per atlas) and per tract
python 0_preprocessing/get_coords_clusters.py \
    --tractogram  /path/to/atlas_100clust_100points.trk \
    --reference   /path/to/ENIGMA_DTI_FA.nii.gz \
    --cluster-csv /path/to/Fiber_atlas_100points_FullName.csv \
    --out-dir     /path/to/coords/
python 0_preprocessing/get_coords_tracts.py \
    --tractogram    /path/to/atlas_100clust_100points.trk \
    --reference     /path/to/ENIGMA_DTI_FA.nii.gz \
    --cluster-csv   /path/to/Fiber_atlas_100points_FullName.csv \
    --cluster2tract /path/to/cluster_tract_summary.csv \
    --out-dir       /path/to/coords/tracts/
```

### 1. Discovery GWAS (Phases 1–6)

Each step is a SLURM array over the 430 fiber clusters.

```bash
sbatch 1_gwas_discovery/1_read_images.slurm          # also reads the Phase 7 images
sbatch 1_gwas_discovery/2_fpca.slurm
sbatch 1_gwas_discovery/3a_make_ldrs.slurm
Rscript 1_gwas_discovery/3b_select_num_ldrs.R        # set MODALITY at the top of the script
sbatch 1_gwas_discovery/4_relatedness.slurm
sbatch 1_gwas_discovery/5a_ldr_gwas.slurm
sbatch 1_gwas_discovery/5b_preprocess_sumstats.slurm
sbatch 1_gwas_discovery/5c_voxelgwas_screening.slurm
sbatch 1_gwas_discovery/5d_cluster_analysis.slurm
```

### 2. Downstream analyses

Heritability and genetic correlation need the summary statistics of step 5b.

```bash
sbatch 2_downstream_analyses/heritability.slurm
sbatch 2_downstream_analyses/genetic_correlation.slurm

# per-cluster output -> fibers x points matrices
python 2_downstream_analyses/herigc_postprocess/reshape_heritability.py --help          # --column H2, SE or P
python 2_downstream_analyses/herigc_postprocess/reshape_genetic_correlation.py --help   # writes GC, GC_SE and GC_P matrices

# tract-level and network-pair-level aggregation of the genetic correlations
sbatch 2_downstream_analyses/herigc_postprocess/aggregate_gc_tracts.slurm
```

`aggregate_gc_tracts.py` combines the location-level genetic correlations of each tract or network pair with a signed GLS test (correlation among locations of a fiber cluster reconstructed from the fPCA bases and the LDR covariance; clusters treated as independent) and a Cauchy combination test, applies Benjamini–Hochberg FDR across the grid, and classifies each cell as consistent direction, mixed direction or not significant. Its paths are set in `PipelineConfig` at the top of the script, where FA is named `FA` rather than `dti_FA`.

Genetic correlation uses 14 traits: ADHD, bipolar disorder, major depressive disorder, cognitive performance, educational attainment, insomnia, autism spectrum disorder, schizophrenia, Alzheimer's disease, intelligence, depressed affect, depression, neuroticism and worry. Set `EXT_SUMSTATS_DIR` and the file names in the script.

**Mendelian randomization** (`2_downstream_analyses/mr/`, details in its README). Exposures are the LDRs of each fiber cluster; outcomes are the same 14 GWAS as for genetic correlation, read from local summary-statistics files.

```bash
sbatch 2_downstream_analyses/mr/steps1-2_instruments.slurm   # instruments per cluster
cd 2_downstream_analyses/mr
bash submit_mr.sh prep        # instrument list, outcome harmonization
bash submit_mr.sh run         # MR per cluster x outcome (12 estimators)
bash submit_mr.sh agg         # merge results
python step5_tract_level_inference.py --mr-results ... --cluster-tract ... --out-dir ...
```

**Colocalization** (`2_downstream_analyses/coloc/`). Tests whether FA association signals at the spatially replicated loci colocalize with 19 external traits using `coloc.abf`. It needs trait summary statistics in a common format (`*_coloc.txt`; columns SNP, EA, NEA, beta, se, MAF, N) and a trait metadata file giving `"quant"` or the case proportion for each trait (`trait_metadata_example.json`).

```bash
sbatch 2_downstream_analyses/coloc/step1_run.slurm        # representative location per locus
sbatch 2_downstream_analyses/coloc/step2_vgwas.slurm      # single-location GWAS per locus
ls /path/to/coloc/voxelgwas/UKB_dti_FA_LOC*_CLST*_VOX*.txt > /path/to/coloc/tmp/voxelgwas_files.txt
sbatch 2_downstream_analyses/coloc/step3_run.slurm        # merge imaging and trait statistics
Rscript 2_downstream_analyses/coloc/step4_coloc_analysis.R
python 2_downstream_analyses/coloc/step5_process_results.py \
    --loci    /path/to/loci/ALL_TRACTS_locus_summary.csv \
    --coloc   /path/to/coloc/coloc_results.txt \
    --output  /path/to/coloc/coloc_results_replicated.tsv \
    --min-pph4 0.5
```

### 3. Validation GWAS (Phase 7)

Mirrors the discovery GWAS with the Phase 7 images; the fPCA bases of the discovery sample are reused.

```bash
sbatch 3_validation_gwas/1_make_ldrs_val.slurm
# run 3b_select_num_ldrs.R again with PHASE="val"
sbatch 3_validation_gwas/2_relatedness_val.slurm
sbatch 3_validation_gwas/3_ldr_gwas_val.slurm
sbatch 3_validation_gwas/4_preprocess_sumstats_val.slurm
sbatch 3_validation_gwas/5_cluster_analysis_val.slurm
```

### 4. Spatial screening, loci and replication

```bash
sbatch 4_post_screening_replication/step1_run.slurm                    # genome-wide significant SNPs per cluster
sbatch 4_post_screening_replication/step2_reconstruct_voxel_assn.slurm # location-level associations (P < 1e-5)
sbatch 4_post_screening_replication/step3_run.slurm                    # cluster-size screening
sbatch 4_post_screening_replication/step4a_ld_clump_per_cluster.slurm  # clumping per cluster (r² 0.6, then 0.1)
bash   4_post_screening_replication/step4b_ld_clump_per_tract.sh       # clumping per tract
python 4_post_screening_replication/step5_extract_candidate_snps.py --modality dti_FA --roi SCC   # per tract
sbatch 4_post_screening_replication/step6_run.slurm                    # targeted Phase 7 reconstruction
sbatch 4_post_screening_replication/step7_run.slurm                    # cluster-size screening of Phase 7 hits
sbatch 4_post_screening_replication/steps8_to_9_run.slurm              # locus replication and spatial overlap per tract
python 4_post_screening_replication/step10_select_top_voxels_per_locus.py --modality dti_FA --tract SCC --snp rs1234567
Rscript 4_post_screening_replication/step11_aggregate_loci.R /path/to/post_screen/tract dti_FA 1e4 phase7
```

### 5. ABCD replication

See `5_abcd_replication/README.md`.

### 6–8. Multi-omics, enrichment, druggability and spatial closure

See the README in each folder.

---

## Key parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| Location-level significance (FA) | `9.12e-12` | 5×10⁻⁸ ÷ effective number of independent locations |
| Cluster-size screening threshold | `1e-5` | Per-location threshold for cluster-size counting |
| Phase 7 reconstruction threshold | `1e-4` | Threshold for replication in the validation sample |
| ABCD reconstruction threshold | `1e-2` | Nominal threshold in ABCD (no cluster-size screen) |
| LD clumping, stage 1 | r²=0.6, 1000 kb | Broad block detection |
| LD clumping, stages 2–3 | r²=0.1, 500 kb | Independent loci per cluster and per tract |
| Locus merging | 500 kb | Tract-level leads within 500 kb form one locus |
| Replication window | ±250 kb | Window around the lead variant |
| Spatial replication | > 50% | Share of the discovery footprint recovered in the replication sample |
| `N_LDRS` | from `3b_select_num_ldrs.R` | Targets 80–90% of variance and 0.85–0.95 reconstruction correlation |
| MR instruments | r² < 0.01, 10 Mb | Clumping of genome-wide significant LDR variants |

---

## Required input files

| File | Description |
|------|-------------|
| `<k>_Aligned.mat`, `UKB_<phase>.txt` | CBSS along-fiber measurements and subject list per phase |
| `Fiber_atlas_100points_FullName.csv` | Cluster label of each atlas fiber |
| `cluster_tract_summary.csv` | Cluster-to-tract mapping (`cluster`, `tract_assign`) |
| `coords_fibersClust_<N>.txt` | 3D coordinates per cluster (from `get_coords_clusters.py`) |
| `eff_num_voxels.csv` | Effective number of independent locations per cluster |
| `selLDRs.txt` / `selLDRs_val.txt` | Selected number of LDRs (from `3b_select_num_ldrs.R`) |
| Genotype matrix (`.mt`) | Hail matrix tables: imputed (GWAS) and reference-panel (relatedness) |
| External GWAS summary statistics | Genetic correlation, MR and colocalization |
| LD matrices | HapMap3 SNPs, per chromosome, regularized (heritability, genetic correlation) |

---

## White-matter tracts

| Abbreviation | Tract |
|---|---|
| ACR | Anterior corona radiata |
| ALIC | Anterior limb of internal capsule |
| BCC | Body of corpus callosum |
| CGH | Cingulum (hippocampal part) |
| CP | Cerebral peduncle |
| EC | External capsule |
| FX | Fornix |
| GCC | Genu of corpus callosum |
| IFO | Inferior fronto-occipital fasciculus |
| PCR | Posterior corona radiata |
| PLIC | Posterior limb of internal capsule |
| PTR | Posterior thalamic radiation |
| RLIC | Retrolenticular part of internal capsule |
| SCC | Splenium of corpus callosum |
| SCR | Superior corona radiata |
| SLF | Superior longitudinal fasciculus |
| SS | Sagittal stratum |

---

## Citation

If you use this code, please cite the paper, "A spatially resolved genomic and molecular atlas of human white matter", together with the CBSS and HEIG papers:

> Li T, Wang X, Cole M, et al. "Connectome-based spatial statistics enabling large-scale population analyses of human connectome across cohorts." *bioRxiv* (2026). https://doi.org/10.64898/2026.04.09.717492

> Jiang Z, Stein J, Li T, et al. "Computation and resource efficient genome-wide association analysis for large-scale imaging studies." *Nature Communications* 17, 3313 (2026).
