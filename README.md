# White Matter VGWAS Pipeline

Scripts for running voxel-level genome-wide association studies (VGWAS) on white matter diffusion tensor imaging (DTI) data from the UK Biobank, using the [HEIG](https://github.com/Zhiwen-Owen-Jiang/heig) framework. The pipeline covers image preprocessing, functional PCA, LDR-based GWAS, spatial cluster screening, downstream analyses (heritability, genetic correlation), and replication in a held-out validation cohort.

---

## Overview

White matter fiber clusters (430 ROIs) are derived from a common atlas using probabilistic tractography. For each cluster, HEIG computes low-dimensional representations (LDRs) of the imaging phenotype and runs a genome-wide scan. Voxel-level association statistics are then recovered from the LDR summary statistics, screened using a spatial cluster-size null distribution, and replicated in an independent UKB release phase.

The pipeline is designed for SLURM-based HPC clusters.

---

## Dependencies

| Tool | Version | Purpose |
|------|---------|---------|
| [HEIG](https://github.com/Zhiwen-Owen-Jiang/heig) | 1.3+ | Core GWAS engine (fPCA, LDR, VGWAS, heritability, genetic correlation) |
| [PLINK 1.9](https://www.cog-genomics.org/plink/) | 1.90b6+ | LD clumping |
| Python | 3.8+ | Preprocessing and post-processing scripts |
| R | 4.0+ | LDR selection (`3b_select_num_ldrs.R`) |
| dipy | — | Tractogram coordinate extraction |
| numpy, pandas, scipy | — | General data processing |

HEIG additionally requires Java (via Conda) and Apache Spark for the GWAS and cluster analysis steps.

---

## Repository Structure

```
.
├── 0_preprocessing/
│   ├── qc.py                        # Scan raw .mat files for outlier values
│   ├── remove_outliers.py           # Replace outliers in fiber txt files with NA
│   ├── get_coords_clusters.py       # Extract 3D coordinates per fiber cluster
│   └── get_coords_tracts.py         # Aggregate coordinates to tract level
│
├── 1_gwas_discovery/
│   ├── 1_read_images.slurm          # Read images per phase, merge across phases
│   ├── 2_fpca.slurm                 # Functional PCA
│   ├── 3a_make_ldrs.slurm           # Compute Low-Dimensional Representations
│   ├── 3b_select_num_ldrs.R         # Select optimal number of LDRs from logs
│   ├── 4_relatedness.slurm          # LOCO predictions (relatedness correction)
│   ├── 5a_ldr_gwas.slurm            # LDR-level genome-wide association study
│   ├── 5b_preprocess_sumstats.slurm # Convert GWAS output to HEIG sumstats format
│   ├── 5c_voxelgwas_screening.slurm # Genome-wide voxel GWAS (all voxels)
│   └── 5d_cluster_analysis.slurm    # Bootstrap cluster-size null distribution
│
├── 2_downstream_analyses/
│   ├── heritability.slurm           # SNP-based heritability (per voxel)
│   └── genetic_correlation.slurm    # Genetic correlations with 14 external traits
│
├── 3_validation_gwas/
│   ├── 1_make_ldrs_val.slurm        # LDRs for validation cohort (Phase 6)
│   ├── 2_relatedness_val.slurm      # LOCO predictions for validation
│   ├── 3_ldr_gwas_val.slurm         # LDR GWAS in validation
│   ├── 4_preprocess_sumstats_val.slurm
│   └── 5_cluster_analysis_val.slurm # Null distribution for validation screening
│
├── 4_post_screening_replication/
│   ├── step1_get_unique_snps.py     # Extract genome-wide significant SNPs
│   ├── step1_run.slurm
│   ├── step2_reconstruct_voxel_assn.slurm  # Reconstruct voxel associations for sig SNPs
│   ├── step3_screen_snps.py         # Cluster-size screening of discovery hits
│   ├── step3_run.slurm
│   ├── step4a_ld_clump_per_cluster.slurm   # 2-stage LD clumping per cluster
│   ├── step4b_ld_clump_per_tract.sh        # Tract-level LD clumping (stage 3)
│   ├── step5_extract_candidate_snps.py     # Extract SNPs ±250kb of lead loci
│   ├── step6_run_vgwas_tract.sh     # Targeted validation VGWAS per cluster
│   ├── step6_run.slurm
│   ├── step7_screen_replication_vgwas.py   # Cluster-size screen of replication hits
│   ├── step7_run.slurm
│   ├── step8_replication_one_tract.py      # Genetic replication check per tract
│   ├── step9_spatial_overlap_one_tract.py  # Voxel-level spatial overlap
│   ├── steps8_to_9_run.slurm               # Run steps 8–9 for all tracts
│   └── step10_select_top_voxels_per_locus.py  # Per-voxel summary for replicated loci
│
└── 5_utilities/
    ├── agg_vgwas_singlesnp.py       # Whole-brain Z-score matrix for a single SNP
    └── extract_snplist_from_gwas.py # Extract SNP list from GWAS results file
```

---

## Configuration

All scripts use a clearly marked `User-configurable paths` block near the top. Before running, edit these variables in each script:

- `HEIG_DIR` — path to the HEIG installation directory
- `MODALITY` — DTI modality (e.g., `dti_FA`, `dti_MD`, `dti_RD`, `dti_L1`)
- All `*_DIR`, `*_OUT`, and `*_FILE` path variables

Python scripts accept `--main-dir` and other path arguments at the command line, so no editing is required for those.

---

## Running the Pipeline

### Step 0 — Preprocessing

**QC scan** (identify outlier fibers in raw .mat files):
```bash
python qc.py \
    --phase Phase12 \
    --modality dti_FA \
    --work-root /path/to/OtherModality/ \
    --cluster-csv /path/to/Fiber_atlas_100points_FullName.csv \
    --out-prefix /path/to/output/Phase12_dti_FA_badscan
```

**Extract cluster coordinates** (run once per atlas):
```bash
python get_coords_clusters.py \
    --tractogram  /path/to/atlas_100clust_100points.trk \
    --reference   /path/to/ENIGMA_DTI_FA.nii.gz \
    --cluster-csv /path/to/Fiber_atlas_100points_FullName.csv \
    --out-dir     /path/to/coords/
```

**Extract tract coordinates** (requires `cluster_tract_summary.csv`):
```bash
python get_coords_tracts.py \
    --tractogram    /path/to/atlas_100clust_100points.trk \
    --reference     /path/to/ENIGMA_DTI_FA.nii.gz \
    --cluster-csv   /path/to/Fiber_atlas_100points_FullName.csv \
    --cluster2tract /path/to/cluster_tract_summary.csv \
    --out-dir       /path/to/coords/tracts/
```

### Phase 1 — Discovery GWAS

Run steps sequentially. Each submits a SLURM array over all 430 fiber clusters.

```bash
# 1. Read images and merge UKB release phases
sbatch 1_gwas_discovery/1_read_images.slurm

# 2. Functional PCA
sbatch 1_gwas_discovery/2_fpca.slurm

# 3a. Make LDRs
sbatch 1_gwas_discovery/3a_make_ldrs.slurm

# 3b. Select optimal number of LDRs (run locally; edits MODALITY at top of script)
Rscript 1_gwas_discovery/3b_select_num_ldrs.R

# 4. Relatedness / LOCO predictions
sbatch 1_gwas_discovery/4_relatedness.slurm

# 5a. LDR-level GWAS
sbatch 1_gwas_discovery/5a_ldr_gwas.slurm

# 5b. Preprocess summary statistics
sbatch 1_gwas_discovery/5b_preprocess_sumstats.slurm

# 5c. Voxel GWAS screening (whole genome, all voxels)
sbatch 1_gwas_discovery/5c_voxelgwas_screening.slurm

# 5d. Bootstrap cluster-size null distribution
sbatch 1_gwas_discovery/5d_cluster_analysis.slurm
```

### Phase 2 — Downstream Analyses

These can be run in parallel with or after Phase 1 (requires sumstats from step 5b).

```bash
sbatch 2_downstream_analyses/heritability.slurm
sbatch 2_downstream_analyses/genetic_correlation.slurm
```

The genetic correlation script tests 14 traits: ADHD, Bipolar Disorder, MDD, Cognitive Performance, Educational Attainment, Insomnia, ASD, Schizophrenia, Alzheimer's Disease, Intelligence, Depressed Affect, Depression, Neuroticism, and Worry. Update `EXT_SUMSTATS_DIR` and file names as needed.

### Phase 3 — Validation GWAS

Run validation (Phase 6 / held-out cohort) GWAS to enable replication. Mirrors Phase 1 but uses Phase 6 images and validation LDRs. Note that the fPCA bases from Phase 1 are reused.

```bash
sbatch 3_validation_gwas/1_make_ldrs_val.slurm
# Run 3b_select_num_ldrs.R again with PHASE="val"
sbatch 3_validation_gwas/2_relatedness_val.slurm
sbatch 3_validation_gwas/3_ldr_gwas_val.slurm
sbatch 3_validation_gwas/4_preprocess_sumstats_val.slurm
sbatch 3_validation_gwas/5_cluster_analysis_val.slurm
```

### Phase 4 — Post-screening and Replication

```bash
# Step 1: Extract unique genome-wide significant SNPs per cluster
sbatch 4_post_screening_replication/step1_run.slurm

# Step 2: Reconstruct full voxel associations for significant SNPs (p<1e-5)
sbatch 4_post_screening_replication/step2_reconstruct_voxel_assn.slurm

# Step 3: Cluster-size screening using bootstrap null
sbatch 4_post_screening_replication/step3_run.slurm

# Step 4a: 2-stage LD clumping per cluster (r²=0.6 then r²=0.1)
sbatch 4_post_screening_replication/step4a_ld_clump_per_cluster.slurm

# Step 4b: Tract-level LD clumping — merges clusters within each tract
bash 4_post_screening_replication/step4b_ld_clump_per_tract.sh

# Step 5: Extract candidate SNPs ±250kb of each lead locus, per tract
#   (run for each tract, e.g.:)
python 4_post_screening_replication/step5_extract_candidate_snps.py \
    --modality dti_FA --roi SCC

# Step 6: Targeted replication VGWAS in validation cohort
sbatch 4_post_screening_replication/step6_run.slurm

# Step 7: Cluster-size screening of replication hits
sbatch 4_post_screening_replication/step7_run.slurm

# Steps 8–9: Genetic replication + spatial overlap per tract
sbatch 4_post_screening_replication/steps8_to_9_run.slurm

# Step 10: Per-voxel summary for a specific replicated locus
python 4_post_screening_replication/step10_select_top_voxels_per_locus.py \
    --modality dti_FA --tract SCC --snp rs1234567
```

---

## Key Parameters

| Parameter | Value | Description |
|-----------|-------|-------------|
| `SIG_THRESH` (discovery VGWAS) | `9.12e-12` | 5×10⁻⁸ ÷ effective number of independent voxels |
| `SIG_THRESH` (cluster screening) | `1e-5` | Per-voxel threshold for cluster-size counting |
| `SIG_THRESH` (validation VGWAS) | `1e-3` or `1e-4` | Relaxed threshold for replication |
| LD clumping: stage 1 | r²=0.6, 1000 kb | Broad block detection |
| LD clumping: stage 2 | r²=0.1, 500 kb | Independent loci per cluster |
| LD clumping: stage 3 | r²=0.1, 500 kb | Independent loci per tract |
| Replication window | ±250 kb | Window around lead SNP for replication check |
| `N_LDRS` | Selected by `3b_select_num_ldrs.R` | Targets 80–90% variance and 85–95% reconstruction correlation |

---

## Required Input Files

| File | Description |
|------|-------------|
| `fibersClust_<N>.txt` | Per-cluster fiber imaging data (one per UKB phase) |
| `coords_fibersClust_<N>.txt` | 3D coordinates for each fiber cluster (from `get_coords_clusters.py`) |
| `Fiber_atlas_100points_FullName.csv` | Per-fiber cluster assignments from the atlas |
| `cluster_tract_summary.csv` | Cluster-to-tract mapping (`cluster`, `tract_assign` columns) |
| `eff_num_voxels.csv` | Effective number of independent voxels per cluster (computed via HEIG) |
| `selLDRs.txt` / `selLDRs_val.txt` | Selected number of LDRs (output of `3b_select_num_ldrs.R`) |
| Genotype matrix (`.mt`) | Hail matrix table — imputed (GWAS) and reference-panel (relatedness) versions |
| External GWAS sumstats | Required for genetic correlation step |
| LD matrices | HapMap3 SNPs, per-chromosome, regularized (for heritability/genetic correlation) |

---

## White Matter Tracts

The 17 tracts and their abbreviations used in this pipeline:

| Abbreviation | Tract |
|---|---|
| ACR | Anterior Corona Radiata |
| ALIC | Anterior Limb of Internal Capsule |
| BCC | Body of Corpus Callosum |
| CGH | Cingulum (Hippocampus) |
| CP | Cerebral Peduncle |
| EC | External Capsule |
| FX | Fornix |
| GCC | Genu of Corpus Callosum |
| IFO | Inferior Fronto-Occipital Fasciculus |
| PCR | Posterior Corona Radiata |
| PLIC | Posterior Limb of Internal Capsule |
| PTR | Posterior Thalamic Radiation |
| RLIC | Retrolenticular Part of Internal Capsule |
| SCC | Splenium of Corpus Callosum |
| SCR | Superior Corona Radiata |
| SLF | Superior Longitudinal Fasciculus |
| SS | Sagittal Stratum |

---

## Citation

If you use these scripts, please cite the HEIG paper:

> Jiang Z, et al. "Imaging Genetics via High-dimensional Imaging Genetic Analysis." *bioRxiv* (2024). https://doi.org/10.1101/2024.01.26.577471
