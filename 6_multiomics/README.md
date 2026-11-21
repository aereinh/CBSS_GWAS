# Multi-omics association

Associations of the white-matter LDRs with plasma proteins, NMR metabolites and clinical biomarkers in the UK Biobank, reconstructed at location level, with replication in Phase 7, and the spatial concordance between features.

The association steps (0 and 1) call the tabular association module (`main.py`; "Voxelwise analysis for tabular data", v0.0.1) developed by Zhiwen (Owen) Jiang and Hongtu Zhu as a companion to HEIG. The module is not yet publicly available and is not included in this repository; the scripts show how it was called (`/path/to/tabular`). Steps 2–4 need only Python.

| Step | Script | What it does |
|---|---|---|
| 0 | `0_make_mt_<panel>.sh` | Hail matrix tables of the omics features for the discovery and validation samples |
| 1 | `1_assoc_<panel>.slurm` | Per fiber cluster: LDR association, summary statistics and location-level reconstruction (P < 0.05 kept), in the discovery (Phases 1–6) and validation (Phase 7) samples |
| 2 | `2_effective_number_features.py` | Effective number of independent features per panel (Li & Ji) |
| 3 | `3_screen_top_features.py` | Features with significant location-level associations: Bonferroni threshold from the effective numbers of locations and features, and same-sign nominal replication in Phase 7; one table per diffusion measure |
| 4a | `4a_aggregate_selected_features.py` | Location-level statistics of the selected features across clusters |
| 4b | `4b_cross_omics_concordance.py` | Pairwise spatial concordance between features (and lead variants): correlation of location-level statistics and Dice overlap of the significant locations (FDR q < 0.05) |

`<panel>` is `proteomics`, `metabolomics` or `biomarkers`. The scripts as provided are set for FA; change the modality in the paths for MD (`dti_MD`) and AD (`dti_L1`).
