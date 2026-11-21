# Druggability and spatial closure

## Druggability

`step1_drug_target_lookup.py` annotates the genes at the white-matter loci with the Open Targets Platform (release 26.06): tractability, approved and clinical-stage compounds, mechanisms and indications, and a four-level tier per gene. With `--coloc-file` and `--tracts-file` it joins the colocalized traits and enriched tracts of the same lead variant. The script queries the public Open Targets GraphQL API and needs internet access.

    python step1_drug_target_lookup.py --genes-file genes.csv --out drug_target_results
    # or gene tables with Ensembl IDs, one per diffusion measure
    python step1_drug_target_lookup.py --gprofiler-file FA=FA_genelist.csv --gprofiler-file MD=MD_genelist.csv \
        --gprofiler-file AD=AD_genelist.csv --out drug_target_results

## Spatial closure

Tests whether the white-matter footprint of a colocalized locus overlaps the genetic-correlation map of a disease more than expected by chance, with a rotation null that preserves footprint size and along-fiber autocorrelation.

    sbatch step2_run.slurm                                   # 10,000 rotations per modality
    python step3_merge_and_fdr.py --dir /path/to/output/closure

Inputs of `step2_spatial_closure.py` (paths at the top of the script; `--base` overrides the root): the table of colocalized locus × trait pairs, the location-level P-value matrices of each lead variant (fibers × points), the genetic-correlation P-value matrices of each trait, and the atlas tables (`Fiber_atlas_100points_FullName.csv`, `cluster_tract_summary.csv`).
