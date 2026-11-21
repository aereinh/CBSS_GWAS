#!/usr/bin/env python3
"""
3_screen_top_features.py

Runs the discovery-screen / validation-replication summary (Bonferroni
threshold from eff_num_vox x eff_num_feature, nominal same-sign replication)
across every fiber cluster for ONE modality's omics types, and writes that
modality's combined result to its own CSV.

Usage:
    python 3_screen_top_features.py --modality FA
    python 3_screen_top_features.py --modality dti_MD --omics-types proteomics,metabolomics
"""
import os
import time
import argparse
import pandas as pd

main_dir = '/path/to/UKB_PSC/multiomics/'
default_out_dir = '/path/to/UKB_PSC/multiomics/results/'
eff_num_vox = {'FA': 5479.479, 'dti_MD': 2946.59, 'dti_L1': 3895.18}
eff_num_feature = {'proteomics': 1242.9014450393718, 'metabolomics': 33.92531132352527, 'biomarkers': 16.814971013680275}
DEFAULT_N_CLUSTERS = 430
LOG_EVERY = 20

# Fiber clusters with only 1 selected LDR for a given modality -- with a
# single basis function, the per-voxel Z-score back-projection has no
# variance to draw from and is degenerate (Z -> Inf for at least some
# voxels; empirically confirmed to actually produce literal Inf/near-0-P
# for a subset of these -- 16, 36, 152, 228, 238, 240, 409 -- via direct
# inspection of the raw voxel-assoc files). The rest (1, 9, 13, 161, 167,
# 227, 293) didn't blow up to Inf this run but are built from the same
# single, badly underpowered basis function, so they're excluded too
# rather than trusted just because they stayed finite. Add entries for
# other modalities here if the same n_ldrs<=1 issue turns up there.
EXCLUDE_CLUSTERS = {
    'dti_MD': {1, 9, 13, 16, 36, 152, 161, 167, 227, 228, 238, 240, 293, 409},
}


def process_cluster(idx, modality, omics_type, disc_thresh):
    """Same screening/replication logic as before, per fiber cluster.

    Returns None (with a [skip] print) if the cluster is excluded
    (n_ldrs<=1 for this modality -- see EXCLUDE_CLUSTERS) or if the
    voxel-assoc files for this cluster are missing, rather than raising --
    some cluster indices may not exist for a given modality/omics type."""
    if idx in EXCLUDE_CLUSTERS.get(modality, ()):
        print(f"[skip] {modality}/{omics_type} cluster {idx}: excluded (n_ldrs<=1 for {modality}, degenerate Z)")
        return None

    disc_path = f'{main_dir}/{modality}/{omics_type}/disc/voxel-assoc/voxel_assoc_{idx}.txt'
    val_path = f'{main_dir}/{modality}/{omics_type}/val/voxel-assoc/voxel_assoc_{idx}.txt'
    if not (os.path.exists(disc_path) and os.path.exists(val_path)):
        print(f"[skip] {modality}/{omics_type} cluster {idx}: missing disc and/or val voxel-assoc file")
        return None
    disc = pd.read_csv(disc_path, sep='\t')
    val = pd.read_csv(val_path, sep='\t')
    disc_sig = disc[disc["P"] < disc_thresh]
    val_sig = val[val["P"] < 0.05]
    disc_summary = (
        disc_sig.groupby("VARIABLE")["INDEX"]
        .nunique()
        .rename("n_disc")
    )
    val_summary = (
        val_sig.groupby("VARIABLE")["INDEX"]
        .nunique()
        .rename("n_val")
    )
    overlap = (
        disc_sig[["VARIABLE", "INDEX", "Z"]]
        .drop_duplicates(subset=["VARIABLE", "INDEX"])
        .merge(
            val_sig[["VARIABLE", "INDEX", "Z"]]
            .drop_duplicates(subset=["VARIABLE", "INDEX"]),
            on=["VARIABLE", "INDEX"],
            how="inner",
            suffixes=("_disc", "_val")
        )
        .query("Z_disc * Z_val > 0")
        .groupby("VARIABLE")["INDEX"]
        .nunique()
        .rename("n_overlap")
    )
    summary = (
        disc_summary.to_frame()
        .join(val_summary, how="left")
        .join(overlap, how="left")
        .fillna(0)
    )
    summary[["n_disc", "n_val", "n_overlap"]] = (
        summary[["n_disc", "n_val", "n_overlap"]].astype(int)
    )
    summary["overlap_pct"] = 100 * summary["n_overlap"] / summary["n_disc"]
    summary = summary.reset_index()
    summary.insert(0, "cluster", idx)
    return summary


def process_combo(modality, omics_type, combo_i, n_combos, t_pipeline_start, n_clusters):
    """Runs process_cluster for every fiber cluster under one (modality,
    omics_type) combination, tagging the result with MODALITY/TYPE columns."""
    K = eff_num_vox[modality] * eff_num_feature[omics_type]
    disc_thresh = 0.05 / K
    print(f"\n=== [{combo_i}/{n_combos}] {modality} x {omics_type} "
          f"(disc_thresh = 0.05 / {K:.2f} = {disc_thresh:.3e}) ===", flush=True)
    t_combo_start = time.time()
    all_summaries = []
    n_skipped = 0
    for idx in range(1, n_clusters + 1):
        res = process_cluster(idx, modality, omics_type, disc_thresh)
        if res is None:
            n_skipped += 1
        else:
            all_summaries.append(res)
        if idx == 1 or idx % LOG_EVERY == 0 or idx == n_clusters:
            elapsed = time.time() - t_combo_start
            avg = elapsed / idx
            eta = avg * (n_clusters - idx)
            print(f"[progress][{modality}/{omics_type}] {idx}/{n_clusters} clusters done "
                  f"({n_skipped} skipped so far) -- {elapsed:.1f}s elapsed, "
                  f"{avg:.2f}s/cluster avg, ETA ~{eta:.1f}s.", flush=True)
    combo_df = pd.concat(all_summaries, ignore_index=True) if all_summaries else pd.DataFrame()
    if not combo_df.empty:
        combo_df.insert(0, "TYPE", omics_type)
        combo_df.insert(0, "MODALITY", modality)
    t_combo = time.time() - t_combo_start
    t_elapsed_total = time.time() - t_pipeline_start
    avg_per_combo = t_elapsed_total / combo_i
    eta_total = avg_per_combo * (n_combos - combo_i)
    print(f"[progress] {modality} x {omics_type} done in {t_combo:.1f}s -- "
          f"{len(all_summaries)} cluster(s) with data, {n_skipped} skipped, "
          f"{len(combo_df)} rows. Overall (this modality): {combo_i}/{n_combos} combos, "
          f"{t_elapsed_total / 60:.1f} min elapsed, ETA ~{eta_total / 60:.1f} min remaining.",
          flush=True)
    return combo_df


def run_modality(modality, omics_types, out_dir, n_clusters):
    """Runs every requested omics type for ONE modality, saving that
    modality's combined result to its own CSV."""
    os.makedirs(out_dir, exist_ok=True)
    n_combos = len(omics_types)
    t_pipeline_start = time.time()
    print(f"\n##### MODALITY: {modality} (omics types: {omics_types}) #####", flush=True)
    if modality in EXCLUDE_CLUSTERS:
        excluded = sorted(EXCLUDE_CLUSTERS[modality])
        print(f"[info] {modality}: excluding {len(excluded)} cluster(s) with n_ldrs<=1: {excluded}", flush=True)
    combo_results = []
    for combo_i, omics_type in enumerate(omics_types, 1):
        combo_results.append(
            process_combo(modality, omics_type, combo_i, n_combos, t_pipeline_start, n_clusters)
        )
    modality_df = pd.concat(combo_results, ignore_index=True) if combo_results else pd.DataFrame()
    out_path = os.path.join(out_dir, f"top_feats_{modality}_allomics_fiberclust.csv")
    modality_df.to_csv(out_path, index=False)
    print(f"[progress] {modality} done in {(time.time() - t_pipeline_start) / 60:.1f} min -- "
          f"wrote {len(modality_df)} rows to {out_path}", flush=True)
    return modality_df


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--modality", required=True, choices=list(eff_num_vox.keys()),
                         help="Which modality to run (also selects its eff_num_vox entry).")
    parser.add_argument("--omics-types", default=None,
                         help=f"Comma-separated subset of omics types, e.g. proteomics,metabolomics. "
                              f"Default: all of {list(eff_num_feature.keys())}.")
    parser.add_argument("--out-dir", default=default_out_dir,
                         help="Directory to write top_feats_<modality>_allomics_fiberclust.csv into.")
    parser.add_argument("--n-clusters", type=int, default=DEFAULT_N_CLUSTERS,
                         help="Number of fiber clusters to loop over (default 430).")
    args = parser.parse_args()
    omics_types = args.omics_types.split(",") if args.omics_types else list(eff_num_feature.keys())
    unknown = set(omics_types) - set(eff_num_feature.keys())
    if unknown:
        parser.error(f"Unknown omics type(s): {unknown} -- valid: {list(eff_num_feature.keys())}")
    run_modality(args.modality, omics_types, args.out_dir, args.n_clusters)


if __name__ == "__main__":
    main()