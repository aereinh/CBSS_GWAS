#!/usr/bin/env python3
"""
step10_select_top_voxels_per_locus.py

For a given tract and lead SNP, identify all voxels (cluster × voxel index)
that are significant in both discovery and replication within a ±250 kb window
around the locus.  Reports per-voxel minimum p-values and SNP counts.

Run example:
  python step10_select_top_voxels_per_locus.py \
      --modality dti_FA --tract SCC --snp rs1234567
"""

import os
import argparse
import pandas as pd
from collections import defaultdict


def main():
    ap = argparse.ArgumentParser(
        description="Per-voxel summary for a replicated discovery locus."
    )
    ap.add_argument("--modality", required=True)
    ap.add_argument("--tract",    required=True, help="Tract name, e.g. SCC")
    ap.add_argument("--snp",      required=True, help="Lead SNP ID")
    ap.add_argument("--window",   type=int, default=250_000,
                    help="Half-window in bp (default: 250000)")
    ap.add_argument("--null-tag", default="1e4",
                    help="Replication file tag (default: 1e4)")
    ap.add_argument("--main-dir", default="/path/to/output")
    ap.add_argument("--cluster-to-tract", default=None,
                    help="Path to cluster_tract_summary.csv")
    args = ap.parse_args()

    modality = args.modality
    tract    = args.tract
    snp      = args.snp
    window   = args.window
    tag      = args.null_tag

    modality_dir = os.path.join(args.main_dir, f"output_{modality}")
    disc_dir     = os.path.join(modality_dir, "post_screen", "disc", "tract")
    disc_vox_dir = os.path.join(modality_dir, "post_screen", "disc", "fibercluster")
    rep_dir      = os.path.join(modality_dir, "replication", "tract")

    # Cluster IDs for this tract
    c2t_file = args.cluster_to_tract or os.path.join(modality_dir, "cluster_tract_summary.csv")
    c2t      = pd.read_csv(c2t_file)
    clusters = c2t.loc[c2t["tract_assign"] == tract, "cluster"].tolist()
    if not clusters:
        raise ValueError(f"No clusters found for tract '{tract}'")
    print(f"Tract {tract}: {len(clusters)} cluster(s)")

    # Get locus position from spatial overlap file
    overlap_file = os.path.join(
        disc_dir,
        f"UKB_{modality}_phase16_Tract_{tract}_clump_replication_spatial_overlap_{tag}.txt"
    )
    if not os.path.exists(overlap_file):
        raise FileNotFoundError(f"Spatial overlap file not found: {overlap_file}")

    clump = pd.read_csv(overlap_file, sep="\t")
    locus_row = clump[clump["snp"] == snp]
    if locus_row.empty:
        raise ValueError(f"SNP '{snp}' not found in {overlap_file}")

    locus = locus_row.iloc[0]
    if locus["overlap_voxels"] == 0:
        print(f"No spatial overlap for locus {snp} at tract {tract}.")
        return

    chr_ = int(locus["chr"])
    pos  = int(locus["pos"])
    start = pos - window
    end   = pos + window

    # Load discovery and replication voxel hits for all clusters
    disc_all, rep_all = [], []
    for cid in clusters:
        d_path = os.path.join(disc_vox_dir, f"UKB_{modality}_phase16_{cid}_sig_vgwas_screen.txt")
        r_path = os.path.join(rep_dir,      f"UKB_{modality}_phase7_{cid}_sig_vgwas_screen_{tag}.txt")

        if os.path.exists(d_path):
            d = pd.read_csv(d_path, sep="\t")[["INDEX", "CHR", "POS", "SNP", "P", "Z"]].copy()
            d["CLUSTER"] = cid
            disc_all.append(d)
        if os.path.exists(r_path):
            r = pd.read_csv(r_path, sep="\t")[["INDEX", "CHR", "POS", "SNP", "P", "Z"]].copy()
            r["CLUSTER"] = cid
            rep_all.append(r)

    if not disc_all or not rep_all:
        print(f"Missing voxel files for tract {tract}.")
        return

    disc = pd.concat(disc_all, ignore_index=True)
    rep  = pd.concat(rep_all,  ignore_index=True)
    disc["CLUSTER_INDEX"] = disc["CLUSTER"].astype(str) + "-" + disc["INDEX"].astype(str)
    rep["CLUSTER_INDEX"]  = rep["CLUSTER"].astype(str)  + "-" + rep["INDEX"].astype(str)

    # Filter to locus window
    disc_locus = disc[(disc["CHR"] == chr_) & (disc["POS"] > start) & (disc["POS"] < end)]
    rep_locus  = rep[(rep["CHR"]  == chr_) & (rep["POS"]  > start) & (rep["POS"]  < end)]

    if disc_locus.empty or rep_locus.empty:
        print(f"No voxel hits within ±{window} bp of locus {snp}.")
        return

    # Voxels present in both
    overlap_idx = sorted(
        set(disc_locus["CLUSTER_INDEX"]).intersection(set(rep_locus["CLUSTER_INDEX"]))
    )
    if not overlap_idx:
        print(f"No overlapping voxels for locus {snp} at tract {tract}.")
        return

    # Per-voxel summary
    disc_grp = disc_locus.groupby("CLUSTER_INDEX", as_index=False).agg(
        disc_min_p  = ("P", "min"),
        disc_n_snps = ("P", "size"),
    )
    rep_grp = rep_locus.groupby("CLUSTER_INDEX", as_index=False).agg(
        rep_min_p  = ("P", "min"),
        rep_n_snps = ("P", "size"),
    )

    out = (
        pd.DataFrame({"CLUSTER_INDEX": overlap_idx})
        .merge(disc_grp, on="CLUSTER_INDEX", how="left")
        .merge(rep_grp,  on="CLUSTER_INDEX", how="left")
    )
    out.insert(0, "tract",    tract)
    out.insert(1, "snp",      snp)
    out.insert(2, "chr",      chr_)
    out.insert(3, "pos",      pos)

    out_file = os.path.join(
        disc_dir,
        f"UKB_{modality}_phase16_Tract_{tract}_locus_{snp}_replicated_voxels.txt"
    )
    out.to_csv(out_file, sep="\t", index=False)
    print(f"Wrote {len(out)} overlapping voxels to: {out_file}")


if __name__ == "__main__":
    main()
