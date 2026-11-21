#!/usr/bin/env python3
"""
step9_spatial_overlap_one_tract.py

Compute spatial (voxel-level) overlap between discovery and replication GWAS
hits for each lead locus in a tract.

For each discovery locus (from step 8), counts:
  - disc_voxels:     discovery voxels with a significant hit in ±window
  - rep_voxels:      replication voxels with a significant hit in ±window
  - overlap_voxels:  voxels significant in both (same cluster × voxel index)
  - overlap_percent: overlap_voxels / disc_voxels

Run example:
  python step9_spatial_overlap_one_tract.py --modality dti_FA --tract SCC
"""

import os
import argparse
import pandas as pd


def main():
    ap = argparse.ArgumentParser(
        description="Spatial voxel overlap between discovery and replication for one tract."
    )
    ap.add_argument("--modality",  required=True)
    ap.add_argument("--tract",     required=True, help="Tract name, e.g. SCC")
    ap.add_argument("--window",    type=int, default=250_000,
                    help="Half-window in bp (default: 250000)")
    ap.add_argument("--null-tag",  default="1e4",
                    help="Tag used in step7 filenames (default: 1e4)")
    ap.add_argument("--main-dir",  default="/path/to/output")
    ap.add_argument("--cluster-to-tract", default=None,
                    help="Path to cluster_tract_summary.csv")
    args = ap.parse_args()

    modality = args.modality
    tract    = args.tract
    window   = args.window
    tag      = args.null_tag

    modality_dir  = os.path.join(args.main_dir, f"output_{modality}")
    disc_dir      = os.path.join(modality_dir, "post_screen", "disc", "tract")
    disc_vox_dir  = os.path.join(modality_dir, "post_screen", "disc", "fibercluster")
    rep_dir       = os.path.join(modality_dir, "replication", "tract")

    # Load cluster assignments
    c2t_file = args.cluster_to_tract or os.path.join(modality_dir, "cluster_tract_summary.csv")
    c2t      = pd.read_csv(c2t_file)
    clusters = c2t.loc[c2t["tract_assign"] == tract, "cluster"].tolist()
    if not clusters:
        raise ValueError(f"No clusters found for tract '{tract}'")
    print(f"Tract {tract}: {len(clusters)} cluster(s)")

    # Load replication result (step 8) to get locus list
    clump_file = os.path.join(
        disc_dir,
        f"UKB_{modality}_phase16_Tract_{tract}_clump_replication_by_phase7_{tag}.txt"
    )
    if not os.path.exists(clump_file):
        raise FileNotFoundError(f"Replication file not found: {clump_file}")

    clump = pd.read_csv(clump_file, sep="\t")
    clump["CHR"] = clump["CHR"].astype(int)
    clump["BP"]  = clump["BP"].astype(int)

    # Aggregate discovery and replication voxel hits
    disc_all, rep_all = [], []
    for cid in clusters:
        disc_path = os.path.join(
            disc_vox_dir, f"UKB_{modality}_phase16_{cid}_sig_vgwas_screen.txt"
        )
        rep_path = os.path.join(
            rep_dir, f"UKB_{modality}_phase7_{cid}_sig_vgwas_screen_{tag}.txt"
        )
        if os.path.exists(disc_path):
            d = pd.read_csv(disc_path, sep="\t")[["INDEX", "CHR", "POS", "SNP"]].copy()
            d["CLUSTER"] = cid
            d["CHR"] = d["CHR"].astype(int)
            d["POS"] = d["POS"].astype(int)
            disc_all.append(d)
        if os.path.exists(rep_path):
            r = pd.read_csv(rep_path, sep="\t")[["INDEX", "CHR", "POS", "SNP"]].copy()
            r["CLUSTER"] = cid
            r["CHR"] = r["CHR"].astype(int)
            r["POS"] = r["POS"].astype(int)
            rep_all.append(r)

    if not disc_all:
        raise RuntimeError(f"No discovery voxel files found for tract {tract}")
    if not rep_all:
        raise RuntimeError(f"No replication voxel files found for tract {tract}")

    disc = pd.concat(disc_all, ignore_index=True)
    rep  = pd.concat(rep_all,  ignore_index=True)

    # Unique voxel identifier = cluster ID + voxel index
    disc["CLUSTER_INDEX"] = disc["CLUSTER"].astype(str) + "-" + disc["INDEX"].astype(str)
    rep["CLUSTER_INDEX"]  = rep["CLUSTER"].astype(str)  + "-" + rep["INDEX"].astype(str)

    print(f"Discovery rows: {len(disc)}  |  Replication rows: {len(rep)}")

    # Per-locus spatial overlap
    rows = []
    for _, locus in clump.iterrows():
        chr_  = int(locus["CHR"])
        pos   = int(locus["BP"])
        start = pos - window
        end   = pos + window

        disc_idx = set(disc.loc[
            (disc["CHR"] == chr_) & (disc["POS"] > start) & (disc["POS"] < end),
            "CLUSTER_INDEX"
        ])
        rep_idx = set(rep.loc[
            (rep["CHR"] == chr_) & (rep["POS"] > start) & (rep["POS"] < end),
            "CLUSTER_INDEX"
        ])

        overlap = disc_idx & rep_idx
        disc_n  = len(disc_idx)
        ov_pct  = round(len(overlap) / disc_n, 3) if disc_n > 0 else 0.0

        rows.append({
            "chr":             chr_,
            "snp":             locus.get("SNP"),
            "pos":             pos,
            "p":               locus.get("P"),
            "replication":     locus.get("replication"),
            "disc_voxels":     disc_n,
            "rep_voxels":      len(rep_idx),
            "overlap_voxels":  len(overlap),
            "overlap_percent": ov_pct,
        })

    output = pd.DataFrame(rows)
    print(f"Loci with >0% overlap:  {(output['overlap_percent'] > 0).sum()}")
    print(f"Loci with >50% overlap: {(output['overlap_percent'] > 0.5).sum()}")

    out_file = os.path.join(
        disc_dir,
        f"UKB_{modality}_phase16_Tract_{tract}_clump_replication_spatial_overlap_{tag}.txt"
    )
    output.to_csv(out_file, sep="\t", index=False)
    print(f"Wrote: {out_file}")


if __name__ == "__main__":
    main()
