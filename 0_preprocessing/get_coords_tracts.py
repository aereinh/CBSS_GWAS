#!/usr/bin/env python3
"""
get_coords_tracts.py

Aggregate cluster-level coordinate files up to tract-level by reading a
cluster-to-tract assignment CSV and combining all fibers within each tract.

Each output file contains (n_fibers_in_tract * n_fiber_pts) rows:  x  y  z

Run example:
  python get_coords_tracts.py \
      --tractogram       /path/to/atlas_100clust_100points.trk \
      --reference        /path/to/ENIGMA_DTI_FA.nii.gz \
      --cluster-csv      /path/to/Fiber_atlas_100points_FullName.csv \
      --cluster-col      Clust_Dis30mm \
      --cluster2tract    /path/to/cluster_tract_summary.csv \
      --out-dir          /path/to/coords/tracts/
"""

import os
import argparse
import numpy as np
import pandas as pd
from dipy.io.streamline import load_tractogram
from dipy.io.stateful_tractogram import StatefulTractogram, Space


def safe_name(s: str) -> str:
    """Make a string safe for use in filenames."""
    return "".join(c if c.isalnum() or c in ("_", "-") else "_" for c in str(s))


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tractogram",    required=True, help="Path to .trk atlas file")
    ap.add_argument("--reference",     required=True, help="Reference anatomy (.nii.gz)")
    ap.add_argument("--cluster-csv",   required=True, help="CSV with per-fiber cluster assignments")
    ap.add_argument("--cluster-col",   default="Clust_Dis30mm",
                    help="Column name for cluster labels (default: Clust_Dis30mm)")
    ap.add_argument("--cluster2tract", required=True,
                    help="CSV mapping cluster IDs to tract names (columns: cluster, tract_assign)")
    ap.add_argument("--out-dir",       required=True, help="Output directory for tract coord files")
    return ap.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # Load tractogram and convert to RAS mm
    tractogram   = load_tractogram(args.tractogram, args.reference)
    st           = StatefulTractogram(tractogram.streamlines, args.reference, space=Space.RASMM)
    coordinates  = np.asarray(st.streamlines)   # shape: (n_fibers, n_pts, 3)
    print(f"Tractogram shape: {coordinates.shape}", flush=True)

    # Per-fiber cluster labels
    clust_data   = pd.read_csv(args.cluster_csv)
    clust_labs   = clust_data[args.cluster_col].astype(int).to_numpy()
    print(f"Fiber cluster labels: {len(clust_labs)}", flush=True)

    # Cluster -> tract mapping
    clust2tract  = pd.read_csv(args.cluster2tract)
    clust2tract["cluster"]      = clust2tract["cluster"].astype(int)
    clust2tract["tract_assign"] = clust2tract["tract_assign"].fillna("UNASSIGNED")

    tract_to_clusters = (
        clust2tract.groupby("tract_assign")["cluster"].apply(list).to_dict()
    )
    print(f"Number of tract groups: {len(tract_to_clusters)}", flush=True)

    for tract, clusters in tract_to_clusters.items():
        fiber_idx = np.where(np.isin(clust_labs, clusters))[0]
        nfibs = len(fiber_idx)

        if nfibs == 0:
            print(f"[WARN] '{tract}': 0 fibers for clusters={clusters}", flush=True)
            continue

        arr  = coordinates[fiber_idx, :, :].reshape(-1, 3)
        fname = os.path.join(args.out_dir, f"coords_tract_{safe_name(tract)}.txt")
        np.savetxt(fname, arr, fmt="%.6f", delimiter=" ")

        # Remove trailing newline
        with open(fname, "rb+") as f:
            f.seek(0, 2)
            if f.tell() > 0:
                f.seek(f.tell() - 1)
                if f.read(1) == b"\n":
                    f.seek(f.tell() - 1)
                    f.truncate()

        print(f"[OK] {tract}: {len(clusters)} cluster(s), {nfibs} fibers → {fname}", flush=True)

    print("Done.", flush=True)


if __name__ == "__main__":
    main()
