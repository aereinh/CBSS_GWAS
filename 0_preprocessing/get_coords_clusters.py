#!/usr/bin/env python3
"""
get_coords_clusters.py

Extract 3-D (RAS mm) coordinates for each fiber cluster from a tractogram
and write one coordinate file per cluster.

Each output file contains (n_fibers_in_cluster * n_fiber_pts) rows,
one row per point, with columns: x  y  z.

Run example:
  python get_coords_clusters.py \
      --tractogram  /path/to/atlas_100clust_100points.trk \
      --reference   /path/to/ENIGMA_DTI_FA.nii.gz \
      --cluster-csv /path/to/Fiber_atlas_100points_FullName.csv \
      --cluster-col Clust_Dis30mm \
      --out-dir     /path/to/coords/
"""

import os
import argparse
import numpy as np
import pandas as pd
from dipy.io.streamline import load_tractogram
from dipy.io.stateful_tractogram import StatefulTractogram, Space


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tractogram",   required=True, help="Path to .trk atlas file")
    ap.add_argument("--reference",    required=True, help="Reference anatomy (.nii.gz)")
    ap.add_argument("--cluster-csv",  required=True, help="CSV with per-fiber cluster assignments")
    ap.add_argument("--cluster-col",  default="Clust_Dis30mm",
                    help="Column name for cluster labels (default: Clust_Dis30mm)")
    ap.add_argument("--fiber-pts",    type=int, default=100,
                    help="Points per fiber (default: 100)")
    ap.add_argument("--out-dir",      required=True, help="Output directory for coord files")
    return ap.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # Load tractogram and convert to RAS mm
    tractogram = load_tractogram(args.tractogram, args.reference)
    st = StatefulTractogram(tractogram.streamlines, args.reference, space=Space.RASMM)
    coordinates = np.asarray(st.streamlines)   # shape: (n_fibers, n_pts, 3)
    print(f"Tractogram shape: {coordinates.shape}")

    # Load cluster labels (one entry per fiber)
    clust_data = pd.read_csv(args.cluster_csv)
    if args.cluster_col not in clust_data.columns:
        raise ValueError(
            f"Column '{args.cluster_col}' not found. "
            f"Available: {list(clust_data.columns)}"
        )
    clust_labs = clust_data[args.cluster_col]
    n_clust = clust_labs.nunique()
    print(f"Number of clusters: {n_clust}")

    for i in range(n_clust):
        cluster_id = i + 1
        indices_i  = clust_labs[clust_labs == cluster_id].index
        nfibs_i    = len(indices_i)

        arr_i = np.zeros((args.fiber_pts * nfibs_i, 3))
        for j, fib_j in enumerate(indices_i):
            start = j * args.fiber_pts
            end   = start + args.fiber_pts
            arr_i[start:end, :] = coordinates[fib_j, :, :]

        fname_i = os.path.join(args.out_dir, f"coords_fibersClust_{cluster_id}.txt")
        np.savetxt(fname_i, arr_i, fmt="%.6f", delimiter=" ")

        # Remove trailing newline to match expected format
        with open(fname_i, "rb+") as f:
            f.seek(0, 2)
            if f.tell() > 0:
                f.seek(f.tell() - 1)
                if f.read(1) == b"\n":
                    f.seek(f.tell() - 1)
                    f.truncate()

    print(f"Wrote {n_clust} coordinate files to {args.out_dir}")


if __name__ == "__main__":
    main()
