#!/usr/bin/env python3
"""
qc.py

Scan raw fiber-cluster .mat files for outlier values, producing:
  1. A per-fiber CSV listing every fiber that exceeds the threshold.
  2. A per-cluster summary CSV with counts and max values.

Run example:
  python qc.py \
      --phase Phase1 \
      --modality dti_MD \
      --work-root /path/to/OtherModality/ \
      --cluster-csv /path/to/Fiber_atlas_100points_FullName.csv \
      --out-prefix /path/to/output/Phase1_dti_MD_badscan
"""

import os
import argparse
import numpy as np
import pandas as pd
import scipy.io as sio


def parse_args():
    ap = argparse.ArgumentParser()
    ap.add_argument("--phase",      required=True,  help="Phase label, e.g. Phase1")
    ap.add_argument("--modality",   required=True,  help="DTI modality, e.g. dti_MD")
    ap.add_argument("--work-root",  required=True,
                    help="Root directory containing modality subfolders")
    ap.add_argument("--cluster-csv", required=True,
                    help="CSV with per-fiber cluster assignments")
    ap.add_argument("--cluster-col", default="Clust_Dis30mm",
                    help="Column name for cluster labels (default: Clust_Dis30mm)")
    ap.add_argument("--threshold",   type=float, default=100.0,
                    help="Outlier threshold; fibers with max > threshold are flagged (default: 100)")
    ap.add_argument("--n-files",     type=int,   default=609,
                    help="Number of Aligned .mat files expected (default: 609)")
    ap.add_argument("--fibers-per-file", type=int, default=10,
                    help="Fibers per .mat file (default: 10)")
    ap.add_argument("--fiber-pts",   type=int,   default=100,
                    help="Points per fiber (default: 100)")
    ap.add_argument("--out-prefix",  required=True,
                    help="Output file prefix, e.g. /path/to/Phase1_dti_MD_badscan")
    return ap.parse_args()


def main():
    args = parse_args()

    data_dir = os.path.join(
        args.work_root, args.modality, f"step2_{args.modality}_UKB_{args.phase}"
    )
    if not os.path.isdir(data_dir):
        raise FileNotFoundError(f"Data directory not found: {data_dir}")

    # Load per-fiber cluster labels
    clust_df   = pd.read_csv(args.cluster_csv)
    if args.cluster_col not in clust_df.columns:
        raise ValueError(
            f"Column '{args.cluster_col}' not found in {args.cluster_csv}. "
            f"Available columns: {list(clust_df.columns)}"
        )
    clust_labs = clust_df[args.cluster_col].to_numpy()
    n_fibers_total = len(clust_labs)

    expected_total = args.n_files * args.fibers_per_file
    if expected_total != n_fibers_total:
        print(
            f"WARNING: n_files * fibers_per_file = {expected_total} but "
            f"cluster-label length = {n_fibers_total}. Proceeding anyway."
        )

    bad_rows = []

    for file_index in range(1, args.n_files + 1):
        mat_path = os.path.join(data_dir, f"{file_index}_Aligned.mat")
        if not os.path.exists(mat_path):
            print(f"[MISS] {mat_path}")
            continue

        try:
            md = sio.loadmat(mat_path)
        except Exception as e:
            print(f"[ERR ] Cannot load {mat_path}: {e}")
            continue

        if "f_aligned" not in md:
            print(f"[ERR ] {mat_path} missing key 'f_aligned'")
            continue

        fa = md["f_aligned"]
        shape = fa.shape

        # Identify the fiber axis (should be size == fibers_per_file)
        fiber_axes = [ax for ax, s in enumerate(shape) if s == args.fibers_per_file]
        if len(fiber_axes) != 1:
            print(f"[WARN] {mat_path}: cannot uniquely identify fiber axis from shape {shape}. Skipping.")
            continue
        fiber_axis = fiber_axes[0]

        # Move fiber axis to position 0 and compute per-fiber max
        fa2 = np.moveaxis(fa, fiber_axis, 0)   # shape: (n_fibers, ...)
        fiber_max = np.max(fa2, axis=tuple(range(1, fa2.ndim)))

        bad_local = np.where(fiber_max > args.threshold)[0]
        if bad_local.size > 0:
            base_global = (file_index - 1) * args.fibers_per_file
            for local_fib in bad_local.tolist():
                global_fib = base_global + local_fib
                cl = clust_labs[global_fib] if global_fib < n_fibers_total else np.nan
                bad_rows.append({
                    "file_index":           file_index,
                    "mat_path":             mat_path,
                    "fiber_local_0based":   local_fib,
                    "fiber_global_0based":  global_fib,
                    "cluster":              cl,
                    "max_value":            float(fiber_max[local_fib]),
                })

        if file_index % 25 == 0:
            print(f"[{file_index}/{args.n_files}] scanned … bad fibers so far: {len(bad_rows)}")

    bad_df     = pd.DataFrame(bad_rows)
    bad_csv    = f"{args.out_prefix}_bad_fibers.csv"
    summary_csv = f"{args.out_prefix}_bad_clusters.csv"

    if len(bad_df) == 0:
        print("No fibers exceeded threshold.")
        bad_df.to_csv(bad_csv, index=False)
        pd.DataFrame(
            columns=["cluster", "n_bad_fibers", "n_bad_files", "max_over_cluster"]
        ).to_csv(summary_csv, index=False)
    else:
        bad_df.to_csv(bad_csv, index=False)

        bad_df["cluster"] = (
            bad_df["cluster"].astype("object").where(~bad_df["cluster"].isna(), other="NA")
        )
        summ = (
            bad_df
            .groupby("cluster", dropna=False)
            .agg(
                n_bad_fibers      = ("max_value",   "size"),
                n_bad_files       = ("file_index",  pd.Series.nunique),
                max_over_cluster  = ("max_value",   "max"),
            )
            .reset_index()
            .sort_values(["n_bad_fibers", "max_over_cluster"], ascending=[False, False])
        )
        summ.to_csv(summary_csv, index=False)

        print(f"\nDone.")
        print(f"  Bad fibers:   {bad_csv}  ({len(bad_df)} rows)")
        print(f"  Bad clusters: {summary_csv}  ({summ.shape[0]} clusters)")
        print("\nTop affected clusters:")
        print(summ.head(15).to_string(index=False))

    print(f"\nWrote: {bad_csv}")
    print(f"Wrote: {summary_csv}")


if __name__ == "__main__":
    main()
