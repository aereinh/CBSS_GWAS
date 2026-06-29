#!/usr/bin/env python3
"""
agg_vgwas_singlesnp.py

Aggregate voxel-level GWAS Z-scores for a single SNP across all fiber clusters,
producing a whole-brain matrix of shape (n_fibers, n_fiber_pts).

The output .npy file can be loaded for visualization or further analysis
(e.g., mapping Z-scores back to 3D atlas space).

Run example:
  python agg_vgwas_singlesnp.py --snp rs1234567 --out-prefix rs1234567
"""

import os
import argparse
import numpy as np
import pandas as pd
from tqdm import tqdm


# ── Constants (edit to match your setup) ─────────────────────────────────
VGWAS_DIR  = "/path/to/output/voxelgwas/leadsnps"  # Per-cluster VGWAS files
CLUST_FILE = "/path/to/Fiber_atlas_100points_FullName.csv"
CLUST_COL  = "Clust_Dis30mm"
N_FIBERS   = 6090
N_FIBER_PTS = 100
N_CLUST    = 430
# ─────────────────────────────────────────────────────────────────────────


def main():
    parser = argparse.ArgumentParser(
        description="Aggregate per-cluster voxel GWAS Z-scores for a single SNP."
    )
    parser.add_argument("--snp",        required=True, help="SNP ID to extract")
    parser.add_argument("--out-prefix", default=None,
                        help="Output filename prefix (default: SNP ID)")
    args = parser.parse_args()

    snp        = args.snp
    out_prefix = args.out_prefix or snp

    # Load per-fiber cluster labels
    clust_data = pd.read_csv(CLUST_FILE)
    clust_labs = clust_data[CLUST_COL].to_numpy()

    # Output: full (n_fibers × n_pts) Z-score matrix
    z_mat = np.full((N_FIBERS, N_FIBER_PTS), np.nan, dtype=np.float32)

    for cluster_id in tqdm(range(1, N_CLUST + 1), desc="Clusters"):
        fp = os.path.join(VGWAS_DIR, f"UKB_FA_allphases_{cluster_id}.txt")
        if not os.path.exists(fp):
            print(f"[WARN] Missing: {fp}")
            continue

        df = pd.read_csv(fp, sep="\t", usecols=["SNP", "Z"])
        df_snp = df[df["SNP"] == snp]
        if df_snp.empty:
            continue

        fiber_idx    = np.where(clust_labs == cluster_id)[0]
        expected_len = len(fiber_idx) * N_FIBER_PTS

        if df_snp.shape[0] != expected_len:
            raise ValueError(
                f"Cluster {cluster_id}: expected {expected_len} rows for {snp}, "
                f"found {df_snp.shape[0]}."
            )

        z_vals = df_snp["Z"].to_numpy(dtype=np.float32).reshape(len(fiber_idx), N_FIBER_PTS)
        for j, fib_idx in enumerate(fiber_idx):
            z_mat[fib_idx, :] = z_vals[j, :]

        print(f"  Cluster {cluster_id}: filled {len(fiber_idx)} fibers")

    out_path = os.path.join(VGWAS_DIR, f"{out_prefix}_zmat.npy")
    np.save(out_path, z_mat)
    print(f"\nSaved Z-score matrix ({z_mat.shape}) to: {out_path}")
    print(f"Non-NaN entries: {np.sum(~np.isnan(z_mat))}")


if __name__ == "__main__":
    main()
