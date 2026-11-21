#!/usr/bin/env python3
"""
reshape_heritability.py

Assemble the per-cluster heritability output of HEIG (heritability.slurm) into
one fibers x points matrix (6,090 atlas fibers x 100 points), the layout used
for tract- and network-level summaries and for the brain maps.

Usage:
    python reshape_heritability.py \
        --heri-dir    /path/to/output/herigc/heri \
        --modality    dti_FA \
        --cluster-csv /path/to/Fiber_atlas_100points_FullName.csv \
        --out         /path/to/output/herigc/heri/heri_matrix.txt
"""
import argparse
import os

import numpy as np
import pandas as pd


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--heri-dir", required=True, help="Directory with UKB_<modality>_phase16_<cluster>_heri.txt")
    ap.add_argument("--modality", default="dti_FA")
    ap.add_argument("--cluster-csv", required=True, help="Fiber atlas table with the cluster label of each fiber")
    ap.add_argument("--cluster-col", default="Clust_Dis30mm")
    ap.add_argument("--column", default="H2", help="Column of the HEIG output to reshape (H2, SE or P)")
    ap.add_argument("--points-per-fiber", type=int, default=100)
    ap.add_argument("--out", required=True)
    args = ap.parse_args()

    clust_labs = pd.read_csv(args.cluster_csv)[args.cluster_col]
    n_clust = len(np.unique(clust_labs))
    npt = args.points_per_fiber
    mat = np.zeros((len(clust_labs), npt), dtype=np.float32)

    for i in range(n_clust):
        res = pd.read_csv(os.path.join(args.heri_dir, f"UKB_{args.modality}_phase16_{i + 1}_heri.txt"), sep=r"\s+")
        # rows of the HEIG output follow the fibers of the cluster in atlas order, 100 points each
        for j, fiber_idx in enumerate(clust_labs[clust_labs == (i + 1)].index):
            mat[fiber_idx, :] = res[args.column].iloc[j * npt:(j + 1) * npt].values

    np.savetxt(args.out, mat, fmt="%.6f")
    print(f"Saved: {args.out}")


if __name__ == "__main__":
    main()
