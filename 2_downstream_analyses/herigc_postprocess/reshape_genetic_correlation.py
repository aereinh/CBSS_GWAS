#!/usr/bin/env python3
"""
reshape_genetic_correlation.py

Assemble the per-cluster genetic-correlation output of HEIG
(genetic_correlation.slurm) into fibers x points matrices per trait
(6,090 atlas fibers x 100 points): the estimate, its standard error and its
P-value (GC_<trait>_matrix.txt, GC_SE_<trait>_matrix.txt, GC_P_<trait>_matrix.txt).

Usage:
    python reshape_genetic_correlation.py \
        --gc-dir      /path/to/output/herigc/ctgc \
        --modality    dti_FA \
        --cluster-csv /path/to/Fiber_atlas_100points_FullName.csv \
        --out-dir     /path/to/output/result_data/dti_FA
"""
import argparse
import os
from glob import glob

import numpy as np
import pandas as pd

TRAITS = ['ADHD', 'ALZ2', 'ASD', 'bip', 'CP', 'depressed_affect', 'depression', 'EA', 'insomnia',
          'intelligence', 'mdd', 'neuroticism', 'SCZ', 'worry']
COLUMNS = ['GC', 'GC_SE', 'GC_P']   # columns of the HEIG output


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--gc-dir", required=True, help="Directory with one sub-directory per trait")
    ap.add_argument("--modality", default="dti_FA")
    ap.add_argument("--cluster-csv", required=True, help="Fiber atlas table with the cluster label of each fiber")
    ap.add_argument("--cluster-col", default="Clust_Dis30mm")
    ap.add_argument("--traits", default=",".join(TRAITS), help="Comma-separated trait codes")
    ap.add_argument("--points-per-fiber", type=int, default=100)
    ap.add_argument("--out-dir", required=True)
    args = ap.parse_args()

    clust_labs = pd.read_csv(args.cluster_csv)[args.cluster_col]
    n_clust = len(np.unique(clust_labs))
    npt = args.points_per_fiber
    os.makedirs(args.out_dir, exist_ok=True)

    for trait in args.traits.split(","):
        mats = {c: np.zeros((len(clust_labs), npt), dtype=np.float64) for c in COLUMNS}
        for i in range(n_clust):
            files = glob(os.path.join(args.gc_dir, trait, f"UKB_{args.modality}_phase16_{i + 1}_*_gc.txt"))
            if not files:
                print(f"Missing file for cluster {i + 1} of {trait}")
                continue
            res = pd.read_csv(files[0], sep=r"\s+")
            for j, fiber_idx in enumerate(clust_labs[clust_labs == (i + 1)].index):
                for c in COLUMNS:
                    mats[c][fiber_idx, :] = res[c].iloc[j * npt:(j + 1) * npt].values
        for c in COLUMNS:
            out = os.path.join(args.out_dir, f"{c}_{trait}_matrix.txt")
            np.savetxt(out, mats[c], fmt="%.6e")
            print(f"Saved: {out}")


if __name__ == "__main__":
    main()
