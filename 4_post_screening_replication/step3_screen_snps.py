#!/usr/bin/env python3
"""
step3_screen_snps.py

Screen candidate SNPs using a spatial cluster-size criterion.  For each SNP,
count the number of voxels in the fiber cluster with p < 1e-5.  SNPs whose
voxel count exceeds the (1 - alpha)^eff_num quantile of the null cluster-size
distribution (from step 5d cluster analysis) are retained.

This mirrors the approach of spatially correlated cluster-size testing, using
an empirical null distribution from the discovery bootstrap.

Run example:
  python step3_screen_snps.py --modality dti_FA --roi 42
"""

import os
import argparse
import numpy as np
import pandas as pd


def main(args):
    modality_dir = os.path.join(args.main_dir, f"output_{args.modality}")
    cluster_dir  = args.cluster_dir
    out_dir      = os.path.join(modality_dir, "post_screen", "disc", "fibercluster")

    data_file = os.path.join(
        out_dir, f"UKB_{args.modality}_phase16_{args.roi}_vgwas_1e5.txt"
    )
    if not os.path.exists(data_file):
        print(f"Skipping ROI {args.roi}: file not found ({data_file})")
        return

    data = pd.read_csv(data_file, sep="\t")
    data_count = data["SNP"].value_counts()
    print(f"{len(data_count)} unique SNPs at p<1e-5.")

    # Load null cluster-size distribution (multiple bootstrap batches)
    null_chunks = []
    for batch in range(1, args.n_batches + 1):
        null_file = os.path.join(
            cluster_dir, f"UKB_{args.modality}_phase16_{args.roi}_{batch}.txt"
        )
        if os.path.exists(null_file):
            null_chunks.append(np.loadtxt(null_file))
        else:
            print(f"  [WARN] Missing null file: {null_file}")

    if not null_chunks:
        print(f"No null distribution files found for ROI {args.roi}. Skipping.")
        return

    null_count = np.concatenate([np.atleast_1d(x) for x in null_chunks])
    print(f"Null distribution size: {null_count.shape[0]}")

    # Effective number of voxels for this cluster
    eff_num_file = os.path.join(modality_dir, "eff_num_voxels.csv")
    eff_num_df   = pd.read_csv(eff_num_file)
    eff_num_roi  = eff_num_df["eff_num"].iloc[int(args.roi) - 1]

    # Cluster-size threshold: control FWER at alpha=1e-5 across eff_num tests
    quantile = np.quantile(null_count, (1 - 1e-5) ** eff_num_roi)
    print(f"Cluster-size threshold: {quantile:.2f}  "
          f"(null values > threshold: {np.sum(null_count > quantile)})")

    count_screen = data_count[data_count > quantile]
    print(f"{len(count_screen)} SNPs pass cluster-size screening.")

    # Re-read genome-wide hits (p < genome-wide threshold) and filter to screened SNPs
    vgwas_file = os.path.join(
        modality_dir, "voxelgwas",
        f"UKB_{args.modality}_phase16_{args.roi}_all_voxels_whole_genome.txt"
    )
    sig_data = pd.read_csv(vgwas_file, sep="\t")
    sig_data  = sig_data[sig_data["P"] < args.sig_thresh]
    sig_data_screen = sig_data[sig_data["SNP"].isin(count_screen.index)]
    print(f"{sig_data_screen.shape[0]} voxel-SNP associations after screening.")

    out_file = os.path.join(
        out_dir, f"UKB_{args.modality}_phase16_{args.roi}_sig_vgwas_screen.txt"
    )
    sig_data_screen.to_csv(out_file, sep="\t", index=False)
    print(f"Wrote: {out_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Screen voxel GWAS hits using cluster-size null distribution."
    )
    parser.add_argument("--modality",    required=True, help="DTI modality, e.g. dti_FA")
    parser.add_argument("--roi",         required=True, help="Fiber cluster ID")
    parser.add_argument("--main-dir",    default="/path/to/output",
                        help="Root output directory")
    parser.add_argument("--cluster-dir", default="/path/to/output/cluster",
                        help="Directory with null cluster-size txt files")
    parser.add_argument("--sig-thresh",  type=float, default=9.12e-12,
                        help="Genome-wide significance threshold (default: 9.12e-12)")
    parser.add_argument("--n-batches",   type=int, default=5,
                        help="Number of bootstrap batches to load (default: 5)")
    main(parser.parse_args())
