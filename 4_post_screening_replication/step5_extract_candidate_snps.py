#!/usr/bin/env python3
"""
step5_extract_candidate_snps.py

For each tract, extract all SNPs within ±250 kb of every LD-clumped lead
locus (step 4b output).  These candidate SNPs are used in step 6 to run
targeted voxel GWAS in the validation cohort.

Run example:
  python step5_extract_candidate_snps.py --modality dti_FA --roi SCC
"""

import os
import argparse
import pandas as pd


def main(args):
    main_dir     = os.path.join(args.output_root, f"output_{args.modality}")
    sumstats_dir = args.sumstats_dir
    tract_dir    = os.path.join(main_dir, "post_screen", "disc", "tract")

    # Stage-3 clumped loci for this tract
    clump_file = os.path.join(
        tract_dir,
        f"UKB_{args.modality}_phase16_Tract_{args.roi}_clump_r2_0.1_step3.clumped"
    )
    if not os.path.exists(clump_file):
        raise FileNotFoundError(f"Clump file not found: {clump_file}")

    clump = pd.read_csv(clump_file, sep=r"\s+")

    # Reference SNP info file (use any large cluster's .snpinfo as the universe)
    snpinfo_file = os.path.join(
        sumstats_dir, f"UKB_{args.modality}_phase16_{args.snpinfo_cluster}.snpinfo"
    )
    if not os.path.exists(snpinfo_file):
        raise FileNotFoundError(f"SNP info file not found: {snpinfo_file}")

    snpinfo = pd.read_csv(snpinfo_file, sep="\t")

    # Collect all SNPs within ±window bp of each lead locus
    candidate_snps = set()
    for _, row in clump.iterrows():
        window_snps = snpinfo.loc[
            (snpinfo["CHR"] == row["CHR"]) &
            (snpinfo["POS"] >= row["BP"] - args.window) &
            (snpinfo["POS"] <= row["BP"] + args.window),
            "SNP"
        ]
        candidate_snps.update(window_snps.tolist())

    out_file = os.path.join(
        tract_dir,
        f"UKB_{args.modality}_phase16_{args.roi}_candidate_sig_snps.txt"
    )
    with open(out_file, "w") as f:
        for snp in candidate_snps:
            f.write(snp + "\n")

    print(f"Wrote {len(candidate_snps)} candidate SNPs to: {out_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract candidate SNPs in ±window of each lead locus per tract."
    )
    parser.add_argument("--modality",        required=True, help="DTI modality, e.g. dti_FA")
    parser.add_argument("--roi",             required=True, help="Tract name, e.g. SCC")
    parser.add_argument("--output-root",     default="/path/to/output",
                        help="Root output directory")
    parser.add_argument("--sumstats-dir",    default="/path/to/output/sumstats",
                        help="Directory containing .snpinfo files")
    parser.add_argument("--snpinfo-cluster", default="100",
                        help="Cluster ID whose .snpinfo to use as SNP universe (default: 100)")
    parser.add_argument("--window",          type=int, default=250000,
                        help="Half-window size in bp around each lead SNP (default: 250000)")
    main(parser.parse_args())
