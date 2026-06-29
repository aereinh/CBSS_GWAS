#!/usr/bin/env python3
"""
step1_get_unique_snps.py

Extract unique SNP IDs from genome-wide voxel GWAS output that pass the
Bonferroni-corrected significance threshold.  The resulting SNP list is used
in step 2 to reconstruct full voxel-level associations for those SNPs only,
saving substantial compute time.

Run example:
  python step1_get_unique_snps.py --modality dti_FA --roi 42
"""

import os
import argparse
import pandas as pd


# Significance threshold: 5e-8 / effective number of independent voxels
SIG_THRESH = 9.12e-12   # Adjust to match your voxelgwas_screening threshold


def main(args):
    vgwas_dir = os.path.join(args.main_dir, f"output_{args.modality}", "voxelgwas")
    out_dir   = os.path.join(args.main_dir, f"output_{args.modality}",
                             "post_screen", "disc", "fibercluster")
    os.makedirs(out_dir, exist_ok=True)

    vgwas_file = os.path.join(
        vgwas_dir,
        f"UKB_{args.modality}_allphases_{args.roi}_all_voxels_whole_genome.txt"
    )
    if not os.path.exists(vgwas_file):
        raise FileNotFoundError(f"Voxel GWAS file not found: {vgwas_file}")

    data = pd.read_csv(vgwas_file, sep="\t")
    data = data[data["P"] < SIG_THRESH]

    unique_snps = data["SNP"].unique()
    if len(unique_snps) == 0:
        print(f"No significant SNPs found for ROI {args.roi} (threshold={SIG_THRESH}).")
        return

    out_file = os.path.join(
        out_dir, f"UKB_{args.modality}_allphases_{args.roi}_sig_unique_snps.txt"
    )
    pd.DataFrame({"SNP": unique_snps}).to_csv(out_file, header=None, index=None)
    print(f"Wrote {len(unique_snps)} unique SNPs to: {out_file}")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(
        description="Extract unique significant SNPs from voxel GWAS output."
    )
    parser.add_argument("--modality", required=True, help="DTI modality, e.g. dti_FA")
    parser.add_argument("--roi",      required=True, help="Fiber cluster ID")
    parser.add_argument("--main-dir", default="/path/to/output",
                        help="Root output directory containing output_<modality>/ subdirs")
    main(parser.parse_args())
