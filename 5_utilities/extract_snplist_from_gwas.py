#!/usr/bin/env python3
"""
extract_snplist_from_gwas.py

Extract a unique list of SNP IDs from a tab-delimited GWAS results file
(e.g., screened voxel GWAS output).  The SNP list can be used to extract
genotype data or run targeted voxel GWAS with --extract.

Run example:
  python extract_snplist_from_gwas.py \
      --infile  UKB_dti_FA_allphases_42_sig_vgwas_screen.txt \
      --outfile snplist_cluster42.txt
"""

import os
import argparse
import pandas as pd


def main():
    parser = argparse.ArgumentParser(
        description="Extract unique SNP IDs from a GWAS results file."
    )
    parser.add_argument("--infile",  required=True,
                        help="Tab-delimited GWAS results file with a 'SNP' column")
    parser.add_argument("--outfile", required=True,
                        help="Output text file (one SNP per line, no header)")
    args = parser.parse_args()

    if not os.path.exists(args.infile):
        raise FileNotFoundError(f"Input file not found: {args.infile}")

    df = pd.read_csv(args.infile, sep="\t")
    if "SNP" not in df.columns:
        raise ValueError(
            f"Column 'SNP' not found in {args.infile}. "
            f"Available columns: {list(df.columns)}"
        )

    unique_snps = df["SNP"].dropna().unique()
    pd.Series(unique_snps).to_csv(args.outfile, index=False, header=False)
    print(f"Wrote {len(unique_snps)} unique SNPs to: {args.outfile}")


if __name__ == "__main__":
    main()
