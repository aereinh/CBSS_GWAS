#!/usr/bin/env python3
"""
remove_outliers.py

Replace outlier values in a fibersClust whitespace-delimited text file with NA.
The first two columns (IID, FID) are left untouched; all subsequent numeric
columns are checked against an absolute-value threshold.

Run example:
  python remove_outliers.py \
      --in-txt  fibersClust_42.txt \
      --out-txt fibersClust_42_clean.txt \
      --threshold 100
"""

import argparse
import numpy as np


def main():
    parser = argparse.ArgumentParser(
        description="Replace outliers in a fibersClust .txt file with NA."
    )
    parser.add_argument("--in-txt",    required=True, help="Input fibersClust_*.txt file")
    parser.add_argument("--out-txt",   required=True, help="Output (cleaned) file path")
    parser.add_argument("--threshold", type=float, required=True,
                        help="Absolute-value threshold; values exceeding this are set to NA")
    args = parser.parse_args()

    thresh = args.threshold

    with open(args.in_txt, "r") as fin, open(args.out_txt, "w") as fout:
        # Preserve header line unchanged
        header = fin.readline()
        fout.write(header)

        for line in fin:
            parts = line.rstrip().split()
            if len(parts) <= 2:
                fout.write(line)
                continue

            # Keep the first two ID columns as-is
            out = parts[:2]

            for v in parts[2:]:
                try:
                    x = float(v)
                    out.append("NA" if (not np.isfinite(x)) or abs(x) > thresh else v)
                except ValueError:
                    out.append("NA")

            fout.write(" ".join(out) + "\n")


if __name__ == "__main__":
    main()
