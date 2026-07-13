#!/usr/bin/env python3
"""
step1_screen_ldr_sumstats.py

For a single fiber cluster, extract all LDR associations for SNPs that are
genome-wide significant (p < pval_thresh) in any phenotype (LDR index).
Reads HEIG binary sumstats (.sumstats + .snpinfo), computes Z² thresholds
from the chi-squared distribution, and writes a flat TSV of significant rows.
Also writes a SNP list (for use with PLINK --extract in the next step).

Run example:
  python step1_screen_ldr_sumstats.py \
      --cluster-id 42 \
      --prefix /path/to/sumstats/UKB_dti_FA_allphases_42 \
      --output /path/to/screened/UKB_dti_FA_allphases_42.txt
"""

import os
import math
import argparse
import numpy as np
import pandas as pd
import h5py
from scipy.stats import chi2, norm


# ── GWAS reader (HEIG .sumstats format) ─────────────────────────────────────

def read_sumstats(prefix: str):
    snpinfo_path = f"{prefix}.snpinfo"
    sumstats_path = f"{prefix}.sumstats"

    if not os.path.exists(snpinfo_path) or not os.path.exists(sumstats_path):
        raise FileNotFoundError(
            f"Missing .sumstats or .snpinfo at prefix: {prefix}"
        )

    fh = h5py.File(sumstats_path, "r")
    snpinfo = pd.read_csv(snpinfo_path, sep="\t", engine="pyarrow")

    if snpinfo.shape[0] != fh.attrs["n_snps"]:
        raise ValueError("Mismatched SNP count between .snpinfo and .sumstats")

    return GWAS(fh, snpinfo)


class GWAS:
    def __init__(self, fh, snpinfo):
        self.n_snps = fh.attrs["n_snps"]
        self.n_gwas = fh.attrs["n_gwas"]
        self.n_blocks = fh.attrs["n_blocks"]
        self.fh = fh
        self.snpinfo = snpinfo

    def close(self):
        self.fh.close()

    def _n_blocks_for(self, n_gwas: int) -> int:
        return math.ceil(n_gwas / 20)

    def read_z_scores(self, snp_idxs: np.ndarray) -> np.ndarray:
        """Read full Z-score matrix for selected SNP indices.
        Returns array of shape (len(snp_idxs), n_gwas).
        """
        n_blocks = self._n_blocks_for(self.n_gwas)
        z = np.zeros((len(snp_idxs), self.n_gwas), dtype=np.float32)
        remaining = self.n_gwas

        for blk in range(n_blocks):
            n_cols = min(20, remaining)
            start = blk * 20
            print(f"  Reading Z block {blk} [{start}:{start + n_cols}]")
            z[:, start : start + n_cols] = self.fh[f"z{blk}"][:][snp_idxs, :n_cols]
            remaining -= n_cols

        return z

    def read_betas(self, snp_idxs: np.ndarray) -> np.ndarray:
        """Read BETA matrix for selected SNP indices (block generator)."""
        n_blocks = self._n_blocks_for(self.n_gwas)
        remaining = self.n_gwas
        beta = np.zeros((len(snp_idxs), self.n_gwas), dtype=np.float32)

        for blk in range(n_blocks):
            n_cols = min(20, remaining)
            start = blk * 20
            print(f"  Reading BETA block {blk}")
            beta[:, start : start + n_cols] = self.fh[f"beta{blk}"][:][snp_idxs, :n_cols]
            remaining -= n_cols

        return beta


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(
        description="Screen HEIG LDR sumstats for genome-wide significant SNPs."
    )
    parser.add_argument("--cluster-id", type=int, required=True,
                        help="Fiber cluster ID (used only for logging)")
    parser.add_argument("--prefix", required=True,
                        help="Path prefix to .sumstats and .snpinfo (no extension)")
    parser.add_argument("--output", required=True,
                        help="Output TSV file path")
    parser.add_argument("--pval-thresh", type=float, default=5e-8,
                        help="P-value threshold for significance (default: 5e-8)")
    args = parser.parse_args()

    print(f"Cluster {args.cluster_id}: reading GWAS data from {args.prefix}")
    gwas = read_sumstats(args.prefix)

    z2_thresh = chi2.ppf(1 - args.pval_thresh, df=1)
    print(f"  P-value threshold: {args.pval_thresh}  =>  Z² threshold: {z2_thresh:.3f}")

    snp_idxs = np.arange(gwas.n_snps)

    print("  Loading Z-scores ...")
    z = gwas.read_z_scores(snp_idxs)

    sig_mask = np.any(z ** 2 > z2_thresh, axis=1)
    sig_idxs = np.where(sig_mask)[0]

    if len(sig_idxs) == 0:
        print("  No SNPs exceed the significance threshold. Exiting.")
        gwas.close()
        return

    print(f"  Found {len(sig_idxs)} significant SNPs")

    z_sig = z[sig_idxs, :]

    print("  Loading BETA values ...")
    beta_sig = gwas.read_betas(sig_idxs)

    # Build flat long-format output: one row per (SNP × phenotype)
    snp_info = gwas.snpinfo.reset_index(drop=True).iloc[sig_idxs].copy()
    n_pheno = gwas.n_gwas

    snp_info_rep = pd.DataFrame(
        np.repeat(snp_info.values, n_pheno, axis=0),
        columns=snp_info.columns
    )
    snp_info_rep["INDEX"] = np.tile(np.arange(n_pheno), len(sig_idxs))
    snp_info_rep["Z"] = z_sig.flatten()
    snp_info_rep["BETA"] = beta_sig.flatten()
    snp_info_rep["Z2"] = (z_sig ** 2).flatten()
    snp_info_rep["P"] = 2 * norm.sf(np.abs(snp_info_rep["Z"]))

    os.makedirs(os.path.dirname(os.path.abspath(args.output)), exist_ok=True)
    snp_info_rep.to_csv(args.output, sep="\t", index=False)
    print(f"  Wrote {len(snp_info_rep)} rows to: {args.output}")

    # Write unique SNP list for PLINK --extract
    snplist_path = args.output.replace(".txt", "_snplist.txt")
    unique_snps = gwas.snpinfo.reset_index(drop=True).iloc[sig_idxs]["SNP"].drop_duplicates()
    unique_snps.to_csv(snplist_path, index=False, header=False)
    print(f"  Wrote {len(unique_snps)} SNPs to: {snplist_path}")

    gwas.close()
    print("Done.")


if __name__ == "__main__":
    main()
