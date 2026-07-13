#!/usr/bin/env python3
"""
step3_merge_sumstats.py

Merge a single locus-specific voxel GWAS sumstats file with external trait
sumstats for colocalization.  For each trait, outputs a harmonized TSV with
matched alleles and flipped BETA for strand swaps.  Palindromic SNPs (A/T,
C/G) are removed.

The output files are read by step4_coloc_analysis.R.

Run example:
  python step3_merge_sumstats.py \
      --vox-sumstats /path/to/coloc/voxelgwas/UKB_dti_FA_LOC1_CLST363_VOX253.txt \
      --trait-dir    /path/to/coloc/trait_sumstats \
      --out-dir      /path/to/coloc/coloc_sumstats
"""

import os
import re
import glob
import argparse
import pandas as pd
from tqdm import tqdm


# ── Helpers ──────────────────────────────────────────────────────────────────

def is_palindromic(a1: str, a2: str) -> bool:
    return {a1, a2} in ({"A", "T"}, {"C", "G"})


def read_sumstats(fp: str) -> pd.DataFrame:
    return pd.read_csv(fp, sep="\t")


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Merge locus voxel GWAS sumstats with trait sumstats for coloc."
    )
    ap.add_argument("--vox-sumstats", required=True,
                    help="Single voxel GWAS sumstats file (output of step 2). "
                         "Filename must match: UKB_<modality>_LOC<N>_CLST<N>_VOX<N>[.txt]")
    ap.add_argument("--trait-dir", required=True,
                    help="Directory of trait sumstats files (*_coloc.txt). "
                         "Each file must have columns: SNP, EA, NEA, beta, se")
    ap.add_argument("--out-dir", required=True,
                    help="Output directory for per-trait merged files")
    ap.add_argument("--window", type=int, default=250_000,
                    help="Half-window in bp to restrict locus (default: 250000)")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    # Discover trait sumstats
    trait_files = {
        os.path.basename(fp).replace("_coloc.txt", ""): fp
        for fp in glob.glob(os.path.join(args.trait_dir, "*_coloc.txt"))
    }
    if not trait_files:
        raise FileNotFoundError(f"No *_coloc.txt files found in: {args.trait_dir}")

    # Parse locus/cluster/voxel from filename
    fname = os.path.basename(args.vox_sumstats)
    m = re.search(r"LOC(\d+)_CLST(\d+)_VOX(\d+)", fname)
    if not m:
        raise ValueError(
            f"Filename does not match expected pattern LOC<N>_CLST<N>_VOX<N>: {fname}"
        )
    locus_id, clst, vox = int(m.group(1)), int(m.group(2)), int(m.group(3))

    # Load and filter imaging sumstats to ±window of locus centre
    img = read_sumstats(args.vox_sumstats)
    img = img.drop(columns=[c for c in ("Z", "P") if c in img.columns], errors="ignore")

    for c in ("CHR", "POS", "A1", "A2", "BETA", "SE"):
        if c not in img.columns:
            raise ValueError(f"Voxel sumstats missing required column '{c}'")

    img["_CHR"] = img["CHR"].astype(str).str.replace("chr", "", regex=False).str.strip()
    img["_POS"] = pd.to_numeric(img["POS"], errors="coerce")

    chr_center = img["_CHR"].mode().iloc[0]
    pos_center = int(img["_POS"].dropna().median())
    lo, hi = pos_center - args.window, pos_center + args.window

    img = img.loc[
        (img["_CHR"] == chr_center) &
        (img["_POS"] >= lo) &
        (img["_POS"] <= hi)
    ].drop(columns=["_CHR", "_POS"]).copy()

    img["A1"] = img["A1"].astype(str).str.strip().str.upper()
    img["A2"] = img["A2"].astype(str).str.strip().str.upper()
    img = img.drop_duplicates("SNP")

    # Merge with each trait
    for trait, trait_fp in tqdm(trait_files.items(), desc="Traits"):
        out_fp = os.path.join(
            args.out_dir,
            f"LOC{locus_id}_CLST{clst}_VOX{vox}_chr{chr_center}_pos{pos_center}_{trait}.txt"
        )
        if os.path.exists(out_fp):
            continue

        try:
            trait_df = read_sumstats(trait_fp)
        except Exception as e:
            print(f"[WARN] Could not read {trait_fp}: {e}")
            continue

        # Standardize allele columns (EA/NEA -> A1/A2)
        trait_df = trait_df.rename(columns={"EA": "A1", "NEA": "A2", "beta": "BETA", "se": "SE"})
        for c in ("SNP", "A1", "A2", "BETA", "SE"):
            if c not in trait_df.columns:
                print(f"[WARN] Trait {trait} missing column '{c}'. Skipping.")
                break
        else:
            trait_df["A1"] = trait_df["A1"].astype(str).str.strip().str.upper()
            trait_df["A2"] = trait_df["A2"].astype(str).str.strip().str.upper()
            trait_df = trait_df.drop_duplicates("SNP")

            merged = trait_df.merge(img, on="SNP", how="inner", suffixes=("_trait", "_img"))
            if merged.empty:
                continue

            # Remove palindromic SNPs
            pal = merged.apply(
                lambda r: is_palindromic(r["A1_trait"], r["A2_trait"]) or
                          is_palindromic(r["A1_img"],   r["A2_img"]),
                axis=1
            )
            merged = merged.loc[~pal]
            if merged.empty:
                continue

            # Keep only allele-concordant rows (matched or flipped)
            concordant = (
                ((merged["A1_trait"] == merged["A1_img"]) & (merged["A2_trait"] == merged["A2_img"])) |
                ((merged["A1_trait"] == merged["A2_img"]) & (merged["A2_trait"] == merged["A1_img"]))
            )
            merged = merged.loc[concordant].copy()
            if merged.empty:
                continue

            # Flip imaging BETA when alleles are swapped relative to trait
            merged["BETA_img"] = pd.to_numeric(merged["BETA_img"], errors="coerce")
            swap = merged["A1_trait"] == merged["A2_img"]
            merged.loc[swap, "BETA_img"] = -merged.loc[swap, "BETA_img"]

            # Tidy allele columns
            merged = merged.rename(columns={"A1_trait": "A1", "A2_trait": "A2"})
            merged = merged.drop(columns=["A1_img", "A2_img"], errors="ignore")

            merged.to_csv(out_fp, sep="\t", index=False)

    print(f"Merged sumstats written to: {args.out_dir}")


if __name__ == "__main__":
    main()
