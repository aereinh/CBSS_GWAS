#!/usr/bin/env python3
"""
step1_select_best_voxel.py

For a single locus, scan voxel GWAS screening files to find the voxel
(cluster + voxel index) with the minimum p-value within ±window bp of the
locus centre.  Output is a small CSV used by step 2 (targeted voxel GWAS).

Run example:
  python step1_select_best_voxel.py \
      --modality dti_FA \
      --locus-id 1 \
      --vox-dir  /path/to/post_screen/disc/fibercluster \
      --loci-info /path/to/loci/ALL_TRACTS_locus_summary.csv \
      --out-dir  /path/to/coloc/sel_voxels
"""

import os
import re
import glob
import argparse
import pandas as pd
import numpy as np
from typing import Any, Dict, Optional


# ── Constants ────────────────────────────────────────────────────────────────
DEFAULT_WINDOW    = 250_000   # ±250 kb around locus centre
DEFAULT_CHUNKSIZE = 1_000_000

# Column names in the loci-info CSV (override via CLI if yours differ)
LOCUS_ID_COL  = "locus_id"
LOCUS_CHR_COL = "chr"
LOCUS_BP_COL  = "bp"
LEAD_SNP_COL  = "snp"   # optional; used only for output naming


# ── Helpers ──────────────────────────────────────────────────────────────────

def safe_filename(x: Any) -> str:
    return re.sub(r"[^A-Za-z0-9_\-\.]+", "_", str(x))


def find_best_voxel(
    locus_chr: str,
    locus_bp: int,
    vgwas_files: list,
    cluster_pattern: re.Pattern,
    window: int = DEFAULT_WINDOW,
    chunksize: int = DEFAULT_CHUNKSIZE,
) -> Optional[Dict[str, Any]]:
    """Scan voxel GWAS files and return the single best (min P) hit in window."""
    locus_chr = str(locus_chr).replace("chr", "").strip()
    lo = locus_bp - window
    hi = locus_bp + window

    best: Optional[Dict[str, Any]] = None

    for fp in vgwas_files:
        m = cluster_pattern.search(fp)
        if not m:
            continue
        cluster_id = int(m.group(1))

        try:
            chunks = pd.read_csv(
                fp, sep="\t", dtype=str,
                usecols=lambda c: c in ("SNP", "P", "INDEX", "CHR", "POS"),
                chunksize=chunksize,
            )
        except Exception:
            chunks = pd.read_csv(
                fp, sep=None, engine="python", dtype=str,
                usecols=lambda c: c in ("SNP", "P", "INDEX", "CHR", "POS"),
                chunksize=chunksize,
            )

        file_best: Optional[Dict[str, Any]] = None

        for chunk in chunks:
            if "CHR" not in chunk.columns or "POS" not in chunk.columns:
                break

            chr_col = chunk["CHR"].astype(str).str.replace("chr", "", regex=False).str.strip()
            pos_col = pd.to_numeric(chunk["POS"], errors="coerce")
            in_win  = (chr_col == locus_chr) & (pos_col >= lo) & (pos_col <= hi)
            sub     = chunk.loc[in_win].copy()
            if sub.empty:
                continue

            sub["P"] = pd.to_numeric(sub["P"], errors="coerce")
            sub = sub.dropna(subset=["P"])
            if sub.empty:
                continue

            row = sub.loc[sub["P"].idxmin()]
            hit: Dict[str, Any] = {
                "cluster": cluster_id,
                "INDEX":   row.get("INDEX"),
                "CHR":     row.get("CHR"),
                "POS":     row.get("POS"),
                "SNP":     row.get("SNP"),
                "P":       float(row["P"]),
            }
            if file_best is None or hit["P"] < file_best["P"]:
                file_best = hit

        if file_best is not None and (best is None or file_best["P"] < best["P"]):
            best = file_best

    return best


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Find best voxel within ±window of a locus for colocalization."
    )
    ap.add_argument("--modality",      default="dti_FA",
                    help="Imaging modality (e.g. dti_FA, dti_MD)")
    ap.add_argument("--locus-id",      type=int, required=True,
                    help="Locus ID to process")
    ap.add_argument("--vox-dir",       required=True,
                    help="Directory containing voxel GWAS screen files "
                         "(e.g. post_screen/disc/fibercluster/)")
    ap.add_argument("--loci-info",     required=True,
                    help="CSV with columns: locus_id, chr, bp[, snp]")
    ap.add_argument("--out-dir",       required=True,
                    help="Output directory for per-locus CSV")
    ap.add_argument("--window",        type=int, default=DEFAULT_WINDOW,
                    help=f"Half-window in bp (default: {DEFAULT_WINDOW})")
    ap.add_argument("--chunksize",     type=int, default=DEFAULT_CHUNKSIZE)
    ap.add_argument("--locus-id-col",  default=LOCUS_ID_COL)
    ap.add_argument("--locus-chr-col", default=LOCUS_CHR_COL)
    ap.add_argument("--locus-bp-col",  default=LOCUS_BP_COL)
    ap.add_argument("--lead-snp-col",  default=LEAD_SNP_COL,
                    help="Column for lead SNP rsID (optional; for output naming)")
    args = ap.parse_args()

    os.makedirs(args.out_dir, exist_ok=True)

    loci = pd.read_csv(args.loci_info)
    for col in [args.locus_id_col, args.locus_chr_col, args.locus_bp_col]:
        if col not in loci.columns:
            raise ValueError(f"loci-info missing column '{col}'. Found: {loci.columns.tolist()}")

    rows = loci.loc[loci[args.locus_id_col] == args.locus_id]
    if rows.empty:
        raise ValueError(f"locus_id={args.locus_id} not found in {args.loci_info}")
    row = rows.iloc[0]

    locus_chr = row[args.locus_chr_col]
    locus_bp  = int(row[args.locus_bp_col])
    lead_snp  = row.get(args.lead_snp_col) if args.lead_snp_col in loci.columns else None

    # Discover voxel GWAS files
    pattern = os.path.join(args.vox_dir, f"UKB_{args.modality}_*_sig_vgwas_screen.txt")
    rex     = re.compile(rf"UKB_{re.escape(args.modality)}_(\d+)_sig_vgwas_screen\.txt$")
    files   = sorted(glob.glob(pattern))

    if not files:
        raise FileNotFoundError(f"No voxel GWAS files matched: {pattern}")

    best = find_best_voxel(
        locus_chr=str(locus_chr),
        locus_bp=locus_bp,
        vgwas_files=files,
        cluster_pattern=rex,
        window=args.window,
        chunksize=args.chunksize,
    )

    locus_tok = safe_filename(args.locus_id)
    snp_tok   = safe_filename(lead_snp) if (lead_snp is not None and not pd.isna(lead_snp)) else "NA"
    out_fp    = os.path.join(args.out_dir, f"locus_{locus_tok}_lead_{snp_tok}_best_voxel.csv")

    if best is None:
        record = {
            "locus_id": args.locus_id, "lead_snp": lead_snp,
            "locus_chr": locus_chr, "locus_bp": locus_bp, "window_bp": args.window,
            "cluster": None, "INDEX": None, "CHR": None, "POS": None,
            "SNP": None, "P": None, "status": "NOT_FOUND",
        }
    else:
        record = {
            "locus_id": args.locus_id, "lead_snp": lead_snp,
            "locus_chr": locus_chr, "locus_bp": locus_bp, "window_bp": args.window,
            **best, "status": "OK",
        }

    pd.DataFrame([record]).to_csv(out_fp, index=False)
    print(f"Written: {out_fp}")
    print(pd.DataFrame([record]).to_string(index=False))


if __name__ == "__main__":
    main()
