#!/usr/bin/env python3
"""
step5_process_results.py

Merge the coloc [RESULTS] from step 4 with the replicated loci table,
keeping only loci with spatial overlap >= overlap_thresh.  The locus ID is
parsed from the "file" field in the coloc results (pattern: LOC<N>_...).

Run example:
  python step5_process_results.py \
      --loci       /path/to/loci/ALL_TRACTS_locus_summary.csv \
      --coloc      /path/to/coloc/coloc_results.txt \
      --output     /path/to/coloc/coloc_results_replicated.tsv \
      --overlap-thresh 0.5 \
      --min-pph4   0.5
"""

import re
import csv
import argparse
from pathlib import Path
from typing import Dict, Iterable, List, Optional, Tuple


# ── Parsing helpers ───────────────────────────────────────────────────────────

LOCUS_ID_RE = re.compile(r"LOC(?P<locus_id>\d+)", re.IGNORECASE)


def _float(x) -> Optional[float]:
    try:
        s = str(x).strip()
        return float(s) if s and s.lower() != "na" else None
    except Exception:
        return None


def _int(x) -> Optional[int]:
    try:
        s = str(x).strip()
        return int(float(s)) if s and s.lower() != "na" else None
    except Exception:
        return None


def iter_results_block(coloc_path: Path) -> Iterable[Dict[str, str]]:
    """Yield rows from the [RESULTS] block of a coloc results file."""
    in_results = False
    header: Optional[List[str]] = None

    with coloc_path.open("r", newline="") as fh:
        for raw in fh:
            line = raw.rstrip("\n")
            if not in_results:
                if line.strip() == "[RESULTS]":
                    in_results = True
                continue
            if line.startswith("[") and line.strip() != "[RESULTS]":
                break
            if not line or line.startswith("#"):
                continue
            if header is None:
                header = line.split("\t")
                continue
            values = line.split("\t")
            if len(values) == len(header):
                yield dict(zip(header, values))


def load_replicated_loci(
    loci_csv: Path, overlap_col: str, overlap_thresh: float
) -> Tuple[Dict[int, Dict[str, str]], List[str]]:
    loci: Dict[int, Dict[str, str]] = {}
    cols: List[str] = []

    with loci_csv.open("r", newline="") as fh:
        reader = csv.DictReader(fh)
        cols = [c for c in (reader.fieldnames or []) if c]

        if "locus_id" not in cols:
            raise ValueError(f"loci CSV missing 'locus_id' column. Found: {cols}")
        if overlap_col not in cols:
            raise ValueError(f"loci CSV missing '{overlap_col}' column. Found: {cols}")

        for row in reader:
            lid = _int(row.get("locus_id"))
            if lid is None:
                continue
            if (_float(row.get(overlap_col)) or 0.0) < overlap_thresh:
                continue
            loci[lid] = row

    return loci, cols


# ── Main ─────────────────────────────────────────────────────────────────────

def main():
    ap = argparse.ArgumentParser(
        description="Merge coloc [RESULTS] with replicated loci by locus_id."
    )
    ap.add_argument("--loci",          required=True,
                    help="Replicated loci CSV (must have locus_id, overlap column)")
    ap.add_argument("--coloc",         required=True,
                    help="coloc_results.txt from step4_coloc_analysis.R")
    ap.add_argument("--output",        required=True,
                    help="Output TSV path")
    ap.add_argument("--overlap-col",   default="overlap_percent",
                    help="Column in loci CSV for spatial overlap (default: overlap_percent)")
    ap.add_argument("--overlap-thresh", type=float, default=0.5,
                    help="Minimum overlap to consider a locus replicated (default: 0.5)")
    ap.add_argument("--min-pph4",      type=float, default=0.0,
                    help="Minimum PP.H4.abf to include (default: 0, keep all)")
    args = ap.parse_args()

    loci, loci_cols = load_replicated_loci(
        Path(args.loci), args.overlap_col, args.overlap_thresh
    )

    merged_rows: List[Dict[str, str]] = []
    coloc_cols:  Optional[List[str]]  = None
    matched = unmatched = parse_fail = filtered = 0

    for row in iter_results_block(Path(args.coloc)):
        if coloc_cols is None:
            coloc_cols = list(row.keys())

        m = LOCUS_ID_RE.search(row.get("file", ""))
        if not m:
            parse_fail += 1
            continue

        lid = _int(m.group("locus_id"))
        if lid not in loci:
            unmatched += 1
            continue

        if (_float(row.get("PP.H4.abf")) or 0.0) < args.min_pph4:
            filtered += 1
            continue

        merged = {**loci[lid], **row, "locus_id": str(lid)}
        matched += 1
        merged_rows.append(merged)

    if coloc_cols is None:
        raise ValueError("No coloc results found in [RESULTS] block.")

    # Build unified fieldname list
    seen = set()
    fieldnames: List[str] = []
    for nm in loci_cols + coloc_cols:
        if nm and nm not in seen:
            fieldnames.append(nm)
            seen.add(nm)

    Path(args.output).parent.mkdir(parents=True, exist_ok=True)
    with open(args.output, "w", newline="") as fh:
        writer = csv.DictWriter(fh, fieldnames=fieldnames, delimiter="\t", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(merged_rows)

    print(f"Wrote {len(merged_rows)} rows to: {args.output}")
    print(f"Replicated loci (overlap >= {args.overlap_thresh}): {len(loci)}")
    print(f"Matched:         {matched}")
    print(f"Unmatched:       {unmatched}")
    print(f"Parse failures:  {parse_fail}")
    print(f"Filtered (PP.H4 < {args.min_pph4}): {filtered}")


if __name__ == "__main__":
    main()
