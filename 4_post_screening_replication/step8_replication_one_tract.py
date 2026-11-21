#!/usr/bin/env python3
"""
step8_replication_one_tract.py

Assess genetic replication for a single white-matter tract by checking whether
validation (Phase 7) screened hits fall within ±window bp of each discovery
LD-clumped locus.

For each discovery locus, outputs "Yes" (replication) or "No".

Run example:
  python step8_replication_one_tract.py --modality dti_FA --tract SCC
"""

import os
import argparse
import pandas as pd
from collections import defaultdict


def find_loc(positions: list, target: int) -> int:
    """Binary search: return index of largest element <= target."""
    l, r = 0, len(positions) - 1
    while l <= r:
        mid = (l + r) // 2
        if positions[mid] == target:
            return mid
        elif positions[mid] > target:
            r = mid - 1
        else:
            l = mid + 1
    return r


def merge_nearby_loci(clump: pd.DataFrame, min_distance: int = 500_000) -> pd.DataFrame:
    """Remove duplicate loci within min_distance bp, keeping the most significant."""
    clump = clump.sort_values(["CHR", "BP"]).reset_index(drop=True)
    if clump.empty:
        return clump

    keep = [0]
    for i in range(1, len(clump)):
        prev = keep[-1]
        if (clump.loc[i, "CHR"] == clump.loc[prev, "CHR"] and
                clump.loc[i, "BP"] - clump.loc[prev, "BP"] < min_distance):
            # Keep whichever has lower p-value
            if clump.loc[i, "P"] < clump.loc[prev, "P"]:
                keep[-1] = i
        else:
            keep.append(i)
    return clump.iloc[keep].reset_index(drop=True)


def main():
    ap = argparse.ArgumentParser(
        description="Genetic replication check for a single white-matter tract."
    )
    ap.add_argument("--modality",  required=True)
    ap.add_argument("--tract",     required=True,
                    help="Tract name (as in cluster_tract_summary.csv), e.g. SCC")
    ap.add_argument("--window",    type=int, default=250_000,
                    help="Half-window for replication (bp, default: 250000)")
    ap.add_argument("--null-tag",  default="1e4",
                    help="Tag used in step7 screened filenames (default: 1e4)")
    ap.add_argument("--main-dir",  default="/path/to/output")
    ap.add_argument("--cluster-to-tract", default=None,
                    help="Path to cluster_tract_summary.csv. "
                         "Default: <main-dir>/output_<modality>/cluster_tract_summary.csv")
    args = ap.parse_args()

    modality = args.modality
    tract    = args.tract
    window   = args.window
    tag      = args.null_tag

    modality_dir = os.path.join(args.main_dir, f"output_{modality}")
    disc_dir     = os.path.join(modality_dir, "post_screen", "disc", "tract")
    rep_dir      = os.path.join(modality_dir, "replication", "tract")

    # Load cluster-to-tract assignments
    c2t_file = args.cluster_to_tract or os.path.join(modality_dir, "cluster_tract_summary.csv")
    c2t      = pd.read_csv(c2t_file)
    tract_clusters = c2t.loc[c2t["tract_assign"] == tract, "cluster"].tolist()
    if not tract_clusters:
        raise ValueError(f"No clusters found for tract '{tract}' in {c2t_file}")
    print(f"Tract {tract}: {len(tract_clusters)} cluster(s)")

    # Discovery clumped loci
    clump_file = os.path.join(
        disc_dir,
        f"UKB_{modality}_phase16_Tract_{tract}_clump_r2_0.1_step3.clumped"
    )
    if not os.path.exists(clump_file):
        raise FileNotFoundError(f"Clump file not found: {clump_file}")

    clump = pd.read_csv(clump_file, sep=r"\s+")[["CHR", "F", "SNP", "BP", "P"]].copy()
    clump["CHR"] = clump["CHR"].astype(int)
    clump["BP"]  = clump["BP"].astype(int)

    print(f"{len(clump)} loci before deduplication.")
    clump = merge_nearby_loci(clump)
    print(f"{len(clump)} loci after deduplication.")

    # Aggregate replication screened hits across all clusters in the tract
    rep_hits = []
    n_missing = 0
    for cid in tract_clusters:
        rep_file = os.path.join(rep_dir, f"UKB_{modality}_phase7_{cid}_sig_vgwas_screen_{tag}.txt")
        if not os.path.exists(rep_file):
            n_missing += 1
            continue
        rep = pd.read_csv(rep_file, sep="\t")[["CHR", "POS"]].copy()
        rep["CHR"] = rep["CHR"].astype(int)
        rep["POS"] = rep["POS"].astype(int)
        rep_hits.append(rep)

    print(f"Replication files: found {len(rep_hits)}/{len(tract_clusters)} (missing={n_missing})")
    if not rep_hits:
        raise RuntimeError(f"No replication results found for tract {tract}")

    rep_sig = (
        pd.concat(rep_hits, ignore_index=True)
        .drop_duplicates(["CHR", "POS"])
        .sort_values(["CHR", "POS"])
    )

    # Check each discovery locus for replication
    replication = []
    for _, locus in clump.iterrows():
        start = int(locus["BP"]) - window
        end   = int(locus["BP"]) + window
        chr_pos = rep_sig.loc[rep_sig["CHR"] == int(locus["CHR"]), "POS"].tolist()
        # Replicated if any validation hit falls within the window
        lo = find_loc(chr_pos, start)
        hi = find_loc(chr_pos, end)
        replication.append("Yes" if lo < hi else "No")

    clump["replication"] = replication
    print(clump["replication"].value_counts().to_string())

    out_file = os.path.join(
        disc_dir,
        f"UKB_{modality}_phase16_Tract_{tract}_clump_replication_by_phase7_{tag}.txt"
    )
    clump.to_csv(out_file, sep="\t", index=False)
    print(f"Wrote: {out_file}")


if __name__ == "__main__":
    main()
