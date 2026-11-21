#!/usr/bin/env python3
"""
Step 2: locus-level replication for ONE tract.

A discovery locus of the tract is replicated when an ABCD association (location-level P < SIG_THRESH,
discovery candidate SNPs of that tract) lies within +-WINDOW bp of the discovery lead variant, in any
fiber cluster of the tract. As in the original ABCD analysis, no cluster-size screen is applied to ABCD.

Logic is unchanged from the original script; only the paths are taken from rep_config.sh.
Output: REP_DIR/<DISC_PREFIX>_Tract_<TRACT>_clump_replication_by_ABCD_<TAG>.txt
"""
import argparse
import os

import pandas as pd

E = os.environ


def find_loc(num_list, target):
    l, r = 0, len(num_list) - 1
    while l <= r:
        mid = (l + r) // 2
        if num_list[mid] == target:
            return mid
        elif num_list[mid] > target:
            r = mid - 1
        else:
            l = mid + 1
    return r


def merge_nearby_loci(clump: pd.DataFrame) -> pd.DataFrame:
    """Within a tract, keep the most significant lead among leads closer than 500 kb."""
    clump = clump.sort_values(["CHR", "BP"])
    if clump.shape[0] == 0:
        return clump
    keep_idx = [0]
    for i in range(1, clump.shape[0]):
        chr_ = clump.iloc[i]["CHR"]
        pos = clump.iloc[i]["BP"]
        p = clump.iloc[i]["P"]
        while keep_idx:
            prev = keep_idx[-1]
            if (clump.iloc[prev]["CHR"] == chr_) and (pos - clump.iloc[prev]["BP"] < 500000):
                if p < clump.iloc[prev]["P"]:
                    keep_idx.pop()
                else:
                    break
            else:
                keep_idx.append(i)
                break
        if not keep_idx:
            keep_idx.append(i)
    return clump.iloc[keep_idx]


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tract", required=True)
    ap.add_argument("--window", type=int, default=int(E.get("WINDOW", 250000)))
    args = ap.parse_args()
    tract, window = args.tract, args.window
    mod, tag, prefix = E.get("MOD", "FA"), E.get("TAG", "1e2"), E.get("DISC_PREFIX", "UKB_FA_phase16")
    rep_dir, disc_tract_dir = E["REP_DIR"], E["DISC_TRACT_DIR"]

    cluster2tract = pd.read_csv(E["CLUSTER_TRACT_CSV"])
    cluster_ids = cluster2tract.loc[cluster2tract["tract_assign"] == tract, "cluster"].tolist()
    if not cluster_ids:
        raise ValueError(f"No clusters found for tract '{tract}'")
    print(f"Processing {tract} with {len(cluster_ids)} clusters...")

    clump_path = os.path.join(disc_tract_dir, f"{prefix}_Tract_{tract}_clump_r2_0.1_step3.clumped")
    clump = pd.read_csv(clump_path, sep=r"\s+")[["CHR", "F", "SNP", "BP", "P"]].copy()
    clump["CHR"] = clump["CHR"].astype(int)
    clump["BP"] = clump["BP"].astype(int)
    print(f"{len(clump)} loci before merge.")
    clump = merge_nearby_loci(clump)
    print(f"{len(clump)} loci after merge.")

    rep_all, missing, not_run = [], 0, []
    for cid in cluster_ids:
        base = os.path.join(rep_dir, f"ABCD_{mod}_{cid}_vgwas_{tag}")
        if not os.path.exists(base + ".txt"):
            missing += 1
            if not os.path.exists(base + ".done"):
                not_run.append(cid)
            continue
        rep = pd.read_csv(base + ".txt", sep="\t")[["CHR", "POS"]].copy()
        rep["CHR"] = rep["CHR"].astype(int)
        rep["POS"] = rep["POS"].astype(int)
        rep_all.append(rep)
    print(f"Replication cluster files found={len(rep_all)}/{len(cluster_ids)} (missing={missing})")
    if not_run:   # a cluster whose voxel GWAS never ran would silently lower the replication counts
        raise RuntimeError(f"Step 1 has not finished for clusters {not_run} of tract {tract}")
    if not rep_all:
        raise RuntimeError(f"No replication cluster results found for tract {tract}")

    rep_unique = pd.concat(rep_all, ignore_index=True).drop_duplicates(["CHR", "POS"]).sort_values(["CHR", "POS"])
    replication = []
    for _, data in clump.iterrows():
        start, end = int(data["BP"]) - window, int(data["BP"]) + window
        cand = rep_unique.loc[rep_unique["CHR"] == int(data["CHR"]), "POS"].to_list()
        replication.append("Yes" if len(cand) > 0 and find_loc(cand, start) < find_loc(cand, end) else "No")
    clump["replication"] = replication
    print(clump.value_counts("replication"))

    out_path = os.path.join(rep_dir, f"{prefix}_Tract_{tract}_clump_replication_by_ABCD_{tag}.txt")
    clump.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
