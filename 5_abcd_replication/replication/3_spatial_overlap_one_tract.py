#!/usr/bin/env python3
"""
Step 3: spatial overlap for ONE tract.

For each discovery locus of the tract:
  discovery footprint = locations (cluster, index) of the tract with a screened discovery association
                        within +-WINDOW bp of the lead variant;
  ABCD footprint      = locations with an ABCD association (P < SIG_THRESH) in the same window;
  overlap_percent     = |discovery ∩ ABCD| / |discovery|   (a fraction, 0-1).

Logic is unchanged from the original script; only the paths are taken from rep_config.sh.
Output: REP_DIR/<DISC_PREFIX>_Tract_<TRACT>_clump_replication_by_ABCD_spatial_overlap_all_loci_<TAG>.txt
"""
import argparse
import os

import pandas as pd

E = os.environ


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--tract", required=True)
    ap.add_argument("--window", type=int, default=int(E.get("WINDOW", 250000)))
    args = ap.parse_args()
    tract, window = args.tract, args.window
    mod, tag, prefix = E.get("MOD", "FA"), E.get("TAG", "1e2"), E.get("DISC_PREFIX", "UKB_FA_phase16")
    rep_dir, disc_vox_dir = E["REP_DIR"], E["DISC_VOX_DIR"]

    cluster2tract = pd.read_csv(E["CLUSTER_TRACT_CSV"])
    cluster_ids = cluster2tract.loc[cluster2tract["tract_assign"] == tract, "cluster"].tolist()
    if not cluster_ids:
        raise ValueError(f"No clusters found for tract '{tract}'")
    print(f"Tract {tract}: {len(cluster_ids)} clusters")

    clump = pd.read_csv(os.path.join(rep_dir, f"{prefix}_Tract_{tract}_clump_replication_by_ABCD_{tag}.txt"), sep="\t")
    clump["CHR"] = clump["CHR"].astype(int)
    clump["BP"] = clump["BP"].astype(int)

    disc_all, rep_all = [], []
    for cid in cluster_ids:
        disc_path = os.path.join(disc_vox_dir, f"{prefix}_{cid}_sig_vgwas_screen.txt")
        rep_path = os.path.join(rep_dir, f"ABCD_{mod}_{cid}_vgwas_{tag}.txt")
        for path, store in ((disc_path, disc_all), (rep_path, rep_all)):
            if os.path.exists(path):
                d = pd.read_csv(path, sep="\t")[["INDEX", "CHR", "POS", "SNP"]]
                d["CLUSTER"] = cid
                d["CHR"] = d["CHR"].astype(int)
                d["POS"] = d["POS"].astype(int)
                store.append(d)
    if not disc_all:
        raise RuntimeError(f"No discovery voxel files found for tract {tract}")
    if not rep_all:
        raise RuntimeError(f"No replication voxel files found for tract {tract}")
    disc = pd.concat(disc_all, ignore_index=True)
    rep = pd.concat(rep_all, ignore_index=True)
    disc["CLUSTER_INDEX"] = disc["CLUSTER"].astype(str) + "-" + disc["INDEX"].astype(str)
    rep["CLUSTER_INDEX"] = rep["CLUSTER"].astype(str) + "-" + rep["INDEX"].astype(str)
    print(f"Loaded discovery rows: {disc.shape[0]}, replication rows: {rep.shape[0]}")

    rows = []
    for _, locus in clump.iterrows():
        chr_, pos = int(locus["CHR"]), int(locus["BP"])
        start, end = pos - window, pos + window
        disc_idx = set(disc.loc[(disc["CHR"] == chr_) & (disc["POS"] > start) & (disc["POS"] < end), "CLUSTER_INDEX"])
        rep_idx = set(rep.loc[(rep["CHR"] == chr_) & (rep["POS"] > start) & (rep["POS"] < end), "CLUSTER_INDEX"])
        ov = disc_idx & rep_idx
        rows.append({"chr": chr_, "snp": locus.get("SNP"), "pos": pos, "p": locus.get("P"),
                     "disc_voxels": len(disc_idx), "rep_voxels": len(rep_idx), "overlap_voxels": len(ov),
                     "overlap_percent": round(len(ov) / len(disc_idx), 3) if disc_idx else 0.0})
    output = pd.DataFrame(rows)
    print(f">0% overlap: {(output['overlap_percent'] > 0).sum()}")
    print(f">50% overlap: {(output['overlap_percent'] > 0.5).sum()}")
    out_path = os.path.join(rep_dir, f"{prefix}_Tract_{tract}_clump_replication_by_ABCD_spatial_overlap_all_loci_{tag}.txt")
    output.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
