#!/usr/bin/env python3
"""Extract one fiber cluster (100-point aligned values per fiber) for all ABCD subjects
into a PLINK-style text file:  <out-dir>/fibersClust_<cluster>.txt  (FID IID X0 X1 ...).

Same logic as the original 0_format_images-1cluster.py; works for any modality whose
<data-dir> holds <k>_Aligned.mat files (key 'f_aligned') and ABCD.txt (subject IDs).
IDs get the 'NDAR_' prefix so they match the genotype/covariate files.
"""
import argparse
import os

import numpy as np
import pandas as pd
import scipy.io as sio


def run_one_cluster(cluster_index, data_dir, out_dir, subj_id_file, clust_file,
                    col="Clust_Dis30mm", nfibers_per_file=10, nfiber_pts=100,
                    file_pattern="{}_Aligned.mat", f_aligned_key="f_aligned"):
    if cluster_index < 1:
        raise ValueError("cluster_index must be 1-based (>= 1).")
    os.makedirs(out_dir, exist_ok=True)

    subj_ids = np.loadtxt(subj_id_file, dtype=str)
    n = len(subj_ids)
    if n == 0:
        raise ValueError(f"No subject IDs found in {subj_id_file}")

    clust_data = pd.read_csv(clust_file)
    if col not in clust_data.columns:
        raise ValueError(f"Column '{col}' not in {clust_file}")
    clust_labs = clust_data[col]

    # fibers in this cluster (0-based fiber indices), grouped by .mat file
    indices_i = clust_labs[clust_labs == cluster_index].index.to_numpy()
    if len(indices_i) == 0:
        raise ValueError(f"No fibers for cluster {cluster_index} in column '{col}'")
    file_grps = {}
    for fiber_index in indices_i:
        file_grps.setdefault(int(fiber_index) // nfibers_per_file + 1, []).append(int(fiber_index))

    blocks = [np.column_stack((subj_ids, subj_ids))]
    for file_index, fiber_indices in file_grps.items():
        file_path = os.path.join(data_dir, file_pattern.format(file_index))
        if not os.path.exists(file_path):
            raise FileNotFoundError(f"Missing aligned mat file: {file_path}")
        offsets = [(fi % nfibers_per_file) * nfiber_pts for fi in fiber_indices]
        pt_inds = np.concatenate([np.arange(o, o + nfiber_pts) for o in offsets])

        print(f"[cluster {cluster_index}] {file_path}: {len(fiber_indices)} fibers")
        mat = sio.loadmat(file_path)
        if f_aligned_key not in mat:
            raise KeyError(f"'{f_aligned_key}' not in {file_path}")
        fiberpts = mat[f_aligned_key].transpose(2, 0, 1).reshape(nfibers_per_file * nfiber_pts, n).T
        blocks.append(fiberpts[:, pt_inds])

    df = pd.DataFrame(np.hstack(blocks))
    df.columns = ["FID", "IID"] + [f"X{i}" for i in range(df.shape[1] - 2)]
    for c in ("FID", "IID"):
        df[c] = df[c].astype(str).apply(lambda x: x if x.startswith("NDAR_") else "NDAR_" + x)

    out_fp = os.path.join(out_dir, f"fibersClust_{cluster_index}.txt")
    df.to_csv(out_fp, sep=" ", index=False)
    print(f"[cluster {cluster_index}] wrote {out_fp} ({df.shape[0]} subjects, {df.shape[1] - 2} values)")
    return out_fp


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--cluster", type=int, required=True, help="1-based fiber cluster index")
    ap.add_argument("--data-dir", required=True, help="folder with <k>_Aligned.mat and ABCD.txt")
    ap.add_argument("--out-dir", required=True)
    ap.add_argument("--clust-file", required=True, help="fiber atlas CSV with cluster labels")
    ap.add_argument("--subj-id-file", default=None, help="default: <data-dir>/ABCD.txt")
    ap.add_argument("--col", default="Clust_Dis30mm")
    args = ap.parse_args()

    run_one_cluster(
        cluster_index=args.cluster,
        data_dir=args.data_dir,
        out_dir=args.out_dir,
        subj_id_file=args.subj_id_file or os.path.join(args.data_dir, "ABCD.txt"),
        clust_file=args.clust_file,
        col=args.col,
    )


if __name__ == "__main__":
    main()
