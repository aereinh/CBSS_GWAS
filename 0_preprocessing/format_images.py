#!/usr/bin/env python3
"""
format_images.py

Convert the CBSS along-fiber measurements of one UK Biobank imaging phase into
one text file per fiber cluster, the input format of HEIG --read-image.

The CBSS pipeline stores each diffusion measure as MATLAB files
`<k>_Aligned.mat` (variable `f_aligned`, 100 points x 10 fibers x subjects),
ten atlas fibers per file, together with the subject list `UKB_<phase>.txt`.
For every fiber cluster of the atlas this script collects the fibers of the
cluster and writes

    <out-dir>/fibersClust_<cluster>.txt     FID IID X0 X1 ...  (100 columns per fiber)

Run once per phase and modality: Phases 1-6 form the discovery sample and
Phase 7 the validation sample (see 1_gwas_discovery/1_read_images.slurm).

Usage:
    python format_images.py \
        --data-dir    /path/to/CBSS/dti_FA/step2_dti_FA_UKB_Phase1 \
        --subject-ids /path/to/CBSS/dti_FA/step2_dti_FA_UKB_Phase1/UKB_Phase1.txt \
        --cluster-csv /path/to/Fiber_atlas_100points_FullName.csv \
        --out-dir     /path/to/Modalities/dti_FA/step2_dti_FA_UKB_Phase1
"""
import argparse
import os

import numpy as np
import pandas as pd
import scipy.io as sio


def parse_args():
    p = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--data-dir", required=True, help="Directory with the <k>_Aligned.mat files of one phase")
    p.add_argument("--subject-ids", required=True, help="Text file with one subject ID per line (order of the .mat files)")
    p.add_argument("--cluster-csv", required=True, help="Fiber atlas table with the cluster label of each fiber")
    p.add_argument("--cluster-col", default="Clust_Dis30mm", help="Column of --cluster-csv with the cluster label")
    p.add_argument("--out-dir", required=True, help="Output directory for fibersClust_<cluster>.txt")
    p.add_argument("--fibers-per-file", type=int, default=10, help="Atlas fibers stored in each .mat file")
    p.add_argument("--points-per-fiber", type=int, default=100, help="Points along each fiber")
    p.add_argument("--clusters", default=None, help="Optional subset of clusters, e.g. 1-50 or 7,19,200")
    return p.parse_args()


def parse_clusters(spec, n_clust):
    if spec is None:
        return list(range(1, n_clust + 1))
    out = []
    for part in spec.split(","):
        if "-" in part:
            a, b = part.split("-")
            out.extend(range(int(a), int(b) + 1))
        else:
            out.append(int(part))
    return out


def main():
    args = parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    subj_ids = np.loadtxt(args.subject_ids, dtype=int)
    n = len(subj_ids)

    clust_labs = pd.read_csv(args.cluster_csv)[args.cluster_col]
    n_clust = len(np.unique(clust_labs))
    npf, npt = args.fibers_per_file, args.points_per_fiber

    for clust in parse_clusters(args.clusters, n_clust):
        # 0-based indices of the fibers of this cluster, grouped by the .mat file that holds them
        indices = clust_labs[clust_labs == clust].index
        file_groups = {}
        for fiber_index in indices:
            file_groups.setdefault(fiber_index // npf + 1, []).append(fiber_index)

        blocks = [np.column_stack((subj_ids, subj_ids))]
        for file_index, fibers in file_groups.items():
            mat = sio.loadmat(os.path.join(args.data_dir, f"{file_index}_Aligned.mat"))
            # f_aligned: points x fibers x subjects -> subjects x (fibers * points)
            pts = mat["f_aligned"].transpose(2, 0, 1).reshape(npf * npt, n).T
            cols = np.concatenate([np.arange(f * npt % (npt * npf), f * npt % (npt * npf) + npt) for f in fibers])
            blocks.append(pts[:, cols])

        df = pd.DataFrame(np.hstack(blocks))
        df.columns = ["FID", "IID"] + [f"X{k}" for k in range(df.shape[1] - 2)]
        df["FID"] = df["FID"].astype(int)
        df["IID"] = df["IID"].astype(int)
        df.to_csv(os.path.join(args.out_dir, f"fibersClust_{clust}.txt"), sep=" ", index=False)
        print(f"cluster {clust}: {len(indices)} fibers, {n} subjects")


if __name__ == "__main__":
    main()
