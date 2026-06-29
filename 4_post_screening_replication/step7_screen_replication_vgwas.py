#!/usr/bin/env python3
"""
step7_screen_replication_vgwas.py

Screen validation voxel GWAS hits using the same cluster-size criterion
applied in the discovery (step 3), but using the validation null distribution
(from 5_cluster_analysis_val.slurm).

A SNP passes screening if the number of voxels with p < sig_thresh in the
validation exceeds the (1 - alpha)^eff_num quantile of the validation null.

Run example:
  python step7_screen_replication_vgwas.py \
      --modality dti_FA --cluster 42 --null-tag 1e4
"""

import os
import argparse
import numpy as np
import pandas as pd


MAIN_DIR = "/path/to/output"    # Override with --main-dir


def load_null_counts(cluster_id: int, tag: str, cluster_val_dir: str,
                     batches: range = range(1, 6)) -> np.ndarray:
    chunks, missing = [], []
    for b in batches:
        fp = os.path.join(cluster_val_dir, f"UKB_{{}}_phase6_{cluster_id}_{b}.txt")
        # Filled in main() once modality is known
        if os.path.exists(fp) and os.path.getsize(fp) > 0:
            chunks.append(np.atleast_1d(np.loadtxt(fp)))
        else:
            missing.append(b)

    if not chunks:
        raise FileNotFoundError(
            f"No validation null files for cluster {cluster_id} (tag={tag}, "
            f"expected batches {list(batches)})"
        )
    if missing:
        print(f"  [WARN] Missing null batches for cluster {cluster_id}: {missing}")
    return np.concatenate(chunks)


def tag_to_threshold(tag: str) -> float:
    if not tag.startswith("1e"):
        raise ValueError(f"Unsupported tag format: {tag}")
    exponent = tag[2:]
    if not exponent.isdigit():
        raise ValueError(f"Invalid exponent in tag: {tag}")
    return 10 ** (-int(exponent))


def main():
    ap = argparse.ArgumentParser(
        description="Screen replication voxel GWAS using cluster-size null."
    )
    ap.add_argument("--modality",          required=True)
    ap.add_argument("--cluster",           type=int, required=True,
                    help="Fiber cluster ID (1–430)")
    ap.add_argument("--null-tag",          default="1e4",
                    help="Null subfolder/file tag (e.g. 1e4, 1e3). "
                         "Also used as the VGWAS sig-thresh tag. (default: 1e4)")
    ap.add_argument("--main-dir",          default=MAIN_DIR)
    ap.add_argument("--cluster-val-dir",   default=None,
                    help="Directory with validation null txt files. "
                         "Default: <main-dir>/output_<modality>/cluster_val/<tag>/")
    args = ap.parse_args()

    modality   = args.modality
    cluster_id = args.cluster
    tag        = args.null_tag
    thresh     = tag_to_threshold(tag)

    modality_dir     = os.path.join(args.main_dir, f"output_{modality}")
    rep_dir          = os.path.join(modality_dir, "replication", "tract")
    cluster_val_dir  = args.cluster_val_dir or os.path.join(
        modality_dir, "cluster_val", tag
    )

    print(f"Cluster {cluster_id} | tag={tag} | threshold={thresh}")

    # Input: validation VGWAS hits at threshold=tag
    sig_file = os.path.join(rep_dir, f"UKB_{modality}_phase6_{cluster_id}_vgwas.txt")
    if not os.path.exists(sig_file):
        print(f"  No validation VGWAS file found: {sig_file}. Skipping.")
        return

    sig_data = pd.read_csv(sig_file, sep="\t")
    if "SNP" not in sig_data.columns:
        raise ValueError(f"Column 'SNP' missing in {sig_file}")

    data_count = sig_data["SNP"].value_counts()
    print(f"{len(data_count)} unique SNPs in validation VGWAS file.")

    # Load validation null distribution
    null_file_pattern = os.path.join(
        cluster_val_dir, f"UKB_{modality}_phase6_{cluster_id}_{{b}}.txt"
    )
    chunks, missing = [], []
    for b in range(1, 6):
        fp = null_file_pattern.format(b=b)
        if os.path.exists(fp) and os.path.getsize(fp) > 0:
            chunks.append(np.atleast_1d(np.loadtxt(fp)))
        else:
            missing.append(b)

    if not chunks:
        raise FileNotFoundError(
            f"No validation null files found for cluster {cluster_id} in {cluster_val_dir}"
        )
    if missing:
        print(f"  [WARN] Missing null batches: {missing}")

    null_count = np.concatenate(chunks)
    print(f"Null distribution size: {null_count.shape[0]}")

    # Effective number of voxels
    eff_num_file = os.path.join(modality_dir, "eff_num_voxels.csv")
    eff_num_df   = pd.read_csv(eff_num_file)
    eff_num_roi  = eff_num_df["eff_num"].iloc[cluster_id - 1]

    quantile = np.quantile(null_count, (1 - thresh) ** eff_num_roi)
    print(f"Cluster-size threshold: {quantile:.2f}  "
          f"(null values > threshold: {np.sum(null_count > quantile)})")

    count_screen    = data_count[data_count > quantile]
    sig_data_screen = sig_data[sig_data["SNP"].isin(count_screen.index)]
    print(f"{len(count_screen)} SNPs pass | {sig_data_screen.shape[0]} associations kept.")

    out_path = os.path.join(
        rep_dir, f"UKB_{modality}_phase6_{cluster_id}_sig_vgwas_screen_{tag}.txt"
    )
    sig_data_screen.to_csv(out_path, sep="\t", index=False)
    print(f"Wrote: {out_path}")


if __name__ == "__main__":
    main()
