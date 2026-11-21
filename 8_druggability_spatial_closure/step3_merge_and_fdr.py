#!/usr/bin/env python3
"""
step3_merge_and_fdr.py

Merge the three per-modality outputs of step2_spatial_closure.py and add
Benjamini-Hochberg q-values across all testable locus x disease-map pairs
(footprint and disease map both non-empty, rotation P available).

Usage:
    python step3_merge_and_fdr.py --dir /path/to/output/closure
"""
import argparse
import os

import numpy as np
import pandas as pd

ap = argparse.ArgumentParser()
ap.add_argument("--dir", default=".", help="directory with closure_<modality>.csv")
args = ap.parse_args()

d = pd.concat([pd.read_csv(os.path.join(args.dir, f"closure_{m}.csv")) for m in ["FA", "dti_MD", "dti_L1"]],
              ignore_index=True)
d.to_csv(os.path.join(args.dir, "closure_all.csv"), index=False)
t = d[(d.n_fp_atlas > 0) & (d.n_gc > 0) & d.p_rot.notna()].copy()
p = t.p_rot.values; n = len(p); o = np.argsort(p)
q = np.empty(n); q[o] = np.minimum.accumulate((p[o] * n / np.arange(1, n + 1))[::-1])[::-1]
t["q_all"] = np.minimum(q, 1)
t.to_csv(os.path.join(args.dir, "closure_tested_pairs_with_q.csv"), index=False)
s = t[t.q_all < 0.05]
print(f"{n} testable pairs; {len(s)} at q < 0.05 in {s.snp.nunique()} loci")
print(s.sort_values("fold", ascending=False)[["modality", "snp", "disease", "fold", "p_rot", "q_all", "overlap_tracts"]].to_string(index=False))
