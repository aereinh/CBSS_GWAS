import os
import pandas as pd
from tqdm import tqdm
from itertools import product

main_dir = '/path/to/UKB_PSC/multiomics/'
types = ['proteomics', 'metabolomics', 'biomarkers']
nroi = 430

all_res = []

# Iterate over the cartesian product of (roi, type) with a single progress bar
pairs = product(range(1, nroi + 1), types)
pairs = tqdm(list(pairs), desc="Scanning selected-feature sumstats", dynamic_ncols=True, leave=False)

for i, typ in pairs:
    input_file = os.path.join(main_dir, typ, 'voxel-assoc', 'sel_feats', f'voxel_assoc_{i}_sel_feats.txt')
    if not os.path.exists(input_file):
        continue

    dat = pd.read_csv(input_file, sep='\t')
    keep_cols = [c for c in ['VARIABLE', 'INDEX', 'Z', 'P'] if c in dat.columns]
    dat = dat[keep_cols].copy()

    dat['type'] = typ
    dat['fiber_clust'] = i
    all_res.append(dat)

if not all_res:
    print("No selected-feature sumstats found — nothing to save.")
else:
    all_sel_feat_sumstats = pd.concat(all_res, ignore_index=True)

    outdir = os.path.join(main_dir, 'results')
    os.makedirs(outdir, exist_ok=True)

    outpath = os.path.join(outdir, 'sel_features_sumstats_allvox.txt')
    all_sel_feat_sumstats.to_csv(outpath, sep='\t', index=False)
    print(f"Saved: {outpath}  | rows={len(all_sel_feat_sumstats)}")
