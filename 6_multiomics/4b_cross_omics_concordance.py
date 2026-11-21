import os
import numpy as np
import pandas as pd
from tqdm import tqdm

# ----------------------
# Config
# ----------------------
main_dir   = '/path/to/UKB_PSC/multiomics/'
omics      = ['proteomics','metabolomics','biomarkers']
nroi       = 430
alpha      = 0.05
subdir     = 'voxel-assoc/sel_feats'
pattern    = 'voxel_assoc_{i}_sel_feats.txt'
min_n      = 1 # sample size cutoff

# ----------------------
# 1) Build the union of VARIABLES (assuming all present everywhere)
#    If truly identical across files, reading a few is enough.
# ----------------------
all_vars = set()
for typ in omics:
    probe = os.path.join(main_dir, typ, subdir, pattern.format(i=1))
    if os.path.exists(probe):
        df = pd.read_csv(probe, sep='\t', usecols=['VARIABLE', 'N'])
        keep = df.loc[df['N'] >= min_n, 'VARIABLE'].unique()
        all_vars.update(keep)
        
# add genetics
gwas_probe = '/path/to/UKB_PSC_GWAS/FA/voxelgwas/leadsnps/UKB_FA_phase16_1.txt'
gwas_df = pd.read_csv(gwas_probe, sep='\t', usecols=['SNP', 'N'])
gwas_keep = gwas_df.loc[gwas_df['N'] >= min_n, 'SNP'].unique()
all_vars.update(gwas_keep)

variables = sorted(all_vars)
N = len(variables)
var2idx = {v:i for i,v in enumerate(variables)}
assert N > 0, "No variables found."

# ----------------------
# 2) Global accumulators (no pairwise-overlap needed)
# ----------------------
sum_x  = np.zeros(N, dtype=np.float64)
sum_x2 = np.zeros(N, dtype=np.float64)
sum_xy = np.zeros((N, N), dtype=np.float64)  # X^T X over all voxels

n_vox_total = 0  # total rows (voxels) aggregated

# Dice (binary significance)
n_sig         = np.zeros(N, dtype=np.int64)      # per-variable significant voxel counts
intersections = np.zeros((N, N), dtype=np.int64) # pairwise AND of significance

# ----------------------
# FDR helper: BH per column (vectorized over columns)
# ----------------------
def fdr_bh_columns(P, q=0.05):
    """
    P: 2D numpy array (n_vox x m_vars) with NaNs allowed.
    Returns boolean mask S with same shape: True = significant at FDR q (per column).
    BH is applied independently to each column using only finite p-values in that column.
    """
    n_vox, m = P.shape
    S = np.zeros_like(P, dtype=bool)
    for j in range(m):
        pj = P[:, j]
        finite = np.isfinite(pj)
        if not finite.any():
            continue
        p = pj[finite]
        m_eff = p.size
        # sort p
        order = np.argsort(p, kind='mergesort')
        p_sorted = p[order]
        # BH thresholds: (i/m)*q for i=1..m
        i = np.arange(1, m_eff + 1, dtype=np.float64)
        thresh = (i / m_eff) * q
        # find largest i with p_(i) <= thresh_i
        passed = p_sorted <= thresh
        if not passed.any():
            continue
        kmax = np.max(np.where(passed)[0])  # index in 0..m_eff-1
        p_star = p_sorted[kmax]
        # mark significant: p <= p_star among the finite ones
        sig_finite = (p <= p_star)
        out_col = np.zeros(n_vox, dtype=bool)
        out_col[finite] = sig_finite
        S[:, j] = out_col
    return S

# ----------------------
# 3) Stream clusters; per cluster, concat all omics and pivot
# ----------------------
for clust_id in tqdm(range(1, nroi+1), desc="Clusters"):
    frames = []
    for typ in omics:
        f = os.path.join(main_dir, typ, subdir, pattern.format(i=clust_id))
        if os.path.exists(f):
            df = pd.read_csv(f, sep='\t', usecols=['INDEX','VARIABLE','Z','P'])
            if not df.empty:
                frames.append(df)
    # add genetics
    gwas_file = f'/path/to/UKB_PSC_GWAS/FA/voxelgwas/leadsnps/UKB_FA_phase16_{clust_id}.txt'
    if os.path.exists(gwas_file):
        gwas_df = pd.read_csv(gwas_file, sep='\t', usecols=['INDEX','SNP','Z','P'])
        if not gwas_df.empty:
            gwas_df = gwas_df.rename(columns={'SNP': 'VARIABLE'})
            frames.append(gwas_df)
            
    if not frames:
        continue

    dfc = pd.concat(frames, ignore_index=True)

    # Wide matrices: rows = INDEX (voxels), columns = VARIABLE
    Mz = dfc.pivot(index='INDEX', columns='VARIABLE', values='Z')
    Mp = dfc.pivot(index='INDEX', columns='VARIABLE', values='P')

    # Keep variables that are present in THIS cluster; preserve global order
    present_vars = [v for v in variables if v in Mz.columns]
    if not present_vars:
        continue
    Mz = Mz[present_vars]
    Mp = Mp[present_vars]
    cols_idx = np.array([var2idx[v] for v in present_vars], dtype=int)

    # Drop rows with any NaN in Z for present vars (to keep correlation clean)
    if Mz.isnull().values.any():
        Mz = Mz.dropna(axis=0, how='any')
        Mp = Mp.loc[Mz.index]
    if Mz.empty:
        continue

    X = Mz.to_numpy(dtype=np.float64)  # (n_vox_local, m_present)
    P = Mp.to_numpy(dtype=np.float64)

    # Update continuous totals
    n_vox_total += X.shape[0]
    sum_x       += np.bincount(cols_idx, weights=X.sum(axis=0), minlength=N)
    sum_x2      += np.bincount(cols_idx, weights=np.square(X).sum(axis=0), minlength=N)
    sum_xy[np.ix_(cols_idx, cols_idx)] += X.T @ X

    # ---------- FDR (BH) per column within this cluster ----------
    S_loc = fdr_bh_columns(P, q=alpha)                  # boolean (n_vox_local x m_present)
    sig_counts = S_loc.sum(axis=0).astype(np.int64)     # per-variable counts this cluster
    n_sig[cols_idx] += sig_counts
    intersections[np.ix_(cols_idx, cols_idx)] += (S_loc.astype(np.int64).T @ S_loc.astype(np.int64))

# ----------------------
# 4) Pearson correlation (global)
# ----------------------
if n_vox_total == 0:
    raise RuntimeError("No voxels aggregated; check the input paths.")

mean_x = sum_x / n_vox_total
var_x  = sum_x2 - n_vox_total * (mean_x ** 2)
var_x[var_x < 0] = 0.0

den   = np.sqrt(np.outer(var_x, var_x))
num   = sum_xy - np.outer(sum_x, sum_x) / n_vox_total
with np.errstate(invalid='ignore', divide='ignore'):
    corr_matrix = np.where(den > 0, num / den, np.nan)

# ----------------------
# 5) Dice matrix from FDR masks
# ----------------------
den_dice = (n_sig.reshape(-1,1) + n_sig.reshape(1,-1)).astype(np.float64)
dice_matrix = np.full((N, N), np.nan, dtype=np.float64)
mask = den_dice > 0
dice_matrix[mask] = 2.0 * intersections[mask] / den_dice[mask]

# ----------------------
# 6) Save
# ----------------------
outdir = os.path.join(main_dir, 'results')
os.makedirs(outdir, exist_ok=True)

pd.DataFrame(corr_matrix, index=variables, columns=variables)\
  .to_csv(os.path.join(outdir, 'pairwise_variable_correlation_allvox.csv'))

pd.DataFrame(dice_matrix, index=variables, columns=variables)\
  .to_csv(os.path.join(outdir, f'pairwise_variable_dice_allvox_FDRq{alpha:.2f}.csv'))

print(f"Done. Voxels used: {n_vox_total}, variables: {N}, total significant voxels across vars: {int(n_sig.sum())}")