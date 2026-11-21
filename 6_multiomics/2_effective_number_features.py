#!/usr/bin/env python3
"""
2_effective_number_features.py

Effective number of independent omics features per panel (proteomics,
metabolomics, clinical biomarkers) by the Li & Ji (2005) procedure, computed
on a positive-semidefinite correlation matrix of exactly the features that
were tested, on the participants that entered the association analysis.

Notes
  * The pairwise-complete correlation matrix is not positive-semidefinite, so
    it is projected to the nearest correlation matrix (Higham 2002,
    alternating projections) before the eigendecomposition; the eigenvalues
    then sum to the number of features and no clipping is needed.
  * The effective number is computed over the tested set. A feature that
    cannot be correlated with the rest (too few overlapping participants) is
    kept and treated as independent (r = 0 with every other feature), which
    counts it as one test.
  * Pass the FID/IID list of the participants used for each panel with
    --ids-<panel>; otherwise all participants with data in the panel file are
    used.
  * The participation ratio (sum lambda)^2 / sum lambda^2 and the number of
    principal components explaining 99.5% of the variance are written for
    comparison only.

Outputs (to --out-dir)
  eff_num_features_liji.csv      one row per panel: M_tested, M_eff_LiJi, M_eff_PR,
                                 M_eff_Gao995, n_subjects, n_isolated, min/max eigenvalue
  eigenvalues_<panel>.txt        eigenvalues of the PSD correlation matrix
  features_<panel>.txt           features counted, with a flag for "isolated" ones

The M_eff_LiJi values enter the Bonferroni threshold of 3_screen_top_features.py
(eff_num_feature).

Usage
  python 2_effective_number_features.py \
      --main-dir /path/to/UKB_PSC/multiomics/ \
      --tested-vars-proteomics  FA/proteomics/disc/voxel-assoc/voxel_assoc_1.txt \
      --tested-vars-metabolomics FA/metabolomics/disc/voxel-assoc/voxel_assoc_1.txt \
      --tested-vars-biomarkers  FA/biomarkers/disc/voxel-assoc/voxel_assoc_1.txt \
      [--ids-proteomics ids_prot.txt --ids-metabolomics ids_metab.txt --ids-biomarkers ids_lab.txt]
"""
import os
import sys
import time
import argparse
import numpy as np
import pandas as pd


# ----------------------------------------------------------------------------
# pairwise-complete Pearson correlation via masked matrix products
# ----------------------------------------------------------------------------
def pairwise_complete_corr(X, min_periods=30):
    """X: n x p float array with NaN for missing. Returns (r, N) where r is the
    p x p pairwise-complete correlation (NaN where N < min_periods) and N the
    p x p matrix of pairwise sample sizes."""
    X = np.asarray(X, dtype=np.float64)
    M = np.isfinite(X).astype(np.float64)
    X0 = np.where(M > 0, X, 0.0)
    N = M.T @ M
    S1 = X0.T @ M            # S1[i,j] = sum_k x_ki over k with both i,j observed
    S2 = (X0 ** 2).T @ M     # S2[i,j] = sum_k x_ki^2 over same k
    Sxy = X0.T @ X0
    with np.errstate(invalid="ignore", divide="ignore"):
        cov = Sxy - S1 * S1.T / N
        var_i = S2 - S1 ** 2 / N
        var_j = var_i.T
        r = cov / np.sqrt(var_i * var_j)
    r[N < min_periods] = np.nan
    np.fill_diagonal(r, 1.0)
    return r, N


# ----------------------------------------------------------------------------
# nearest correlation matrix (Higham 2002, alternating projections)
# ----------------------------------------------------------------------------
def nearest_correlation(A, tol=1e-7, max_iter=200, label=""):
    """Project a symmetric matrix with unit diagonal onto the set of
    correlation matrices (PSD, unit diagonal). Returns the PSD matrix and the
    number of iterations used."""
    A = np.asarray(A, dtype=np.float64)
    A = (A + A.T) / 2
    Y = A.copy()
    dS = np.zeros_like(A)
    n = A.shape[0]
    for it in range(1, max_iter + 1):
        R = Y - dS
        w, V = np.linalg.eigh(R)
        w_pos = np.clip(w, 0, None)
        X = (V * w_pos) @ V.T          # projection onto PSD cone
        dS = X - R
        Y_new = X.copy()
        np.fill_diagonal(Y_new, 1.0)   # projection onto unit diagonal
        diff = np.linalg.norm(Y_new - Y, "fro") / max(np.linalg.norm(Y, "fro"), 1e-12)
        Y = Y_new
        if diff < tol:
            break
    w = np.linalg.eigvalsh((Y + Y.T) / 2)
    if w.min() < -1e-8:
        # final tiny clean-up: clip and rescale the trace to n so eigenvalues sum to M
        V = np.linalg.eigh((Y + Y.T) / 2)[1]
        w = np.clip(w, 0, None)
        w *= n / w.sum()
        Y = (V * w) @ V.T
        np.fill_diagonal(Y, 1.0)
    print(f"[{label}] nearest-correlation projection: {it} iteration(s), "
          f"min eigenvalue after projection = {np.linalg.eigvalsh(Y).min():.2e}", flush=True)
    return Y, it


# ----------------------------------------------------------------------------
# effective numbers from an eigenvalue spectrum
# ----------------------------------------------------------------------------
def li_ji(eigvals):
    """Li & Ji (2005) Heredity 95:221: M_eff = sum_i [ I(lambda_i >= 1) + (lambda_i - floor(lambda_i)) ]."""
    lam = np.clip(np.asarray(eigvals, dtype=np.float64), 0, None)
    return float(np.sum((lam >= 1).astype(float) + (lam - np.floor(lam))))


def participation_ratio(eigvals):
    """The estimator used in the submitted version: (sum lambda)^2 / sum lambda^2."""
    lam = np.clip(np.asarray(eigvals, dtype=np.float64), 0, None)
    return float(lam.sum() ** 2 / np.sum(lam ** 2))


def gao_995(eigvals):
    """Gao et al. (2008) Genet Epidemiol 32:361: number of PCs explaining 99.5% of variance."""
    lam = np.sort(np.clip(np.asarray(eigvals, dtype=np.float64), 0, None))[::-1]
    cum = np.cumsum(lam) / lam.sum()
    return int(np.searchsorted(cum, 0.995) + 1)


# ----------------------------------------------------------------------------
# per-panel driver
# ----------------------------------------------------------------------------
def effective_number_panel(df, label, tested_vars=None, ids=None, n_id_cols=2,
                           min_periods=30, out_dir=None):
    t0 = time.time()
    id_cols = list(df.columns[:n_id_cols])
    feat = df.iloc[:, n_id_cols:]

    # 1. restrict to the participants used in the PheWAS, if given
    if ids is not None:
        keep = df[id_cols[-1]].astype(str).isin(set(map(str, ids)))
        print(f"[{label}] restricting to {int(keep.sum())} of {len(df)} participants listed in the ids file")
        feat = feat.loc[keep]

    # 2. restrict to the features that were actually tested
    if tested_vars is not None:
        tested = [c for c in feat.columns if str(c) in tested_vars]
        missing = sorted(set(tested_vars) - set(map(str, feat.columns)))
        if missing:
            print(f"[{label}] [warn] {len(missing)} tested variable(s) not found in the panel file: {missing[:10]}")
        print(f"[{label}] {len(tested)} tested features (of {feat.shape[1]} columns in the panel file)")
        feat = feat[tested]
    elif label == "metabolomics":
        inst1 = [c for c in feat.columns if str(c).endswith("-1.0")]
        if inst1:
            print(f"[{label}] no tested-variable list given: dropping {len(inst1)} instance-1 (repeat-visit) "
                  f"columns, which were not tested; {feat.shape[1] - len(inst1)} baseline measures remain")
            feat = feat.drop(columns=inst1)

    # 3. drop features that could not have been tested (all missing / constant)
    std = feat.std(skipna=True)
    bad = std[(std == 0) | std.isna()].index.tolist()
    if bad:
        print(f"[{label}] dropping {len(bad)} all-missing/zero-variance feature(s): {bad[:10]}")
        feat = feat.drop(columns=bad)
    cols = np.array([str(c) for c in feat.columns])
    n_subj, M = feat.shape
    print(f"[{label}] {n_subj} participants x {M} tested features -- pairwise-complete correlation...", flush=True)

    # 4. pairwise-complete correlation; unresolvable pairs -> r = 0 (treated as independent)
    r, N = pairwise_complete_corr(feat.to_numpy(dtype=np.float64), min_periods=min_periods)
    off = ~np.eye(M, dtype=bool)
    n_unres = int(np.isnan(r[off]).sum() // 2)
    isolated = np.array([np.isnan(r[i, off[i]]).all() for i in range(M)])
    if n_unres:
        print(f"[{label}] {n_unres} feature pair(s) with < {min_periods} overlapping participants set to r = 0; "
              f"{int(isolated.sum())} feature(s) isolated (counted as one independent test each): "
              f"{list(cols[isolated][:10])}")
    r = np.where(np.isnan(r), 0.0, r)

    # 5. project to the nearest correlation matrix (PSD, unit diagonal)
    w_raw = np.linalg.eigvalsh(r)
    print(f"[{label}] raw pairwise-complete matrix: min eigenvalue = {w_raw.min():.3f}, "
          f"{int((w_raw < -1e-8).sum())} negative", flush=True)
    r_psd, _ = nearest_correlation(r, label=label)
    w = np.linalg.eigvalsh(r_psd)[::-1]
    w = np.clip(w, 0, None)

    # 6. effective numbers
    res = dict(panel=label, n_subjects=n_subj, M_tested=M, n_isolated=int(isolated.sum()),
               M_eff_LiJi=li_ji(w), M_eff_PR_old=participation_ratio(w), M_eff_Gao995=gao_995(w),
               lambda_max=float(w[0]), lambda_min=float(w[-1]), sum_lambda=float(w.sum()),
               min_pairwise_N=int(N[off].min()), median_pairwise_N=float(np.median(N[off])))
    print(f"[{label}] DONE in {time.time() - t0:.1f}s -- M = {M}, M_eff (Li & Ji) = {res['M_eff_LiJi']:.1f}, "
          f"participation ratio (old) = {res['M_eff_PR_old']:.1f}, 99.5%-variance PCs = {res['M_eff_Gao995']}", flush=True)

    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
        np.savetxt(os.path.join(out_dir, f"eigenvalues_{label}.txt"), w, fmt="%.6f")
        pd.DataFrame({"feature": cols, "isolated": isolated.astype(int)}).to_csv(
            os.path.join(out_dir, f"features_{label}.txt"), sep="\t", index=False)
    return res


def read_tested_vars(path):
    """Unique VARIABLE values of one discovery voxel-assoc file = the tested feature set."""
    if path is None:
        return None
    v = pd.read_csv(path, sep="\t", usecols=["VARIABLE"])["VARIABLE"].astype(str).unique()
    return set(v)


def read_ids(path):
    if path is None:
        return None
    ids = pd.read_csv(path, sep=None, engine="python", header=None).iloc[:, -1].astype(str)
    return ids.tolist()


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--main-dir", default="/path/to/UKB_PSC/multiomics/")
    ap.add_argument("--out-dir", default=None, help="default: <main-dir>/results/eff_num")
    ap.add_argument("--min-periods", type=int, default=30)
    for p in ["proteomics", "metabolomics", "biomarkers"]:
        ap.add_argument(f"--tested-vars-{p}", default=None,
                        help=f"a discovery voxel-assoc file for {p}; its unique VARIABLE values define the tested set")
        ap.add_argument(f"--ids-{p}", default=None,
                        help=f"file of participant IDs (FID IID or IID) used in the {p} PheWAS")
    args = ap.parse_args()
    out_dir = args.out_dir or os.path.join(args.main_dir, "results", "eff_num")

    files = {"proteomics": "proteomics_dat.txt", "metabolomics": "metabolomics_dat.txt",
             "biomarkers": "biomarkers_dat.txt"}
    rows = []
    for label, fn in files.items():
        print(f"\n=== {label} ===", flush=True)
        df = pd.read_csv(os.path.join(args.main_dir, "data", fn), sep="\t")
        tv = read_tested_vars(getattr(args, f"tested_vars_{label}"))
        ids = read_ids(getattr(args, f"ids_{label}"))
        rows.append(effective_number_panel(df, label, tested_vars=tv, ids=ids,
                                           min_periods=args.min_periods, out_dir=out_dir))
    res = pd.DataFrame(rows)
    os.makedirs(out_dir, exist_ok=True)
    res.to_csv(os.path.join(out_dir, "eff_num_features_liji.csv"), index=False)
    print("\n" + res.to_string(index=False))
    print("\nPaste into eff_num_feature{} of 1_get_top_features_all_omics_allmodalities.py and 2_make_atlas.py:")
    print("eff_num_feature = {" + ", ".join(f"'{r['panel']}': {r['M_eff_LiJi']:.4f}" for r in rows) + "}")


if __name__ == "__main__":
    main()
