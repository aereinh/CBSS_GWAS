#!/usr/bin/env python3
"""
aggregate_gc_tracts.py

Aggregate the location-level (fiber x along-fiber point) genetic-correlation
estimates to tracts or functional network pairs, for every trait x diffusion
measure x unit cell.

Two tests are computed per cell:

  1. Signed test: a precision-weighted GLS combination of the location-level
     r_g estimates. The covariance accounts for the correlation among
     locations of the same fiber cluster (reconstructed from the fPCA bases
     and the LDR covariance) and treats different clusters as independent.
     It gives a combined z / P-value and a signed aggregate r_g.
  2. Cauchy (ACAT) combination of the same locations' P-values: blind to
     sign and robust to unknown correlation. It is used as a heterogeneity
     companion to the signed test.

Classification per cell after Benjamini-Hochberg FDR across the full grid:
  - signed test significant                   -> "consistent_direction"
  - signed not significant, Cauchy significant -> "mixed_direction"
  - neither significant                        -> "not_significant"

Inputs (paths in PipelineConfig below)
  GC_<trait>_matrix.txt, GC_SE_<trait>_matrix.txt, GC_P_<trait>_matrix.txt
      fibers x points matrices from reshape_genetic_correlation.py
  Fiber_atlas_100points_FullName.csv, cluster_tract_summary.csv
  fPCA bases (*_bases_top*) and LDR covariance (*_ldr_cov_top*) per cluster,
  selLDRs.txt

Assumptions about the input files (--validate-only prints diagnostics for
one measure / cluster / trait without running the full analysis)
  (a) Bases files are (locations, components) with columns ordered by
      descending eigenvalue, as written by HEIG; the leading n_ldrs x n_ldrs
      block of the LDR covariance is used, assuming the same component order.
  (b) Locations within a cluster's bases file are fiber-major (100 points of
      the first fiber, then of the next), in the column order of the GC
      matrices. If this did not hold, only the reconstructed correlation
      among locations would be affected, not the location-level statistics.
  (c) The LDR covariance is proportional to the covariance of the LDR-level
      estimates; only its correlation shape is used. Setting
      use_cluster_correlation = False falls back to independence.
  (d) Cluster IDs in the Clust_Dis30mm column match the index in the file
      names of the bases, LDR covariance and selLDRs files.
  (e) Fiber_index is 1..N contiguous, in the row order of the GC matrices.

Usage
  python aggregate_gc_tracts.py --validate-only
  python aggregate_gc_tracts.py --out tract_assign_gc_results.csv
  python aggregate_gc_tracts.py --group-col network_pair --out network_pair_gc_results.csv
"""

import os
import glob
import time
import argparse
import warnings
from dataclasses import dataclass, field

import numpy as np
import pandas as pd
from scipy.stats import norm
from scipy.linalg import solve, LinAlgError


# ============================================================================
# Config
# ============================================================================

@dataclass
class PipelineConfig:
    modalities: tuple = ("FA", "dti_MD", "dti_L1")
    phenotypes: tuple = (
        "ADHD", "ASD", "CP", "depression", "insomnia", "mdd", "SCZ",
        "ALZ2", "bip", "depressed_affect", "EA", "intelligence",
        "neuroticism", "worry",
    )
    n_points_per_fiber: int = 100

    gc_dir_template: str = "/path/to/UKB_PSC/result_data/{modality}_result/"
    gc_file_template: str = "GC_{phenotype}_matrix.txt"
    gc_se_file_template: str = "GC_SE_{phenotype}_matrix.txt"
    gc_p_file_template: str = "GC_P_{phenotype}_matrix.txt"

    fiber_cluster_csv: str = "/path/to/UKB_PSC/Fiber_atlas_100points_FullName.csv"
    fiber_index_col: str = "Fiber_index"
    cluster_col: str = "Clust_Dis30mm"
    cluster_tract_csv: str = "/path/to/UKB_PSC/cluster_tract_summary.csv"
    cluster_tract_cluster_col: str = "cluster"
    cluster_tract_tract_col: str = "tract_assign"   # which cluster_tract_csv column defines the
                                                     # grouping unit -- e.g. "tract_assign" or
                                                     # "network_pair". Whatever's chosen here is
                                                     # also the output column name (see
                                                     # build_results_table()).

    bases_dir_template: dict = field(default_factory=lambda: {
        "FA": "/path/to/UKB_PSC/output/output_FA/fpca/",
        "dti_MD": "/path/to/UKB_PSC/output/output_dti_MD/fpca/",
        "dti_L1": "/path/to/UKB_PSC/output/output_dti_L1/fpca/",
    })
    ldrcov_dir_template: dict = field(default_factory=lambda: {
        "FA": "/path/to/UKB_PSC/output/output_FA/ldrs/",
        "dti_MD": "/path/to/UKB_PSC/output/output_dti_MD/ldrs/",
        "dti_L1": "/path/to/UKB_PSC/output/output_dti_L1/ldrs/",
    })
    bases_filename_glob: str = "UKB_{modality}_phase16_{idx}_bases_top*"
    ldrcov_filename_glob: str = "UKB_{modality}_phase16_{idx}_ldr_cov_top*"
    sel_ldrs_template: str = "/path/to/UKB_PSC/output/output_{modality}/selLDRs.txt"

    use_cluster_correlation: bool = True
    corr_shrinkage: float = 0.10
    fdr_alpha: float = 0.05
    min_fibers_per_tract: int = 1
    tracts: tuple = None   # None -> use every group found in cluster_tract_csv's cluster_tract_tract_col
                            # column (whatever that's set to); else restrict to these label values

    output_dir: str = "/path/to/UKB_PSC/result_data/agg/"


# ============================================================================
# Loading: fiber -> cluster -> tract mapping
# ============================================================================

def load_fiber_tract_mapping(cfg: PipelineConfig) -> pd.DataFrame:
    """One row per fiber: [fiber_index, cluster, tract]."""
    fiber_clust = pd.read_csv(cfg.fiber_cluster_csv)
    missing = {cfg.fiber_index_col, cfg.cluster_col} - set(fiber_clust.columns)
    if missing:
        raise ValueError(f"{cfg.fiber_cluster_csv} missing columns {missing} "
                          f"(has: {list(fiber_clust.columns)})")
    fiber_clust = fiber_clust[[cfg.fiber_index_col, cfg.cluster_col]].rename(
        columns={cfg.fiber_index_col: "fiber_index", cfg.cluster_col: "cluster"})

    clust_tract = pd.read_csv(cfg.cluster_tract_csv)
    missing = {cfg.cluster_tract_cluster_col, cfg.cluster_tract_tract_col} - set(clust_tract.columns)
    if missing:
        raise ValueError(f"{cfg.cluster_tract_csv} missing columns {missing} "
                          f"(has: {list(clust_tract.columns)})")
    clust_tract = clust_tract[[cfg.cluster_tract_cluster_col, cfg.cluster_tract_tract_col]].rename(
        columns={cfg.cluster_tract_cluster_col: "cluster", cfg.cluster_tract_tract_col: "tract"})

    merged = fiber_clust.merge(clust_tract, on="cluster", how="left")
    n_unmapped = merged["tract"].isna().sum()
    if n_unmapped:
        warnings.warn(f"{n_unmapped}/{len(merged)} fibers' clusters have no tract assignment.")
    return merged.sort_values("fiber_index").reset_index(drop=True)


def build_cluster_fiber_order(mapping: pd.DataFrame) -> dict:
    """{cluster_id: ascending-sorted list of fiber_index in that cluster}.
    ASSUMPTION (b): this order matches the bases file's row layout."""
    return {c: sorted(g["fiber_index"].tolist()) for c, g in mapping.groupby("cluster")}


def validate_inputs(cfg: PipelineConfig, mapping: pd.DataFrame, n_fibers_in_gc_matrix: int):
    print(f"[validate] mapping: {len(mapping)} fibers, {mapping['cluster'].nunique()} clusters, "
          f"{mapping['tract'].nunique()} {cfg.cluster_tract_tract_col} group(s)")
    if len(mapping) != n_fibers_in_gc_matrix:
        warnings.warn(f"Mapping has {len(mapping)} fibers, GC matrix has {n_fibers_in_gc_matrix} "
                       f"rows -- MUST match for the row-alignment assumption (e) to hold.")
    fi = np.sort(mapping["fiber_index"].unique())
    contiguous = np.array_equal(fi, np.arange(fi.min(), fi.min() + len(fi)))
    print(f"[validate] fiber_index range {fi.min()}-{fi.max()}, contiguous: {contiguous} "
          f"(assumption (e) -- if False, fiber_index_to_row must be rebuilt from an explicit "
          f"row-number column instead of assumed ascending order)")
    if mapping["fiber_index"].duplicated().any():
        warnings.warn(f"{mapping['fiber_index'].duplicated().sum()} duplicate fiber_index values.")


# ============================================================================
# Loading: GC / SE / P matrices
# ============================================================================

def load_gc_triplet(cfg: PipelineConfig, modality: str, phenotype: str):
    d = cfg.gc_dir_template.format(modality=modality)
    paths = {
        "gc": os.path.join(d, cfg.gc_file_template.format(phenotype=phenotype)),
        "se": os.path.join(d, cfg.gc_se_file_template.format(phenotype=phenotype)),
        "p": os.path.join(d, cfg.gc_p_file_template.format(phenotype=phenotype)),
    }
    for p in paths.values():
        if not os.path.exists(p):
            print(f"[skip] {modality} x {phenotype}: missing {p}")
            return None
    gc, se, p_ = (np.loadtxt(paths[k]) for k in ("gc", "se", "p"))
    if not (gc.shape == se.shape == p_.shape):
        raise ValueError(f"{modality} x {phenotype}: shape mismatch GC{gc.shape} SE{se.shape} P{p_.shape}")
    if gc.shape[1] != cfg.n_points_per_fiber:
        raise ValueError(f"{modality} x {phenotype}: expected {cfg.n_points_per_fiber} cols, got {gc.shape[1]}")
    return gc, se, p_


# ============================================================================
# Loading: per-cluster bases / ldr-cov / n_ldrs (cached across phenotypes --
# these are imaging-side only, independent of the y2/external trait)
# ============================================================================

def load_sel_ldrs(cfg: PipelineConfig, modality: str) -> dict:
    path = cfg.sel_ldrs_template.format(modality=modality)
    df = pd.read_csv(path)
    if df.shape[1] < 2:
        raise ValueError(f"{path}: expected >=2 columns (idx, n_ldrs), got {df.shape[1]}")
    idx_col, nldr_col = df.columns[0], df.columns[1]
    return dict(zip(df[idx_col].astype(int), df[nldr_col].astype(int)))


def _find_one(pattern: str) -> str:
    matches = glob.glob(pattern)
    if len(matches) != 1:
        raise FileNotFoundError(f"Expected exactly 1 file matching {pattern}, found {len(matches)}: {matches}")
    return matches[0]


def _load_array(path: str) -> np.ndarray:
    """Load a bases/ldr-cov matrix, auto-detecting numpy binary (.npy/.npz)
    vs plain-text formats.

    heig apparently writes these out with np.save rather than as
    whitespace-delimited text: real files start with the byte 0x93, which
    is the first byte of numpy's .npy magic string b'\\x93NUMPY' (not valid
    UTF-8 on its own, which is why np.loadtxt raised UnicodeDecodeError).
    """
    with open(path, "rb") as fh:
        magic = fh.read(6)
    if magic == b"\x93NUMPY":
        return np.load(path, allow_pickle=False)
    if magic[:2] == b"PK":  # zip container -> likely .npz
        loaded = np.load(path, allow_pickle=False)
        if hasattr(loaded, "files"):  # NpzFile: expect exactly one array inside
            keys = list(loaded.files)
            if len(keys) != 1:
                raise ValueError(f"{path}: .npz contains {len(keys)} arrays {keys}, expected exactly 1 "
                                  f"-- update _load_array() to select the right key.")
            return loaded[keys[0]]
        return loaded
    # fall back to plain text (covers the case where some files really are text)
    return np.loadtxt(path)


def load_bases(cfg, modality, idx, n_targets_expected, n_ldrs) -> np.ndarray:
    """Bases are saved by HEIG as (locations, components), with columns
    ordered by descending eigenvalue, so the first n_ldrs columns are the
    selected LDRs."""
    d = cfg.bases_dir_template[modality]
    path = _find_one(os.path.join(d, cfg.bases_filename_glob.format(modality=modality, idx=idx)))
    B = _load_array(path)
    if B.ndim == 1:
        B = B.reshape(-1, 1)
    if B.shape[0] != n_targets_expected:
        raise ValueError(f"cluster {idx} ({modality}): bases shape {B.shape} -- expected "
                          f"{n_targets_expected} rows (targets), per heig's FPCA save convention "
                          f"(bases.shape == (n_voxels, n_top)).")
    if B.shape[1] < n_ldrs:
        raise ValueError(f"cluster {idx} ({modality}): bases has only {B.shape[1]} retained "
                          f"components, fewer than n_ldrs={n_ldrs} selected -- can't select that many LDRs.")
    return B[:, :n_ldrs]


def load_ldr_cov(cfg, modality, idx, n_ldrs) -> np.ndarray:
    """Same story as load_bases(): the ldr-cov file holds the covariance
    among ALL retained components (K x K, K = the "top*" count), not just
    the n_ldrs actually selected. Since selected LDRs are the first n_ldrs
    components in ranked order, their covariance is exactly the leading
    n_ldrs x n_ldrs principal submatrix -- so we slice rather than requiring
    an exact shape match."""
    d = cfg.ldrcov_dir_template[modality]
    path = _find_one(os.path.join(d, cfg.ldrcov_filename_glob.format(modality=modality, idx=idx)))
    S = _load_array(path)
    if S.ndim != 2 or S.shape[0] != S.shape[1]:
        raise ValueError(f"cluster {idx} ({modality}): ldr-cov shape {S.shape} is not square")
    if S.shape[0] < n_ldrs:
        raise ValueError(f"cluster {idx} ({modality}): ldr-cov shape {S.shape} has fewer than "
                          f"n_ldrs={n_ldrs} retained components -- can't select that many LDRs.")
    return S[:n_ldrs, :n_ldrs]


class ClusterBasisCache:
    LOG_EVERY = 20   # print a [cache] line every this-many new disk loads, so a slow first
                      # combo shows live progress instead of a silent multi-minute gap.

    def __init__(self, cfg: PipelineConfig):
        self.cfg = cfg
        self._sel_ldrs = {}
        self._bases = {}
        self._ldrcov = {}
        self.n_loads = 0   # count of (modality, cluster) pairs actually read from disk so far --
                            # used for progress reporting (a big jump here means a "cold" modality).
        self._t_start = time.time()

    def n_ldrs(self, modality: str, idx) -> int:
        if modality not in self._sel_ldrs:
            self._sel_ldrs[modality] = load_sel_ldrs(self.cfg, modality)
        return self._sel_ldrs[modality][int(idx)]

    def get(self, modality: str, idx, fiber_order: list):
        key = (modality, idx)
        if key not in self._bases:
            n_ldrs = self.n_ldrs(modality, idx)
            n_targets = len(fiber_order) * self.cfg.n_points_per_fiber
            self._bases[key] = load_bases(self.cfg, modality, idx, n_targets, n_ldrs)
            self._ldrcov[key] = load_ldr_cov(self.cfg, modality, idx, n_ldrs)
            self.n_loads += 1
            if self.n_loads == 1 or self.n_loads % self.LOG_EVERY == 0:
                elapsed = time.time() - self._t_start
                print(f"[cache] {self.n_loads} cluster file-pair(s) loaded from disk so far "
                      f"(most recent: {modality} cluster {idx}) -- {elapsed:.1f}s elapsed since "
                      f"first load, {elapsed / self.n_loads:.2f}s/load avg.", flush=True)
        return self._bases[key], self._ldrcov[key]


# ============================================================================
# Cross-point correlation reconstruction + GLS combination
# ============================================================================

def cluster_correlation_block(B_sub: np.ndarray, Sigma_LDR: np.ndarray, shrinkage: float) -> np.ndarray:
    """Correlation matrix implied by the shared FPCA basis, shrunk toward
    identity. Shrinkage is necessary, not cosmetic: Sigma_LDR is only
    n_ldrs x n_ldrs, so B_sub @ Sigma_LDR @ B_sub.T is rank <= n_ldrs --
    singular whenever a block has more points than retained components
    (the typical case). Shrinkage reflects real idiosyncratic point-level
    noise beyond what the truncated basis captures, and makes the block
    invertible."""
    Omega = B_sub @ Sigma_LDR @ B_sub.T
    d = np.sqrt(np.clip(np.diag(Omega), 1e-12, None))
    R = Omega / np.outer(d, d)
    np.fill_diagonal(R, 1.0)
    n = R.shape[0]
    return (1 - shrinkage) * R + shrinkage * np.eye(n)


def cluster_block_gls(df_cluster_valid: pd.DataFrame, fiber_pos: dict, B_full: np.ndarray,
                       Sigma_LDR: np.ndarray, n_points: int, shrinkage: float):
    """df_cluster_valid: rows for ONE cluster within ONE tract, valid points
    only, columns fiber_index/point_idx/rg/se. Returns this block's
    contribution to the tract-level GLS numerator (1^T Sigma^-1 r) and
    denominator (1^T Sigma^-1 1)."""
    fiber_pos_arr = df_cluster_valid["fiber_index"].map(fiber_pos).to_numpy()
    if np.isnan(fiber_pos_arr).any():
        raise KeyError("Some fibers in this cluster/tract block aren't in cluster_fiber_order -- "
                       "mapping inconsistency, stop and check inputs.")
    point_idx_arr = df_cluster_valid["point_idx"].to_numpy()
    basis_row = (fiber_pos_arr * n_points + point_idx_arr).astype(int)

    B_sub = B_full[basis_row, :]
    R = cluster_correlation_block(B_sub, Sigma_LDR, shrinkage)

    se_vec = df_cluster_valid["se"].to_numpy()
    rg_vec = df_cluster_valid["rg"].to_numpy()
    Sigma_rg = (se_vec[:, None] * R) * se_vec[None, :]   # D_se @ R @ D_se

    ones = np.ones(len(rg_vec))
    try:
        x_r = solve(Sigma_rg, rg_vec, assume_a="pos")
        x_1 = solve(Sigma_rg, ones, assume_a="pos")
    except LinAlgError:
        ridge = 1e-8 * np.mean(se_vec ** 2) * np.eye(len(rg_vec))
        x_r = solve(Sigma_rg + ridge, rg_vec, assume_a="pos")
        x_1 = solve(Sigma_rg + ridge, ones, assume_a="pos")

    return float(ones @ x_r), float(ones @ x_1)


def cauchy_combine(pvals: np.ndarray) -> float:
    """ACAT / Cauchy combination -- sign-blind by construction, robust to
    unknown correlation among the inputs (this is exactly why it needs no
    Sigma). Used here as the heterogeneity companion to the signed test."""
    pvals = np.clip(np.asarray(pvals, dtype=float), 1e-300, 1 - 1e-16)
    stat = np.mean(np.tan((0.5 - pvals) * np.pi))
    return float(np.clip(0.5 - np.arctan(stat) / np.pi, 0.0, 1.0))


def bh_fdr(pvals: np.ndarray) -> np.ndarray:
    """Plain Benjamini-Hochberg, no extra dependency needed."""
    pvals = np.asarray(pvals, dtype=float)
    n = len(pvals)
    order = np.argsort(pvals)
    ranked = pvals[order]
    q = ranked * n / np.arange(1, n + 1)
    q = np.minimum.accumulate(q[::-1])[::-1]
    out = np.empty(n)
    out[order] = np.clip(q, 0, 1)
    return out


# ============================================================================
# Per-cell aggregation
# ============================================================================

def build_tract_long_frame(fiber_ids, fiber_index_to_row, fiber_to_cluster, GC, SE, P, n_points) -> pd.DataFrame:
    rows = []
    for f in fiber_ids:
        r = fiber_index_to_row[f]
        rows.append(pd.DataFrame({
            "fiber_index": f,
            "point_idx": np.arange(n_points),
            "cluster": fiber_to_cluster[f],
            "rg": GC[r, :],
            "se": SE[r, :],
            "p": P[r, :],
        }))
    df = pd.concat(rows, ignore_index=True)
    df["valid"] = np.isfinite(df["rg"]) & np.isfinite(df["se"]) & (df["se"] > 0) & np.isfinite(df["p"])
    return df


def aggregate_tract_cell(fiber_ids, fiber_index_to_row, fiber_to_cluster, GC, SE, P,
                          cluster_fiber_order, basis_cache, modality, cfg: PipelineConfig) -> dict:
    df = build_tract_long_frame(fiber_ids, fiber_index_to_row, fiber_to_cluster, GC, SE, P, cfg.n_points_per_fiber)
    df_valid = df[df["valid"]]
    n_points_total = len(df_valid)

    base = dict(n_fibers=len(fiber_ids), n_points=n_points_total, n_clusters=0,
                signed_rg=np.nan, signed_se=np.nan, signed_z=np.nan, signed_p=np.nan, cauchy_p=np.nan)
    if n_points_total == 0:
        return base

    cauchy_p = cauchy_combine(df_valid["p"].to_numpy())

    if not cfg.use_cluster_correlation:
        w = 1.0 / df_valid["se"].to_numpy() ** 2
        num, den = float((w * df_valid["rg"].to_numpy()).sum()), float(w.sum())
    else:
        num, den = 0.0, 0.0
        for c, df_c in df_valid.groupby("cluster"):
            fiber_order_list = cluster_fiber_order[c]
            fiber_pos = {fid: i for i, fid in enumerate(fiber_order_list)}
            B_full, Sigma_LDR = basis_cache.get(modality, c, fiber_order_list)
            num_c, den_c = cluster_block_gls(df_c, fiber_pos, B_full, Sigma_LDR, cfg.n_points_per_fiber, cfg.corr_shrinkage)
            num += num_c
            den += den_c

    rg_agg = num / den
    se_agg = 1.0 / np.sqrt(den)
    z_agg = rg_agg / se_agg
    p_agg = float(2 * (1 - norm.cdf(abs(z_agg))))

    base.update(n_clusters=int(df_valid["cluster"].nunique()), signed_rg=rg_agg, signed_se=se_agg,
                signed_z=z_agg, signed_p=p_agg, cauchy_p=cauchy_p)
    return base


# ============================================================================
# Driver
# ============================================================================

def build_results_table(cfg: PipelineConfig) -> pd.DataFrame:
    mapping = load_fiber_tract_mapping(cfg)
    cluster_fiber_order = build_cluster_fiber_order(mapping)
    fiber_to_cluster = dict(zip(mapping["fiber_index"], mapping["cluster"]))
    tract_to_fibers = mapping.dropna(subset=["tract"]).groupby("tract")["fiber_index"].apply(list).to_dict()
    if cfg.tracts is not None:
        missing = set(cfg.tracts) - set(tract_to_fibers)
        if missing:
            warnings.warn(f"cfg.tracts requested labels not found in the mapping file: {missing}")
        tract_to_fibers = {t: fibers for t, fibers in tract_to_fibers.items() if t in cfg.tracts}
        print(f"[run] restricted to {len(tract_to_fibers)} {cfg.cluster_tract_tract_col} group(s): "
              f"{list(tract_to_fibers)}")
    basis_cache = ClusterBasisCache(cfg)

    # ---- progress tracking ---------------------------------------------
    # Each modality x phenotype combo prints one [progress] line when done,
    # with its own wall time, how many NEW cluster bases/ldr-cov file-pairs
    # it had to read from disk (vs. reused from cache), running elapsed
    # time, and a rough ETA. The first phenotype of each modality is
    # expected to be much slower than the rest (it's the one paying the
    # disk-I/O cost for that modality's ~430 clusters) -- n_new_loads makes
    # that visible directly instead of having to infer it from timing alone.
    total_combos = len(cfg.modalities) * len(cfg.phenotypes)
    combo_i = 0
    t_pipeline_start = time.time()

    rows = []
    first_pass = True
    for modality in cfg.modalities:
        for phenotype in cfg.phenotypes:
            combo_i += 1
            t_combo_start = time.time()

            loaded = load_gc_triplet(cfg, modality, phenotype)
            if loaded is None:
                print(f"[progress] {combo_i}/{total_combos} {modality} x {phenotype}: skipped (missing GC/SE/P files)")
                continue
            GC, SE, P = loaded
            if first_pass:
                validate_inputs(cfg, mapping, GC.shape[0])
                first_pass = False

            # ASSUMPTION (e): fiber_index is 1..N contiguous, ascending order
            # matches GC/SE/P row order. Rebuild from an explicit row column
            # instead if validate_inputs() flagged this as False.
            fiber_index_to_row = {int(fid): i for i, fid in enumerate(sorted(mapping["fiber_index"].unique()))}

            n_loads_before = basis_cache.n_loads
            print(f"[run] {modality} x {phenotype}: {len(tract_to_fibers)} {cfg.cluster_tract_tract_col} group(s)")
            for tract, fiber_ids in tract_to_fibers.items():
                if len(fiber_ids) < cfg.min_fibers_per_tract:
                    continue
                res = aggregate_tract_cell(fiber_ids, fiber_index_to_row, fiber_to_cluster, GC, SE, P,
                                            cluster_fiber_order, basis_cache, modality, cfg)
                res.update(tract=tract, phenotype=phenotype, modality=modality)
                rows.append(res)

            t_combo = time.time() - t_combo_start
            t_elapsed = time.time() - t_pipeline_start
            n_new_loads = basis_cache.n_loads - n_loads_before
            avg_per_combo = t_elapsed / combo_i
            eta_min = avg_per_combo * (total_combos - combo_i) / 60
            print(f"[progress] {combo_i}/{total_combos} done -- {modality} x {phenotype} took {t_combo:.1f}s "
                  f"({n_new_loads} new cluster file-pair(s) loaded from disk this combo). "
                  f"Elapsed {t_elapsed / 60:.1f} min, avg {avg_per_combo:.1f}s/combo, "
                  f"ETA ~{eta_min:.1f} min remaining.")

    print(f"[progress] all {combo_i} combos finished in {(time.time() - t_pipeline_start) / 60:.1f} min total.")

    df = pd.DataFrame(rows)

    df["signed_q"], df["cauchy_q"] = np.nan, np.nan
    m_signed = df["signed_p"].notna()
    m_cauchy = df["cauchy_p"].notna()
    df.loc[m_signed, "signed_q"] = bh_fdr(df.loc[m_signed, "signed_p"].to_numpy())
    df.loc[m_cauchy, "cauchy_q"] = bh_fdr(df.loc[m_cauchy, "cauchy_p"].to_numpy())

    df["classification"] = "not_significant"
    df.loc[df["signed_q"] < cfg.fdr_alpha, "classification"] = "consistent_direction"
    mixed = (df["classification"] == "not_significant") & (df["cauchy_q"] < cfg.fdr_alpha)
    df.loc[mixed, "classification"] = "mixed_direction"

    cols = ["tract", "phenotype", "modality", "n_fibers", "n_points", "n_clusters",
            "signed_rg", "signed_se", "signed_z", "signed_p", "signed_q",
            "cauchy_p", "cauchy_q", "classification"]
    # Output column takes the NAME of whichever grouping column was actually used
    # (cfg.cluster_tract_tract_col, e.g. "tract_assign" or "network_pair") rather
    # than a generic "tract" label, so the CSV is self-documenting about which
    # mapping produced it.
    return df[cols].rename(columns={"tract": cfg.cluster_tract_tract_col})


# ============================================================================
# Validation-only diagnostic (run this first)
# ============================================================================

def run_validation_only(cfg: PipelineConfig):
    mapping = load_fiber_tract_mapping(cfg)
    print(mapping.head())
    cluster_fiber_order = build_cluster_fiber_order(mapping)

    modality = cfg.modalities[0]
    phenotype = cfg.phenotypes[0]
    loaded = load_gc_triplet(cfg, modality, phenotype)
    if loaded is None:
        print(f"Could not load GC/SE/P for {modality} x {phenotype} -- fix paths first.")
        return
    GC, SE, P = loaded
    print(f"[validate] GC/SE/P shape: {GC.shape} (expect (n_fibers, {cfg.n_points_per_fiber}))")
    validate_inputs(cfg, mapping, GC.shape[0])

    first_cluster = sorted(cluster_fiber_order.keys())[0]
    fiber_order_list = cluster_fiber_order[first_cluster]
    print(f"[validate] cluster {first_cluster}: {len(fiber_order_list)} fibers")

    cache = ClusterBasisCache(cfg)
    n_ldrs = cache.n_ldrs(modality, first_cluster)
    print(f"[validate] n_ldrs for cluster {first_cluster} ({modality}): {n_ldrs}")
    B, Sigma_LDR = cache.get(modality, first_cluster, fiber_order_list)
    print(f"[validate] bases shape: {B.shape} (expect ({len(fiber_order_list) * cfg.n_points_per_fiber}, {n_ldrs}))")
    print(f"[validate] ldr-cov shape: {Sigma_LDR.shape} (expect ({n_ldrs}, {n_ldrs}))")

    # Spot-check: reconstructed correlation between two adjacent points on
    # the SAME fiber should typically be higher than between points on
    # different fibers, if assumption (b)'s ordering is right and the
    # underlying signal is spatially smooth. Not a proof, but a red flag if
    # it looks backwards.
    R_full = cluster_correlation_block(B, Sigma_LDR, cfg.corr_shrinkage)
    same_fiber_corr = R_full[0, 1]          # points 0 and 1 of the first fiber (if ordering is right)
    diff_fiber_corr = R_full[0, cfg.n_points_per_fiber]  # point 0 of fiber 1 vs point 0 of fiber 2
    print(f"[validate] reconstructed corr, same-fiber-adjacent-points: {same_fiber_corr:.3f}")
    print(f"[validate] reconstructed corr, different-fiber-same-position: {diff_fiber_corr:.3f}")
    print("[validate] if same-fiber is NOT clearly higher than different-fiber, "
          "assumption (b) (row ordering) is probably wrong -- stop and check with the data owner "
          "of the bases files before trusting the full pipeline.")


# ============================================================================
if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--validate-only", action="store_true")
    parser.add_argument("--out", default=None, help="override output CSV path")
    parser.add_argument("--modalities", default=None,
                         help="comma-separated subset, e.g. FA or FA,dti_MD")
    parser.add_argument("--phenotypes", default=None,
                         help="comma-separated subset, e.g. ADHD,SCZ")
    parser.add_argument("--group-col", default=None,
                         help="which cluster_tract_csv column defines the grouping unit, "
                              "e.g. tract_assign (default) or network_pair. This also becomes "
                              "the output CSV's grouping column name.")
    parser.add_argument("--tracts", default=None,
                         help="comma-separated subset of group labels (from whichever column "
                              "--group-col selects), e.g. CST.R,CST.L or specific network_pair values")
    args = parser.parse_args()

    cfg = PipelineConfig()
    if args.modalities:
        cfg.modalities = tuple(args.modalities.split(","))
    if args.phenotypes:
        cfg.phenotypes = tuple(args.phenotypes.split(","))
    if args.group_col:
        cfg.cluster_tract_tract_col = args.group_col
    if args.tracts:
        cfg.tracts = tuple(args.tracts.split(","))

    if args.validate_only:
        run_validation_only(cfg)
    else:
        df = build_results_table(cfg)
        os.makedirs(cfg.output_dir, exist_ok=True)
        out_path = args.out or os.path.join(cfg.output_dir, "tract_gc_results.csv")
        df.to_csv(out_path, index=False)
        print(f"Wrote {len(df)} rows to {out_path}")
        print(df["classification"].value_counts())