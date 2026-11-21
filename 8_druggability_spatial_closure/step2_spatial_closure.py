"""
step2_spatial_closure.py

Spatial closure: does the white-matter footprint of a colocalized locus overlap
the genetic-correlation (r_g) map of a disease more than expected by chance?

Definitions
  * Locus footprint: locations with P < 9.12e-12 for the lead variant
    (5e-8 / M_eff, M_eff = 5479.5), the threshold used for the locus atlases.
  * Disease map: locations with r_g P < 9.12e-6 (0.05 / M_eff).
  * P == 0 entries are set to NaN before thresholding (a coding artefact at
    1,833 MD locations on 19 fibers that is identical across SNP files).
  * Fold enrichment = observed overlap / overlap expected from the two
    fractions; its ceiling (1 / footprint fraction) is also reported.
  * Rotation null: each draw circularly shifts every fiber's footprint along
    its 100 points by a random offset and permutes fiber identity within
    tract, which preserves footprint size and along-fiber autocorrelation.
    The empirical one-sided P and the null mean / 95th percentile are reported.
  * All SNPs of the colocalization table are tested against all disease maps;
    the primary test of a locus is the map of its own colocalized trait.
  * Sensitivity analyses: top-k locations by -log10 P (k = 500, 2000, 10000)
    instead of the footprint threshold; BH-FDR 0.05 thresholds on both sides;
    and a threshold-free disease side (mean -log10 P of the r_g map inside the
    footprint, gc_mean_nlp_fp, against the same rotation null), for disease
    maps without any location below the map threshold.

Usage
  python step2_spatial_closure.py --coloc-table coloc_hits_strong_all_pairs.tsv \
      --modalities FA --n-rot 10000 --seed 1 --out closure_FA.csv

  --snps rs1,rs2 restricts to given SNPs; --no-null skips the rotation null;
  --base overrides the data root (BASE below).
"""

import argparse
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd

# ---------------------------------------------------------------- constants

BASE = Path("/path/to/UKB_PSC")
M_EFF = 5479.479
P_FOOTPRINT = 5e-8 / M_EFF        # 9.12e-12
P_GC = 0.05 / M_EFF               # 9.12e-6
N_POINTS = 100                    # points per fiber
N_FIBERS = 6090
TOPK = (500, 2000, 10000)

DISEASES = [
    "EA", "ADHD", "insomnia", "ALZ2", "intelligence", "ASD", "mdd", "bip",
    "neuroticism", "CP", "SCZ", "depressed_affect", "worry", "depression",
]

# coloc trait label -> r_g panel code (None = not in the panel)
TRAIT_TO_GC = {
    "AD": "ALZ2", "SCZ": "SCZ", "EA": "EA", "intelligence": "intelligence",
    "CP": "CP", "BD": "bip", "ADHD": "ADHD", "insomina": "insomnia",
    "insomnia": "insomnia", "MDD": "mdd", "ASD": "ASD",
}


# ---------------------------------------------------------------- paths

def voxelgwas_path(modality, snp):
    return BASE / f"output/output_{modality}/voxelgwas/colocsnps/agg/{snp}_pmat.npy"


def gc_path(modality, disease):
    return BASE / f"herigc/result_data/{modality}_result/GC_P_{disease}_matrix.txt"


def fiber_tract_path():
    return BASE / "output/output_FA/voxelgwas/voxel_to_tract.csv"


# ---------------------------------------------------------------- loading

def load_pmap(path):
    a = np.asarray(np.load(path), dtype=float)
    if a.ndim == 1:
        a = a.reshape(N_FIBERS, N_POINTS)
    a = a.copy()
    a[a == 0] = np.nan            # artefact: P exactly 0 is not an association
    return a


def load_gc_pmap(path):
    try:
        a = np.loadtxt(path)
    except ValueError:
        a = pd.read_csv(path, sep=None, engine="python", header=None).to_numpy(dtype=float)
    a = np.asarray(a, dtype=float)
    if a.ndim == 1:
        a = a.reshape(N_FIBERS, N_POINTS)
    a = a.copy()
    a[a == 0] = np.nan
    return a


def load_fiber_tracts(path):
    if not path.exists():
        print(f"WARNING: {path} not found; rotation null will permute fibers "
              f"across the whole skeleton instead of within tract", file=sys.stderr)
        return np.array(["ALL"] * N_FIBERS)
    df = pd.read_csv(path).sort_values("fiber_idx").reset_index(drop=True)
    t = df["tract"].astype(str).to_numpy()
    assert len(t) == N_FIBERS, f"expected {N_FIBERS} fibers, got {len(t)}"
    return t


def bh_mask(p, alpha=0.05):
    sig = np.zeros(p.shape, dtype=bool)
    v = np.isfinite(p)
    if not v.any():
        return sig
    pv = p[v]
    n = pv.size
    order = np.argsort(pv)
    ranked = pv[order]
    passed = ranked <= (np.arange(1, n + 1) / n) * alpha
    if passed.any():
        cutoff = ranked[np.max(np.where(passed)[0])]
        sig[v] = pv <= cutoff
    return sig


# ---------------------------------------------------------------- rotation null

class RotationNull:
    """Circular shift along each fiber + permutation of fibers within tract."""

    def __init__(self, fiber_tracts, n_draws, seed=1):
        self.rng = np.random.default_rng(seed)
        self.n_draws = n_draws
        self.groups = [np.where(fiber_tracts == t)[0] for t in np.unique(fiber_tracts)]
        self.row = np.arange(N_FIBERS)[:, None]
        self.col = np.arange(N_POINTS)[None, :]

    def draw(self, mask2d):
        shifts = self.rng.integers(0, N_POINTS, size=N_FIBERS)
        idx = (self.col - shifts[:, None]) % N_POINTS
        shifted = mask2d[self.row, idx]
        perm = np.arange(N_FIBERS)
        for g in self.groups:
            perm[g] = self.rng.permutation(g)
        return shifted[perm]


def null_overlaps(null, fp2d, gc_masks, valid2d, gc_nlp=None):
    """Return (overlap counts (n_draws, G), mean -log10 P inside footprint (n_draws, G) or None)."""
    out = np.zeros((null.n_draws, len(gc_masks)), dtype=np.int64)
    gc_stack = np.stack([g & valid2d for g in gc_masks])          # (G, 6090, 100)
    nlp_out = None
    if gc_nlp is not None:
        nlp_stack = np.stack(gc_nlp)                                # (G, 6090, 100), NaN outside valid
        nlp_out = np.zeros((null.n_draws, len(gc_masks)), dtype=float)
    for b in range(null.n_draws):
        f = null.draw(fp2d) & valid2d
        out[b] = (gc_stack & f).sum(axis=(1, 2))
        if nlp_out is not None:
            nlp_out[b] = np.nanmean(nlp_stack[:, f], axis=1)
    return out, nlp_out


# ---------------------------------------------------------------- main

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--coloc-table", required=True,
                    help="coloc_hits_strong_all_pairs.tsv (columns snp, trait, PP.H4.abf, ...)")
    ap.add_argument("--modalities", default="FA", help="comma-separated: FA,dti_MD,dti_L1")
    ap.add_argument("--snps", default=None, help="comma-separated subset of SNPs")
    ap.add_argument("--n-rot", type=int, default=1000)
    ap.add_argument("--no-null", action="store_true")
    ap.add_argument("--seed", type=int, default=1)
    ap.add_argument("--out", default="closure_results.csv")
    ap.add_argument("--base", default=None, help="override BASE directory")
    args = ap.parse_args()

    global BASE
    if args.base:
        BASE = Path(args.base)

    modalities = [m.strip() for m in args.modalities.split(",") if m.strip()]

    coloc = pd.read_csv(args.coloc_table, sep="\t")
    coloc = coloc[["snp", "trait", "PP.H4.abf"]].drop_duplicates()
    snp_traits = (coloc.sort_values("PP.H4.abf", ascending=False)
                       .groupby("snp")
                       .apply(lambda d: list(zip(d.trait, d["PP.H4.abf"])))
                       .to_dict())
    snps = sorted(snp_traits)
    if args.snps:
        keep = set(s.strip() for s in args.snps.split(","))
        snps = [s for s in snps if s in keep]
    print(f"{len(snps)} SNPs from coloc table; modalities {modalities}")

    fiber_tracts = load_fiber_tracts(fiber_tract_path())
    null = None if args.no_null else RotationNull(fiber_tracts, args.n_rot, args.seed)

    rows, errors = [], []
    for mod in modalities:
        # disease maps for this modality
        gc_maps = {}
        for d in DISEASES:
            try:
                gc_maps[d] = load_gc_pmap(gc_path(mod, d))
            except Exception as e:
                errors.append(f"GC LOAD FAILED {mod} {d}: {e}")
        if not gc_maps:
            errors.append(f"no GC maps loaded for {mod}; skipping modality")
            continue
        dz = list(gc_maps)
        gc_valid = np.all([np.isfinite(gc_maps[d]) for d in dz], axis=0)
        gc_sig = [gc_maps[d] < P_GC for d in dz]
        gc_fdr = [bh_mask(gc_maps[d]) for d in dz]
        # threshold-free disease side: -log10 P of the r_g map (NaN where invalid)
        gc_nlp = [np.where(gc_valid, -np.log10(np.clip(gc_maps[d], 1e-300, 1)), np.nan) for d in dz]
        print(f"[{mod}] r_g maps: " + ", ".join(
            f"{d}={int(m.sum())}" for d, m in zip(dz, gc_sig)))

        for snp in snps:
            t0 = time.time()
            try:
                fp = load_pmap(voxelgwas_path(mod, snp))
            except Exception as e:
                errors.append(f"FOOTPRINT LOAD FAILED {mod} {snp}: {e}")
                continue
            if fp.shape != (N_FIBERS, N_POINTS):
                errors.append(f"SHAPE {mod} {snp}: {fp.shape}")
                continue

            valid = np.isfinite(fp) & gc_valid
            n_valid = int(valid.sum())
            fp_sig = (fp < P_FOOTPRINT) & valid
            fp_fdr = bh_mask(fp) & valid
            n_fp = int(fp_sig.sum())

            # top-k masks by -log10 P
            order = np.argsort(np.where(valid, fp, np.inf), axis=None)
            topk_masks = {}
            for k in TOPK:
                m = np.zeros(fp.size, dtype=bool)
                m[order[:k]] = True
                topk_masks[k] = m.reshape(fp.shape) & valid

            # observed overlaps
            masks = {"atlas": fp_sig, "fdr": fp_fdr}
            masks.update({f"top{k}": m for k, m in topk_masks.items()})

            # null overlaps (one set of draws per footprint mask, reused across diseases)
            null_obs, null_nlp = {}, None
            if null is not None:
                for name, m in masks.items():
                    if name == "fdr":
                        continue
                    if m.sum() == 0:
                        null_obs[name] = None
                    elif name == "atlas":
                        null_obs[name], null_nlp = null_overlaps(null, m, gc_sig, valid, gc_nlp)
                    else:
                        null_obs[name], _ = null_overlaps(null, m, gc_sig, valid)

            traits = snp_traits[snp]
            trait_codes = {TRAIT_TO_GC.get(t) for t, _ in traits} - {None}

            for j, d in enumerate(dz):
                g = gc_sig[j] & valid
                n_gc = int(g.sum())
                ov = int((fp_sig & g).sum())
                exp = n_fp * n_gc / n_valid if n_valid else np.nan
                fold = ov / exp if exp > 0 else np.nan
                row = {
                    "modality": mod, "snp": snp,
                    "coloc_traits": ";".join(f"{t}:{h4:.3f}" for t, h4 in traits),
                    "disease": d,
                    "primary": d in trait_codes,
                    "n_valid": n_valid, "n_fp_atlas": n_fp, "n_gc": n_gc,
                    "overlap": ov, "expected": exp, "fold": fold,
                    "fold_ceiling": (n_valid / n_fp) if n_fp else np.nan,
                    "frac_fp": n_fp / n_valid if n_valid else np.nan,
                    "frac_gc_inside_fp": ov / n_gc if n_gc else np.nan,
                    "dice": 2 * ov / (n_fp + n_gc) if (n_fp + n_gc) else np.nan,
                }
                # tracts where the overlap lives
                if ov:
                    fib = np.where((fp_sig & g).any(axis=1))[0]
                    tr, cnt = np.unique(fiber_tracts[fib], return_counts=True)
                    row["overlap_tracts"] = ";".join(
                        f"{t}:{c}" for t, c in sorted(zip(tr, cnt), key=lambda x: -x[1]))
                else:
                    row["overlap_tracts"] = ""
                # rotation null
                if null is not None:
                    no = null_obs.get("atlas")
                    if no is None or n_gc == 0:
                        row.update(null_mean=np.nan, null_p95=np.nan, p_rot=np.nan, z_rot=np.nan)
                    else:
                        nv = no[:, j]
                        row["null_mean"] = float(nv.mean())
                        row["null_p95"] = float(np.percentile(nv, 95))
                        row["p_rot"] = (1 + int((nv >= ov).sum())) / (1 + len(nv))
                        sd = nv.std()
                        row["z_rot"] = (ov - nv.mean()) / sd if sd > 0 else np.nan
                # threshold-free disease side: mean -log10 P of the r_g map inside the footprint
                if n_fp:
                    obs_nlp = float(np.nanmean(gc_nlp[j][fp_sig]))
                    row["gc_mean_nlp_fp"] = obs_nlp
                    row["gc_mean_nlp_all"] = float(np.nanmean(gc_nlp[j]))
                    if null_nlp is not None:
                        nv = null_nlp[:, j]
                        nv = nv[np.isfinite(nv)]
                        if nv.size:
                            row["null_mean_nlp"] = float(nv.mean())
                            row["p_rot_nlp"] = (1 + int((nv >= obs_nlp).sum())) / (1 + nv.size)
                # top-k sensitivity
                for k in TOPK:
                    mk = topk_masks[k]
                    ovk = int((mk & g).sum())
                    expk = mk.sum() * n_gc / n_valid if n_valid else np.nan
                    row[f"fold_top{k}"] = ovk / expk if expk > 0 else np.nan
                    if null is not None and null_obs.get(f"top{k}") is not None and n_gc:
                        nv = null_obs[f"top{k}"][:, j]
                        row[f"p_rot_top{k}"] = (1 + int((nv >= ovk).sum())) / (1 + len(nv))
                # sensitivity: BH-FDR thresholds on both sides
                gf = gc_fdr[j] & valid
                n_fpf, n_gcf = int(fp_fdr.sum()), int(gf.sum())
                ovf = int((fp_fdr & gf).sum())
                expf = n_fpf * n_gcf / n_valid if n_valid else np.nan
                row.update(n_fp_fdr=n_fpf, n_gc_fdr=n_gcf, overlap_fdr=ovf,
                           fold_fdr=ovf / expf if expf > 0 else np.nan)
                rows.append(row)
            print(f"  {snp:12s} footprint={n_fp:6d}  ({time.time() - t0:.1f}s)")

    df = pd.DataFrame(rows)
    if df.empty:
        print("no rows produced; see errors")
    else:
        # BH across primary tests on the rotation P
        if "p_rot" in df:
            pri = df["primary"] & df["p_rot"].notna()
            df["q_rot_primary"] = np.nan
            if pri.any():
                p = df.loc[pri, "p_rot"].to_numpy()
                n = p.size
                order = np.argsort(p)
                q = np.empty(n)
                q[order] = np.minimum.accumulate((p[order] * n / np.arange(1, n + 1))[::-1])[::-1]
                df.loc[pri, "q_rot_primary"] = np.minimum(q, 1)
        df = df.sort_values(["primary", "fold"], ascending=[False, False])
        df.to_csv(args.out, index=False, float_format="%.6g")
        print(f"\nwrote {len(df)} rows to {args.out}")
        pd.set_option("display.width", 250)
        pd.set_option("display.max_columns", None)
        show = ["modality", "snp", "disease", "primary", "n_fp_atlas", "n_gc",
                "overlap", "fold", "fold_ceiling"]
        if "p_rot" in df:
            show += ["p_rot", "q_rot_primary"]
        if "p_rot_nlp" in df:
            show += ["gc_mean_nlp_fp", "null_mean_nlp", "p_rot_nlp"]
        print("\nPrimary tests (locus vs own coloc trait):")
        print(df[df.primary][show].to_string(index=False))
        print("\nTop 15 exploratory by fold (n_gc >= 50, n_fp >= 50):")
        ex = df[(~df.primary) & (df.n_gc >= 50) & (df.n_fp_atlas >= 50)]
        print(ex.head(15)[show].to_string(index=False))

    if errors:
        Path(args.out).with_suffix(".errors.txt").write_text("\n".join(errors))
        print(f"\n{len(errors)} errors logged to {Path(args.out).with_suffix('.errors.txt')}")


if __name__ == "__main__":
    main()
