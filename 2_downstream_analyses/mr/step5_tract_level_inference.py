#!/usr/bin/env python3
"""
step5_tract_level_inference.py

Tract-level Mendelian randomization inference from the fiber-cluster (LDR)
level TwoSampleMR results of step 4.

Rule
  1. Keep cluster-level results with at least 5 instruments (nsnp >= 5).
  2. For each tract x outcome x estimator, combine the cluster-level P-values
     with the Cauchy combination test (ACAT, equal weights) and compute the
     -log10(P)-weighted fraction of fits whose sign(b) agrees with the
     fixed-effect IVW estimate of the same cluster x LDR x outcome.
  3. Benjamini-Hochberg FDR of the ACAT P-values across all tract x outcome
     cells, separately for each estimator.
  4. An estimator is "concordant-significant" in a cell if q < 0.05 and the
     weighted sign concordance is 1.
  5. A tract x outcome pair is reported if fixed-effect IVW has q < 0.05 and
     at least 3 of the 12 estimators (fixed-effect IVW included) are
     concordant-significant.

Inputs
  --mr-results     mr_all_results.csv from step4_agg_res.R
  --cluster-tract  cluster_tract_summary.csv (columns cluster, tract_assign)

Outputs (in --out-dir)
  mr_tract_outcome_method.csv   one row per tract x outcome x estimator
  mr_tract_outcome_pairs.csv    one row per tract x outcome, with the flag
  mr_significant_pairs.csv      the pairs passing the rule

Usage:
    python step5_tract_level_inference.py \
        --mr-results    /path/to/mr/agg/mr_all_results.csv \
        --cluster-tract /path/to/cluster_tract_summary.csv \
        --out-dir       /path/to/mr/agg
"""
import argparse
import os

import numpy as np
import pandas as pd

IVW = 'Inverse variance weighted (fixed effects)'
NSNP_MIN, Q_THR, SIGN_THR, N_METHOD_THR = 5, 0.05, 1.0, 3

METHOD_ORDER = [IVW, 'Inverse variance weighted', 'IVW radial', 'Simple median', 'Weighted median', 'Penalised weighted median',
                'Simple mode', 'Simple mode (NOME)', 'Weighted mode', 'Weighted mode (NOME)', 'MR Egger', 'MR GRIP']
METHOD_LABEL = {IVW: 'IVW (fixed effects)',
                'Inverse variance weighted': 'IVW (random effects; TwoSampleMR default)', 'IVW radial': 'IVW (radial)',
                'Simple median': 'Simple median', 'Weighted median': 'Weighted median',
                'Penalised weighted median': 'Penalised weighted median', 'Simple mode': 'Simple mode',
                'Simple mode (NOME)': 'Simple mode (NOME)', 'Weighted mode': 'Weighted mode',
                'Weighted mode (NOME)': 'Weighted mode (NOME)', 'MR Egger': 'MR-Egger', 'MR GRIP': 'MR-GRIP'}
TRACT_ORDER = ['ACR', 'ALIC', 'BCC', 'CGH', 'CP', 'EC', 'FX', 'GCC', 'IFO', 'PCR', 'PLIC', 'PTR', 'RLIC', 'SCC', 'SCR', 'SLF', 'SS']
OUTCOME_ORDER = ['ADHD', 'ALZ2', 'ASD', 'bip', 'CP', 'depressed_affect', 'depression', 'EA', 'insomnia', 'intelligence',
                 'mdd', 'neuroticism', 'SCZ', 'worry']


def cct(p):
    """Cauchy combination test (ACAT) with equal weights, as GRAB::CCT."""
    p = np.asarray(p, float)
    if (p == 0).any():
        return 0.0
    if (p == 1).any():
        return 1.0
    small = p < 1e-16
    stat = np.sum(np.where(small, (1 / np.where(small, p, 1)) / np.pi, np.tan((0.5 - p) * np.pi))) / len(p)
    if stat > 1e15:
        return (1 / stat) / np.pi
    # upper tail of the standard Cauchy; arctan(1/stat) keeps precision in the tail
    return float(np.arctan(1 / stat) / np.pi) if stat > 1 else float(0.5 - np.arctan(stat) / np.pi)


def bh(p):
    """Benjamini-Hochberg, identical to R p.adjust(method='fdr'); NaN left as NaN and not counted."""
    p = np.asarray(p, float); out = np.full(p.shape, np.nan); ok = ~np.isnan(p); pv = p[ok]; n = len(pv)
    o = np.argsort(pv)[::-1]; q = np.minimum(1, np.minimum.accumulate(n / np.arange(n, 0, -1) * pv[o]))
    tmp = np.empty(n); tmp[o] = q; out[ok] = tmp; return out


def main():
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument('--mr-results', required=True)
    ap.add_argument('--cluster-tract', required=True)
    ap.add_argument('--out-dir', required=True)
    args = ap.parse_args()
    os.makedirs(args.out_dir, exist_ok=True)

    # ---------------- load ----------------
    mr = pd.read_csv(args.mr_results)
    assert set(mr.method) == set(METHOD_ORDER), set(mr.method) ^ set(METHOD_ORDER)
    ta = pd.read_csv(args.cluster_tract)[['cluster', 'tract_assign']]
    ta = ta.rename(columns={'cluster': 'fiber_id', 'tract_assign': 'tract'})
    n_raw = len(mr)
    mr = mr[mr.nsnp >= NSNP_MIN].merge(ta, on='fiber_id', how='left')
    assert set(mr.trait) == set(OUTCOME_ORDER)
    ref = mr[mr.method == IVW][['fiber_id', 'exposure', 'outcome', 'b']].rename(columns={'b': 'b_ivw'})
    mr = mr.merge(ref, on=['fiber_id', 'exposure', 'outcome'], how='left')
    mr['agree'] = (np.sign(mr.b) == np.sign(mr.b_ivw)).astype(float)
    mr['w'] = -np.log10(mr.pval.clip(lower=1e-300))
    mr_t = mr.dropna(subset=['tract'])  # clusters not assigned to a JHU tract are excluded

    # ---------------- cell level (tract x outcome x estimator) ----------------
    def summ(x):
        sign = np.sign(x.b.values)
        return pd.Series({
            'n_clusters_contributing': int(x.fiber_id.nunique()),
            'n_cluster_LDR_exposures': len(x),
            'n_clusters_sig_nominal': int((x.pval < 0.05).sum()),
            'acat_p': cct(x.pval.values),
            'mean_sign_concordance_w': np.average(x.agree, weights=x.w) if x.w.sum() > 0 else np.nan,
            'frac_positive_b': float((sign > 0).mean()),
            'median_b': float(np.median(x.b)),
            'median_nsnp': float(np.median(x.nsnp)),
            'total_nsnp': int(x.nsnp.sum()),
        })

    cells = mr_t.groupby(['tract', 'trait', 'method']).apply(summ, include_groups=False).reset_index()
    cells['q_grid'] = cells.groupby('method').acat_p.transform(bh)
    n_cells = cells.groupby('method').size().unique()
    assert len(n_cells) == 1
    n_grid = int(n_cells[0])
    cells['fdr_significant'] = cells.q_grid < Q_THR
    cells['concordant_significant'] = cells.fdr_significant & (cells.mean_sign_concordance_w >= SIGN_THR)
    cells['direction'] = np.where(cells.frac_positive_b > 0.5, '+', np.where(cells.frac_positive_b < 0.5, '-', '0'))
    cells['direction_consistent'] = cells.mean_sign_concordance_w >= SIGN_THR

    # full tract x outcome x estimator scaffold (tracts without MR results are listed with empty values)
    full = pd.MultiIndex.from_product([TRACT_ORDER, OUTCOME_ORDER, METHOD_ORDER], names=['tract', 'trait', 'method']).to_frame(index=False)
    cells = full.merge(cells, on=['tract', 'trait', 'method'], how='left')

    # ---------------- pair level ----------------
    def pair(x):
        ivw = x[x.method == IVW].iloc[0]
        return pd.Series({
            'ivw_fixed_acat_p': ivw.acat_p, 'ivw_fixed_q': ivw.q_grid,
            'ivw_fixed_direction': ivw.direction, 'ivw_fixed_median_b': ivw.median_b,
            'n_methods_fdr_sig': int(x.fdr_significant.sum()) if x.acat_p.notna().any() else np.nan,
            'n_methods_concordant_sig': int(x.concordant_significant.sum()) if x.acat_p.notna().any() else np.nan,
            'concordant_methods': '; '.join(METHOD_LABEL[m] for m in METHOD_ORDER
                                            if x.loc[x.method == m, 'concordant_significant'].fillna(False).iloc[0]),
            'n_clusters_contributing': ivw.n_clusters_contributing,
            'n_cluster_LDR_exposures': ivw.n_cluster_LDR_exposures,
            'instruments_median_per_cluster': ivw.median_nsnp,
            'instruments_total_across_clusters': ivw.total_nsnp,
        })

    pairs = cells.groupby(['tract', 'trait']).apply(pair, include_groups=False).reset_index()
    pairs['significant_flag'] = ((pairs.ivw_fixed_q < Q_THR) & (pairs.n_methods_concordant_sig >= N_METHOD_THR)).fillna(False)
    pairs = pairs.sort_values(['significant_flag', 'ivw_fixed_q'], ascending=[False, True])
    sig = pairs[pairs.significant_flag]

    cells_out = cells.copy()
    cells_out['method'] = cells_out.method.map(METHOD_LABEL)
    cells_out.to_csv(os.path.join(args.out_dir, 'mr_tract_outcome_method.csv'), index=False)
    pairs.to_csv(os.path.join(args.out_dir, 'mr_tract_outcome_pairs.csv'), index=False)
    sig.drop(columns=['significant_flag']).to_csv(os.path.join(args.out_dir, 'mr_significant_pairs.csv'), index=False)

    print(f'{n_raw:,} fits read; {len(mr):,} with nsnp >= {NSNP_MIN}')
    print(f'grid: {n_grid} tract x outcome cells; significant pairs: {len(sig)}')
    print(sig.groupby('trait').size().reindex(OUTCOME_ORDER).fillna(0).astype(int).to_string())


if __name__ == '__main__':
    main()
