#!/usr/bin/env python3
"""
step1_drug_target_lookup.py

Druggability annotation of the genes at the white-matter GWAS loci with the
Open Targets Platform.

For each gene the script
  - resolves the Ensembl ID,
  - retrieves tractability buckets (small molecule, antibody, other modalities)
    and the drugs / clinical candidates recorded for the target, with mechanism
    of action and indications,
  - assigns one of four tiers: approved-drug target, clinical-stage target,
    tractable (no compound yet), not currently tractable,
  - optionally joins the colocalized traits and the enriched tracts of the
    same lead variant (--coloc-file, --tracts-file).

The tier counts are descriptive: no matched-background enrichment test is run.

Usage
-----
  pip install requests pandas matplotlib

  # connectivity and output-format check with a small built-in gene list
  python step1_step1_drug_target_lookup.py --demo

  # gene list
  python step1_step1_drug_target_lookup.py --genes-file my_genes.csv --out results

`my_genes.csv` is either a single column of gene symbols (header optional) or
a table with a column named gene, gene_symbol, effector_gene or symbol, and
optionally a `locus` or `network_pair` column that is carried to the output.
A g:Profiler / FUMA gene table with Ensembl IDs can be given with
--gprofiler-file (see --help).

Outputs (under --out, default "drug_target_results")
  gene_drug_results.csv            one row per gene x drug hit (or one row per
                                   gene without a drug)
  tier_summary.csv                 counts and percentages per tier
  druggable_disease_tract_table.csv  approved / clinical-stage targets with
                                   colocalized traits and tracts (if requested)
  fig_tractability_bar.png, fig_approved_hits_table.png

The script queries the public Open Targets GraphQL API
(https://api.platform.opentargets.org/api/v4/graphql; no key needed) and
therefore needs internet access. Results in the paper were obtained with
Open Targets Platform release 26.06.
"""
from __future__ import annotations

import argparse
import json
import os
import re
import sys
import time

import pandas as pd
import requests

OT_API = "https://api.platform.opentargets.org/api/v4/graphql"

SEARCH_QUERY = """
query SearchGene($q: String!) {
  search(queryString: $q, entityNames: ["target"], page: {index: 0, size: 5}) {
    hits {
      id
      name
      entity
    }
  }
}
"""

TARGET_QUERY = """
query TargetInfo($ensemblId: String!) {
  target(ensemblId: $ensemblId) {
    id
    approvedSymbol
    tractability {
      label
      modality
      value
    }
    drugAndClinicalCandidates {
      count
      rows {
        maxClinicalStage
        drug {
          id
          name
          drugType
          mechanismsOfAction {
            rows {
              mechanismOfAction
            }
          }
        }
        diseases {
          disease {
            id
            name
          }
        }
      }
    }
  }
}
"""
# NOTE: Open Targets Platform release 26.06 renamed `Target.knownDrugs` to
# `Target.drugAndClinicalCandidates` and reshaped it: `phase`, `status`,
# `mechanismOfAction` (string), `drug.isApproved` and `disease` were replaced by
# `maxClinicalStage`, `drug.mechanismsOfAction.rows[].mechanismOfAction` and
# `diseases` (a list). The queries below follow the 26.06 schema. SEARCH_QUERY
# (gene-symbol lookup, used for --demo / --genes-file) may need the same kind
# of update on later releases; _post_graphql() prints the API's
# "Cannot query field ..." message when a field has been renamed.

DEMO_GENES = ["APOE", "PCSK9", "NR3C1", "GRIN2A", "MAPT", "HTR2A", "DRD2", "SORT1"]


def _post_graphql(query, variables, *, max_retries=3, timeout=20):
    """POSTs a GraphQL query. IMPORTANT: many GraphQL servers (including Open
    Targets') return HTTP 400 -- not 200 with an "errors" array -- for a
    query-validation failure (e.g. an unknown field), per the GraphQL-over-
    HTTP spec. So we must read and surface the JSON body on a 400 BEFORE
    treating it as a generic HTTP error, or the real diagnostic ("Cannot
    query field X on type Y") gets thrown away and replaced with an opaque
    "400 Client Error: Bad Request"."""
    last_err = None
    for attempt in range(max_retries):
        try:
            resp = requests.post(
                OT_API, json={"query": query, "variables": variables}, timeout=timeout
            )
            # try to parse JSON regardless of status code -- a 400 from a
            # GraphQL validation error still has a useful JSON body
            try:
                payload = resp.json()
            except ValueError:
                payload = None

            if payload is not None and "errors" in payload:
                raise RuntimeError(
                    f"GraphQL error (HTTP {resp.status_code}) for {variables}:\n"
                    f"{json.dumps(payload['errors'], indent=2)[:1500]}\n"
                    "This is almost always an unknown/renamed field -- the message above "
                    "names it directly. Fix it in SEARCH_QUERY/TARGET_QUERY and rerun."
                )
            resp.raise_for_status()
            return payload["data"]
        except (requests.RequestException, RuntimeError) as e:
            last_err = e
            time.sleep(1.5 * (attempt + 1))
    raise last_err


def resolve_ensembl_id(gene_symbol):
    data = _post_graphql(SEARCH_QUERY, {"q": gene_symbol})
    hits = data.get("search", {}).get("hits", [])
    for h in hits:
        if h.get("entity") == "target" and h.get("name", "").upper() == gene_symbol.upper():
            return h["id"]
    # fall back to first target hit if no exact symbol match
    for h in hits:
        if h.get("entity") == "target":
            return h["id"]
    return None


def clinical_stage(row):
    """Numeric coercion of a drugAndClinicalCandidates row's maxClinicalStage,
    used only for the `phase` column of the output and not on its own for the
    tier classification (see classify_tier). The field has been a numeric
    string ("4"), a free-text label ("Approved", "Phase III") or missing,
    depending on the Open Targets API version, so numbers, numeric strings and
    recognizable labels are all parsed and anything else returns None."""
    v = row.get("maxClinicalStage")
    if v is None:
        return None
    if isinstance(v, (int, float)):
        return float(v)
    s = str(v).strip()
    if not s:
        return None
    try:
        return float(s)
    except ValueError:
        pass
    s_low = s.lower()
    if "approv" in s_low:
        return 4.0
    m = re.search(r"phase\s*(iv|iii|ii|i|4|3|2|1)\b", s_low)
    if m:
        roman = {"iv": 4, "iii": 3, "ii": 2, "i": 1, "4": 4, "3": 3, "2": 2, "1": 1}
        return float(roman[m.group(1)])
    if "preclinical" in s_low:
        return 0.0
    return None


APPROVED_TRACT_LABELS = {"Approved Drug"}
CLINICAL_TRACT_LABELS = {"Advanced Clinical", "Phase 1 Clinical"}


def classify_tier(tractability_rows, candidate_rows):
    """PRIMARY signal: the target's own `tractability` buckets -- these flag
    'Approved Drug' / 'Advanced Clinical' / 'Phase 1 Clinical' as simple
    stable booleans per modality, independent of whatever shape
    drugAndClinicalCandidates.rows.maxClinicalStage happens to be in this
    API version (see clinical_stage() docstring -- that field's format is
    NOT reliable and should not be trusted alone; confirmed live: DPP4 and
    EPHA3 both carry known approved drugs by name -- Sitagliptin,
    Vildagliptin, Alogliptin, etc. for DPP4 -- but maxClinicalStage came back
    null/unparseable for every single row of both). clinical_stage() is kept
    only as a secondary/fallback signal in case a future API version makes
    it reliable again."""
    tract_labels = {t.get("label") for t in tractability_rows if t.get("value")}
    stages = [clinical_stage(r) for r in candidate_rows]

    if (tract_labels & APPROVED_TRACT_LABELS) or any(s == 4 for s in stages):
        return "approved-drug target"
    if (tract_labels & CLINICAL_TRACT_LABELS) or any((s or 0) >= 1 for s in stages):
        return "clinical-stage target"
    is_tractable = any(
        (t.get("modality") in ("SM", "AB", "PR", "OC")) and t.get("value") for t in tractability_rows
    )
    if is_tractable:
        return "tractable (no compound yet)"
    return "not currently tractable"


def tractability_evidence(tractability_rows):
    """Turn the tractability bucket list into a short 'modality:evidence' summary
    string, e.g. 'SM:High-Quality Ligand; AB:GO CC high conf'. This is the extra
    detail behind a 'tractable (no compound yet)' tier -- there's no known drug,
    but this says WHY Open Targets thinks the target class is druggable at all
    (crystal structure with a ligand, a well-characterized binding pocket,
    favorable subcellular location for an antibody, etc.), which is exactly the
    kind of evidence worth reporting for a gene with no compound yet."""
    hits = [t for t in tractability_rows if t.get("value")]
    if not hits:
        return None
    return "; ".join(f"{t.get('modality')}:{t.get('label')}" for t in hits)


def run_lookup(genes, loci_map=None, ensembl_map=None, modality_map=None, *, sleep_s=0.3, verbose=True):
    """genes: list of gene symbols. loci_map: optional {gene: locus_label}.
    ensembl_map: optional {gene: ensembl_id} -- when given, skips the symbol
    search step entirely (faster AND more accurate: no ambiguity from
    Open Targets' fuzzy symbol search). modality_map: optional
    {gene: "AD;FA;MD"-style label} carried straight through to output."""
    rows = []
    loci_map = loci_map or {}
    ensembl_map = ensembl_map or {}
    modality_map = modality_map or {}
    for i, gene in enumerate(genes):
        if verbose:
            print(f"[{i+1}/{len(genes)}] {gene} ...", end=" ", flush=True)
        locus = loci_map.get(gene)
        modality = modality_map.get(gene)
        try:
            ensembl_id = ensembl_map.get(gene) or resolve_ensembl_id(gene)
            if ensembl_id is None:
                if verbose:
                    print("NOT FOUND in Open Targets search")
                rows.append({"gene": gene, "locus": locus, "modality": modality, "ensembl_id": None,
                             "tier": "not found", "drug": None, "phase": None,
                             "mechanism": None, "indication": None, "tractability_evidence": None})
                continue

            data = _post_graphql(TARGET_QUERY, {"ensemblId": ensembl_id})
            target = data.get("target") or {}
            tractability = target.get("tractability") or []
            candidates = (target.get("drugAndClinicalCandidates") or {}).get("rows") or []
            tier = classify_tier(tractability, candidates)
            evidence = tractability_evidence(tractability)

            if not candidates:
                rows.append({"gene": gene, "locus": locus, "modality": modality, "ensembl_id": ensembl_id,
                             "tier": tier, "drug": None, "phase": None,
                             "mechanism": None, "indication": None, "tractability_evidence": evidence})
            else:
                # one row per known drug/indication pair, so the same gene can
                # appear multiple times if it has several drugs -- collapse
                # with tier_summary.csv / drop_duplicates(subset='gene') downstream
                # when only gene-level counts are needed.
                for r in candidates:
                    drug = r.get("drug") or {}
                    moa_rows = (drug.get("mechanismsOfAction") or {}).get("rows") or []
                    mechanism = moa_rows[0].get("mechanismOfAction") if moa_rows else None
                    diseases = r.get("diseases") or []
                    indication = (diseases[0].get("disease") or {}).get("name") if diseases else None
                    rows.append({
                        "gene": gene, "locus": locus, "modality": modality, "ensembl_id": ensembl_id,
                        "tier": tier,
                        "drug": drug.get("name"),
                        "phase": clinical_stage(r),
                        "mechanism": mechanism,
                        "indication": indication,
                        "tractability_evidence": evidence,
                    })
            if verbose:
                print(tier)
        except Exception as e:
            if verbose:
                print(f"ERROR: {e}")
            rows.append({"gene": gene, "locus": locus, "modality": modality, "ensembl_id": None,
                         "tier": "error", "drug": None, "phase": None,
                         "mechanism": None, "indication": str(e)[:200], "tractability_evidence": None})
        time.sleep(sleep_s)
    return pd.DataFrame(rows)


# ---------------------------------------------------------------------------
# g:Profiler SNP->gene export parsing (the format this actually got used
# with: a blank first row, then id/chr/start/end/strand/gene_ids/gene_names/
# <variant-consequence columns>, one row per lead SNP, gene_ids and
# gene_names comma-separated when a SNP maps to >1 gene, and some SNPs
# entirely unmapped -- chr/start/end == -1, gene_ids blank).
# ---------------------------------------------------------------------------

def load_gprofiler_genelist(path, modality):
    """One g:Profiler export -> long df with one row per (SNP, gene)."""
    df = pd.read_csv(path, skiprows=1)
    df = df.dropna(subset=["gene_ids"])
    rows = []
    for _, r in df.iterrows():
        locus = r["id"]
        gene_ids = [g.strip() for g in str(r["gene_ids"]).split(",")]
        gene_names = [g.strip() for g in str(r["gene_names"]).split(",")]
        for i, gid in enumerate(gene_ids):
            # g:Profiler sometimes can't resolve a symbol and just repeats
            # the Ensembl ID in gene_names -- keep it (still queryable by
            # ensembl_id), but it won't look like a "real" gene symbol.
            gname = gene_names[i] if i < len(gene_names) else gid
            rows.append({"locus": locus, "ensembl_id": gid, "gene_symbol": gname, "modality": modality})
    return pd.DataFrame(rows)


def combine_gprofiler_genelists(paths_by_modality):
    """{modality: path} -> one row per unique Ensembl gene, with all loci and
    modalities it was implicated by joined (";"-separated) for traceability.
    This is the UNION across modalities -- a gene counts once even if it shows
    up under more than one modality's lead-SNP list."""
    dfs = [load_gprofiler_genelist(p, m) for m, p in paths_by_modality.items()]
    all_df = pd.concat(dfs, ignore_index=True)

    # prefer a real symbol over a row where g:Profiler couldn't resolve one
    # (and just repeated the Ensembl ID) when collapsing duplicates below
    all_df["is_real_symbol"] = all_df["gene_symbol"] != all_df["ensembl_id"]
    all_df = all_df.sort_values("is_real_symbol", ascending=False)

    grouped = all_df.groupby("ensembl_id").agg(
        gene_symbol=("gene_symbol", "first"),
        loci=("locus", lambda x: ";".join(sorted(set(x)))),
        modalities=("modality", lambda x: ";".join(sorted(set(x)))),
    ).reset_index()
    return grouped


def load_gene_list(path):
    df = pd.read_csv(path, sep=None, engine="python")
    gene_col = None
    for cand in ("gene", "gene_symbol", "effector_gene", "symbol"):
        if cand in df.columns:
            gene_col = cand
            break
    if gene_col is None:
        # single-column file, possibly no header (pandas may have consumed
        # the first gene as a header) -- reread with no header as a fallback
        if df.shape[1] == 1:
            genes = pd.read_csv(path, header=None)[0].astype(str).tolist()
            return genes, {}
        raise ValueError(
            f"Could not find a gene column in {path} -- columns are {list(df.columns)}. "
            "Expected one named gene/gene_symbol/effector_gene/symbol, or a single-column file."
        )
    genes = df[gene_col].dropna().astype(str).str.strip().tolist()
    loci_map = {}
    locus_col = next((c for c in ("locus", "network_pair") if c in df.columns), None)
    if locus_col:
        loci_map = dict(zip(df[gene_col].astype(str).str.strip(), df[locus_col]))
    # dedupe, preserve order
    seen = set()
    uniq = []
    for g in genes:
        if g not in seen:
            seen.add(g)
            uniq.append(g)
    return uniq, loci_map


def summarize(df):
    per_gene = df.drop_duplicates(subset="gene")[["gene", "tier"]]
    counts = per_gene["tier"].value_counts()
    order = ["approved-drug target", "clinical-stage target", "tractable (no compound yet)",
             "not currently tractable", "not found", "error"]
    counts = counts.reindex([o for o in order if o in counts.index]).fillna(0).astype(int)
    total = len(per_gene)
    summary = pd.DataFrame({
        "tier": counts.index,
        "n_genes": counts.values,
        "pct_of_total": (counts.values / total * 100).round(1) if total else 0,
    })
    return summary, total


# ---------------------------------------------------------------------------
# Disease-colocalization / tract-enrichment join -- "for the druggable genes,
# which diseases are they related to, and which tracts are enriched" (the
# drug-target evidence and the disease/tract evidence both come from the
# SAME lead SNP, so joining on `locus` == the g:Profiler `id`/lead-SNP column
# is exactly the right key -- no separate coloc/tract lookup needed).
# ---------------------------------------------------------------------------

def load_coloc_prioritized(path):
    """A *_coloc_prioritized.tsv (tab-separated): one row per locus with the
    disease/trait it colocalized with. Key columns: snp (lead SNP), trait,
    trait_type (quant/cc), PP.H4.abf (colocalization posterior probability --
    the standard coloc.abf evidence strength for "same causal variant"),
    is_novel."""
    df = pd.read_csv(path, sep="\t")
    keep = [c for c in ("snp", "trait", "trait_type", "PP.H4.abf", "is_novel", "locus_id") if c in df.columns]
    out = df[keep].rename(columns={"PP.H4.abf": "coloc_pp_h4"})
    return out


def load_tract_enrichment(path):
    """A *_all.csv (comma-separated, quoted, first column an unnamed row
    index): one row per locus with the white-matter tracts it showed signal
    in. Key columns: snp (lead SNP), tracts (comma-separated tract
    abbreviations), n_tracts."""
    df = pd.read_csv(path)
    keep = [c for c in ("snp", "tracts", "n_tracts") if c in df.columns]
    return df[keep]


def build_disease_tract_table(gene_drug_df, coloc_df=None, tracts_df=None,
                               tiers=("approved-drug target", "clinical-stage target")):
    """One row per (gene, locus) among the druggable tiers (approved/clinical
    by default -- pass tiers=(..., "tractable (no compound yet)") to widen
    it), collapsing that gene's possibly-many drug rows into a short summary,
    then left-joining in the colocalized disease/trait and the enriched
    tracts for that SAME lead SNP. Genes with no coloc/tract match still
    appear (the join is left, not inner) with those columns blank -- that's
    expected if a gene's locus wasn't in the colocalization table (not
    every locus colocalizes with something) or wasn't in scope for *_all.csv."""
    d = gene_drug_df[gene_drug_df["tier"].isin(tiers)].copy()
    if d.empty:
        return pd.DataFrame()
    rows = []
    for (gene, locus), g in d.groupby(["gene", "locus"], sort=False):
        drugs = sorted(set(x for x in g["drug"].dropna()))
        mechs = sorted(set(x for x in g["mechanism"].dropna()))
        inds = sorted(set(x for x in g["indication"].dropna()))
        example_drugs = ", ".join(drugs[:3]) + (f" (+{len(drugs) - 3} more)" if len(drugs) > 3 else "")
        rows.append({
            "gene": gene, "locus": locus, "tier": g["tier"].iloc[0],
            "n_drugs": len(drugs),
            "example_drugs": example_drugs if drugs else None,
            "mechanism": "; ".join(mechs) if mechs else None,
            "indications": "; ".join(inds) if inds else None,
        })
    summary = pd.DataFrame(rows)
    if coloc_df is not None and len(coloc_df):
        summary = summary.merge(coloc_df, left_on="locus", right_on="snp", how="left", suffixes=("", "_coloc"))
    if tracts_df is not None and len(tracts_df):
        summary = summary.merge(tracts_df, left_on="locus", right_on="snp", how="left", suffixes=("", "_tract"))
    keep = ["gene", "locus", "tier", "n_drugs", "example_drugs", "mechanism", "indications",
            "trait", "trait_type", "coloc_pp_h4", "is_novel", "tracts", "n_tracts"]
    return summary[[c for c in keep if c in summary.columns]]


def make_disease_tract_table_figure(table, out_path, title="Druggable genes: colocalized disease/trait + enriched WM tracts"):
    """Render build_disease_tract_table()'s output as a presentation-ready
    image table (matches fig_approved_hits_table.png's styling). matplotlib's
    table `wrap=True` does NOT reliably wrap long cell text (drug lists, long
    tract-abbreviation lists) -- it just clips at the column edge. So we
    pre-wrap with textwrap ourselves (inserting literal newlines, which
    matplotlib DOES render correctly) and size each row's height to its own
    tallest cell, rather than relying on the auto-wrap flag."""
    import textwrap
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    if table.empty:
        print("[note] no druggable-tier genes to put in the disease/tract table -- skipping figure")
        return

    def _fmt_pp(v):
        return f"{v:.2f}" if pd.notna(v) else "—"

    # (header, wrap-width-in-chars, relative column-width-weight)
    col_spec = [
        ("gene", "Gene", 10, 0.7),
        ("locus", "Lead SNP", 12, 0.9),
        ("example_drugs", "Drug(s)", 28, 1.8),
        ("trait", "Colocalized trait", 16, 1.0),
        ("coloc_pp_h4", "Coloc PP.H4", 10, 0.7),
        ("tracts", "Enriched tracts", 28, 1.8),
    ]
    col_spec = [c for c in col_spec if c[0] in table.columns]

    t = table.copy()
    if "coloc_pp_h4" in t.columns:
        t["coloc_pp_h4"] = t["coloc_pp_h4"].apply(_fmt_pp)

    def _wrap(val, width):
        s = "—" if pd.isna(val) else str(val)
        # tracts/drug lists are comma-joined with no space (e.g. "ACR,EC,FX,...")
        # -- textwrap has no whitespace to break on there, so it falls back to
        # slicing mid-word (e.g. "SCR" -> "S" / "CR"). Give it a break point.
        s = s.replace(",", ", ").replace("  ", " ")
        return "\n".join(textwrap.wrap(s, width=width, break_long_words=False, break_on_hyphens=False)) or s

    wrapped_rows = []
    for _, r in t.iterrows():
        wrapped_rows.append([_wrap(r[col], width) for col, _, width, _ in col_spec])

    n_rows = len(wrapped_rows)
    # row height driven by the tallest (most-wrapped) cell in that row
    row_lines = [max(cell.count("\n") + 1 for cell in row) for row in wrapped_rows]
    fig_h = 1.3 + sum(0.28 * n for n in row_lines) + 0.35  # + header row
    fig, ax = plt.subplots(figsize=(14, fig_h), dpi=150)
    ax.axis("off")
    tbl = ax.table(cellText=wrapped_rows, colLabels=[h for _, h, _, _ in col_spec],
                    cellLoc="left", colLoc="left", loc="center")
    tbl.auto_set_font_size(False)
    tbl.set_fontsize(9)

    total_w = sum(w for *_, w in col_spec)
    for col_idx, (_, _, _, weight) in enumerate(col_spec):
        for row_idx in range(n_rows + 1):  # +1 for header
            tbl[row_idx, col_idx].set_width(weight / total_w)

    for row_idx in range(1, n_rows + 1):
        n_lines = row_lines[row_idx - 1]
        for col_idx in range(len(col_spec)):
            tbl[row_idx, col_idx].set_height(0.10 * n_lines + 0.05)
    for col_idx in range(len(col_spec)):
        tbl[0, col_idx].set_height(0.14)

    for (row, col), cell in tbl.get_celld().items():
        cell.set_edgecolor("#DDDDDD")
        cell.PAD = 0.02
        if row == 0:
            cell.set_text_props(weight="bold", color="white", va="center")
            cell.set_facecolor("#3B7DD8")
        else:
            cell.set_text_props(va="center")
            cell.set_facecolor("#F7F9FB" if row % 2 == 0 else "white")
    ax.set_title(title, fontsize=12, loc="left", pad=12)
    fig.tight_layout()
    fig.savefig(out_path, bbox_inches="tight")
    plt.close(fig)


def make_plots(df, summary, out_dir):
    import matplotlib
    matplotlib.use("Agg")
    import matplotlib.pyplot as plt

    # Palette: ordered categorical, warm-to-cool by "how actionable" the tier
    # is (approved = most actionable). Kept to 4 fixed colors, fixed order,
    # never cycled -- consistent with the dataviz skill's categorical rule.
    tier_colors = {
        "approved-drug target": "#3B7DD8",
        "clinical-stage target": "#5FA8D3",
        "tractable (no compound yet)": "#9FB8C8",
        "not currently tractable": "#C9C9C9",
        "not found": "#E8E8E8",
        "error": "#E8E8E8",
    }

    # --- Chart 1: tractability tier bar chart -------------------------------
    plot_summary = summary[~summary["tier"].isin(["not found", "error"])]
    fig, ax = plt.subplots(figsize=(7.5, 4.2), dpi=150)
    colors = [tier_colors.get(t, "#999999") for t in plot_summary["tier"]]
    bars = ax.barh(plot_summary["tier"], plot_summary["n_genes"], color=colors, height=0.6)
    for b, n, pct in zip(bars, plot_summary["n_genes"], plot_summary["pct_of_total"]):
        ax.text(b.get_width() + max(plot_summary["n_genes"]) * 0.02, b.get_y() + b.get_height() / 2,
                f"{n} ({pct}%)", va="center", fontsize=10, color="#333333")
    ax.set_xlabel("Number of genes")
    ax.set_title("Druggability tiers at the white-matter GWAS loci", fontsize=12, loc="left")
    ax.spines[["top", "right"]].set_visible(False)
    ax.invert_yaxis()
    fig.tight_layout()
    fig.savefig(os.path.join(out_dir, "fig_tractability_bar.png"))
    plt.close(fig)

    # --- Chart 2: approved-drug hit table -----------------------------------
    approved = df[df["tier"] == "approved-drug target"].dropna(subset=["drug"])
    approved = approved.drop_duplicates(subset=["gene", "drug", "indication"])
    if len(approved):
        approved = approved.sort_values(["gene", "drug"])
        cell_text = approved[["gene", "drug", "indication"]].fillna("—").values.tolist()
        n_rows = len(cell_text)
        fig_h = max(2.0, 0.35 * n_rows + 1.0)
        fig, ax = plt.subplots(figsize=(9, fig_h), dpi=150)
        ax.axis("off")
        table = ax.table(
            cellText=cell_text,
            colLabels=["Gene", "Approved drug", "Indication"],
            cellLoc="left", colLoc="left", loc="center",
        )
        table.auto_set_font_size(False)
        table.set_fontsize(10)
        table.scale(1, 1.4)
        for (row, col), cell in table.get_celld().items():
            cell.set_edgecolor("#DDDDDD")
            if row == 0:
                cell.set_text_props(weight="bold", color="white")
                cell.set_facecolor("#3B7DD8")
            else:
                cell.set_facecolor("#F7F9FB" if row % 2 == 0 else "white")
        ax.set_title("Genes at the loci already targeted by an approved drug", fontsize=12, loc="left", pad=12)
        fig.tight_layout()
        fig.savefig(os.path.join(out_dir, "fig_approved_hits_table.png"), bbox_inches="tight")
        plt.close(fig)
    else:
        print("[note] no approved-drug hits found -- skipping fig_approved_hits_table.png "
              "(this alone is a valid, presentable preliminary result: 0/N at the approved tier)")


def main(argv=None):
    ap = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    ap.add_argument("--genes-file", help="CSV/TSV with a gene column (see module docstring)")
    ap.add_argument("--gprofiler-file", action="append", default=[],
                     metavar="MODALITY=PATH",
                     help="A g:Profiler SNP->gene export (blank first row, id/chr/.../gene_ids/"
                          "gene_names columns), tagged with its modality, e.g. "
                          "--gprofiler-file AD=AD_genelist.csv. Repeatable -- pass one per "
                          "modality; genes are unioned across all of them, with loci/modalities "
                          "carried through for traceability. Skips the symbol-search step "
                          "entirely since gene_ids already has the Ensembl ID.")
    ap.add_argument("--demo", action="store_true", help=f"Run on a small built-in demo gene list: {DEMO_GENES}")
    ap.add_argument("--coloc-file", help="A *_coloc_prioritized.tsv (tab-separated: locus_id, snp, trait, "
                                          "trait_type, PP.H4.abf, is_novel, ...). If given, joins the "
                                          "colocalized disease/trait onto each druggable gene's lead SNP.")
    ap.add_argument("--tracts-file", help="A *_all.csv (comma-separated: locus_id, snp, tracts, n_tracts, ...). "
                                           "If given, joins the enriched white-matter tracts onto each "
                                           "druggable gene's lead SNP.")
    ap.add_argument("--out", default="drug_target_results")
    ap.add_argument("--sleep", type=float, default=0.3, help="Seconds between API calls (politeness delay)")
    args = ap.parse_args(argv)

    if not args.demo and not args.genes_file and not args.gprofiler_file:
        ap.error("pass --genes-file <path>, --gprofiler-file MODALITY=PATH (repeatable), or --demo")

    os.makedirs(args.out, exist_ok=True)

    ensembl_map, modality_map = {}, {}
    if args.demo:
        genes, loci_map = DEMO_GENES, {}
        print(f"[demo mode] using built-in gene list: {genes}")
    elif args.gprofiler_file:
        paths_by_modality = {}
        for spec in args.gprofiler_file:
            modality, _, path = spec.partition("=")
            if not path:
                ap.error(f"--gprofiler-file expects MODALITY=PATH, got: {spec!r}")
            paths_by_modality[modality] = path
        combined = combine_gprofiler_genelists(paths_by_modality)
        print(f"[gprofiler] {len(combined)} unique gene(s) across {list(paths_by_modality)} "
              f"(union of {sum(len(pd.read_csv(p, skiprows=1)) for p in paths_by_modality.values())} "
              f"lead-SNP rows)")
        genes = combined["gene_symbol"].tolist()
        loci_map = dict(zip(combined["gene_symbol"], combined["loci"]))
        ensembl_map = dict(zip(combined["gene_symbol"], combined["ensembl_id"]))
        modality_map = dict(zip(combined["gene_symbol"], combined["modalities"]))
        combined.to_csv(os.path.join(args.out, "combined_gene_list.csv"), index=False)
    else:
        genes, loci_map = load_gene_list(args.genes_file)
        print(f"[genes] {len(genes)} unique gene(s) loaded from {args.genes_file}")

    df = run_lookup(genes, loci_map, ensembl_map, modality_map, sleep_s=args.sleep)
    df.to_csv(os.path.join(args.out, "gene_drug_results.csv"), index=False)

    summary, total = summarize(df)
    summary.to_csv(os.path.join(args.out, "tier_summary.csv"), index=False)

    print("\n=== Summary (mirrors the paper's '~10% of effector genes target an "
          "approved drug' framing) ===")
    print(summary.to_string(index=False))
    n_approved = summary.loc[summary["tier"] == "approved-drug target", "n_genes"]
    if len(n_approved) and total:
        n = int(n_approved.iloc[0])
        print(f"\n{n} of {total} genes ({n/total*100:.1f}%) are targeted by an approved drug "
              f"(preliminary, unadjusted -- no matched-background enrichment test yet).")

    make_plots(df, summary, args.out)

    if args.coloc_file or args.tracts_file:
        coloc_df = load_coloc_prioritized(args.coloc_file) if args.coloc_file else None
        tracts_df = load_tract_enrichment(args.tracts_file) if args.tracts_file else None
        table = build_disease_tract_table(df, coloc_df, tracts_df)
        if table.empty:
            print("\n[note] no approved/clinical-stage genes to build the disease/tract table from -- "
                  "widen tiers= in build_disease_tract_table() to include 'tractable' genes too.")
        else:
            table.to_csv(os.path.join(args.out, "druggable_disease_tract_table.csv"), index=False)
            make_disease_tract_table_figure(
                table, os.path.join(args.out, "fig_druggable_disease_tract_table.png"))
            print(f"\n=== Druggable genes x colocalized disease/trait x enriched tracts "
                  f"({len(table)} gene-locus row(s)) ===")
            print(table.to_string(index=False))
            print(f"\nWrote: {args.out}/druggable_disease_tract_table.csv, "
                  f"{args.out}/fig_druggable_disease_tract_table.png")

    print(f"\nWrote: {args.out}/gene_drug_results.csv, {args.out}/tier_summary.csv, "
          f"{args.out}/fig_tractability_bar.png"
          + (f", {args.out}/fig_approved_hits_table.png" if os.path.exists(
                os.path.join(args.out, "fig_approved_hits_table.png")) else ""))


if __name__ == "__main__":
    main()
