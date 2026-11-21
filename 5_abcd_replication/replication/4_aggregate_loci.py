#!/usr/bin/env python3
"""
Step 4: merge the per-tract ABCD replication files into loci and write every ABCD number the
manuscript reports. Same procedure as 4_post_screening_replication/step11_aggregate_loci.R, with the
comparison against UK Biobank Phase 7.

Definitions (as in the manuscript)
  locus                 tract-level discovery leads within MERGE_BP (500 kb) of each other, chained
                        along the chromosome, are one locus; the lead is the most significant of them
  locus-level replication  the locus is replicated in >= 1 of its tracts (step 2)
  overlap fraction      sum over tracts of overlapping locations / sum over tracts of discovery locations
  spatially replicated  overlap fraction > OVERLAP_CUT (0.5), strict.  Counts with >= are also printed,
                        because the earlier tables mixed the two (157 used >, "77 in both" used >=).

Outputs (SUMMARY_DIR)
  ALL_TRACTS_tract_locus_replication_merged_<MOD>_<TAG>_ABCD.csv   one row per tract x locus
  ALL_TRACTS_locus_summary_500kb_merged_<MOD>_<TAG>_ABCD.csv       one row per locus (same columns as before)
  abcd_replication_summary.txt       the counts for the text, Table S5-S7 summary row and Supp. Fig. 24 caption
  table_S5_ABCD_columns.csv          the six ABCD_* columns of Supplementary Table 5, by locus
  loci_fa_phenogram.txt              PhenoGram input for Supp. Fig. 24 (CHR POS PHENOTYPE)
  replicated_loci_by_tract.csv       loci per tract: all, replicated in Phase 7, in ABCD, in both
  tract_enrichment_fisher.csv        Fisher test per tract for the ABCD-replicated and both-replicated sets
  trait_group_by_tract.csv           colocalized-trait group x tract among loci replicated in both cohorts
"""
import argparse
import glob
import math
import os
import re

import numpy as np
import pandas as pd

E = os.environ
ap = argparse.ArgumentParser()
ap.add_argument("--rep-dir", default=E.get("REP_DIR"))
ap.add_argument("--summary-dir", default=E.get("SUMMARY_DIR"))
ap.add_argument("--ref-dir", default=E.get("REF_DIR"))
ap.add_argument("--mod", default=E.get("MOD", "FA"))
ap.add_argument("--tag", default=E.get("TAG", "1e2"))
ap.add_argument("--prefix", default=E.get("DISC_PREFIX", "UKB_FA_phase16"))
ap.add_argument("--merge-bp", type=int, default=int(E.get("MERGE_BP", 500000)))
ap.add_argument("--cut", type=float, default=float(E.get("OVERLAP_CUT", 0.5)))
ap.add_argument("--inclusive", action="store_true", help="use >= cut instead of > cut for 'spatially replicated'")
args = ap.parse_args()
os.makedirs(args.summary_dir, exist_ok=True)
CUT = args.cut
rep_ok = (lambda x: x >= CUT) if args.inclusive else (lambda x: x > CUT)
CUT_TXT = f"{'>=' if args.inclusive else '>'}{CUT:.0%}"

# ---------------------------------------------------------------- 1-3. read the per-tract files
def tract_of(path):
    return re.search(r"Tract_([^_]+)_clump", os.path.basename(path)).group(1)

gen_files = sorted(glob.glob(os.path.join(args.rep_dir, f"{args.prefix}_Tract_*_clump_replication_by_ABCD_{args.tag}.txt")))
spat_files = sorted(glob.glob(os.path.join(args.rep_dir, f"{args.prefix}_Tract_*_clump_replication_by_ABCD_spatial_overlap_all_loci_{args.tag}.txt")))
if not gen_files:
    raise SystemExit(f"No locus-level replication files in {args.rep_dir}")
if len(gen_files) != len(spat_files):
    raise SystemExit(f"{len(gen_files)} replication files but {len(spat_files)} spatial-overlap files: steps 2-3 incomplete")

gen = []
for f in gen_files:
    d = pd.read_csv(f, sep="\t")
    d.columns = [c.upper() for c in d.columns]
    gen.append(pd.DataFrame({"TRACT": tract_of(f), "CHR": d["CHR"].astype(int), "BP": d["BP"].astype(int),
                             "SNP": d["SNP"].astype(str), "P": d["P"].astype(float),
                             "replication_gen": d["REPLICATION"].astype(str)}))
gen = pd.concat(gen, ignore_index=True)

spat = []
for f in spat_files:
    d = pd.read_csv(f, sep="\t")
    d.columns = [c.lower() for c in d.columns]
    spat.append(pd.DataFrame({"TRACT": tract_of(f), "CHR": d["chr"].astype(int), "BP": d["pos"].astype(int),
                              "disc_voxels": d["disc_voxels"], "rep_voxels": d["rep_voxels"],
                              "overlap_voxels": d["overlap_voxels"], "overlap_percent": d["overlap_percent"]}))
spat = pd.concat(spat, ignore_index=True)

# ---------------------------------------------------------------- 4. tract-level table
tl = gen.merge(spat, on=["TRACT", "CHR", "BP"], how="left", validate="one_to_one")
n_unmatched = int(tl["disc_voxels"].isna().sum())
for c in ("disc_voxels", "rep_voxels", "overlap_voxels"):
    tl[c] = tl[c].fillna(0).astype(int)
tl["overlap_percent"] = tl["overlap_percent"].fillna(0.0)

# ---------------------------------------------------------------- 5. loci across tracts (chained 500 kb merge)
tl = tl.sort_values(["CHR", "BP"], kind="mergesort").reset_index(drop=True)
ids, cur, last_chr, last_bp = [], 0, None, None
for chr_, bp in zip(tl["CHR"], tl["BP"]):
    if chr_ != last_chr or (bp - last_bp) > args.merge_bp:
        cur += 1
    ids.append(cur)
    last_chr, last_bp = chr_, bp
tl["locus_id"] = ids
g = tl.groupby("locus_id")
tl["locus_chr"] = g["CHR"].transform("first")
tl["locus_bp_min"] = g["BP"].transform("min")
tl["locus_bp_max"] = g["BP"].transform("max")
tl["locus_bp_mid"] = g["BP"].transform("mean").round().astype(int)

# ---------------------------------------------------------------- 6. locus-level summary
rows = []
for lid, d in tl.groupby("locus_id", sort=True):
    lead = d.loc[d["P"].idxmin()]
    yes = d["replication_gen"] == "Yes"
    disc, rep, ov = int(d["disc_voxels"].sum()), int(d["rep_voxels"].sum()), int(d["overlap_voxels"].sum())
    rows.append({"locus_id": lid, "snp": lead["SNP"], "chr": int(lead["CHR"]), "bp": int(lead["BP"]), "p": d["P"].min(),
                 "n_tracts": d["TRACT"].nunique(), "tracts": ",".join(sorted(d["TRACT"].unique())),
                 "gen_rep_any": bool(yes.any()), "n_tracts_gen_rep": d.loc[yes, "TRACT"].nunique(),
                 "disc_voxels_sum": disc, "rep_voxels_sum": rep, "overlap_voxels_sum": ov,
                 "overlap_percent": ov / disc if disc > 0 else np.nan})
loc = pd.DataFrame(rows).sort_values(["chr", "bp"], kind="mergesort").reset_index(drop=True)

suffix = f"{args.mod}_{args.tag}_ABCD"
out_tl = os.path.join(args.summary_dir, f"ALL_TRACTS_tract_locus_replication_merged_{suffix}.csv")
out_loc = os.path.join(args.summary_dir, f"ALL_TRACTS_locus_summary_500kb_merged_{suffix}.csv")
tl.to_csv(out_tl, index=False)
loc.to_csv(out_loc, index=False)

# ---------------------------------------------------------------- 7. counts for the manuscript
L = []
def say(s=""):
    print(s); L.append(s)
def cnt(mask, n):
    k = int(mask.sum()); return f"{k} of {n} ({100 * k / n:.1f}%)"

n = len(loc)
ovl = loc["overlap_percent"].fillna(0)
say(f"ABCD replication of {args.mod} discovery loci  (location-level P < {E.get('SIG_THRESH', '0.01')}, +-{E.get('WINDOW', '250000')} bp, merge {args.merge_bp} bp)")
say(f"results: {args.rep_dir}")
say(f"tract-level loci: {len(tl)} in {tl['TRACT'].nunique()} tracts; distinct loci: {n}")
if n_unmatched:
    say(f"WARNING: {n_unmatched} tract-level loci had no spatial-overlap row (set to 0)")
say("")
say(f"Locus-level replication (variant within the window passes in >= 1 tract): {cnt(loc['gen_rep_any'], n)}")
say(f"Spatial overlap > 0:      {cnt(ovl > 0, n)}")
say(f"Spatial overlap > 50%:    {cnt(ovl > 0.5, n)}      <- criterion stated in the Results")
say(f"Spatial overlap >= 50%:   {cnt(ovl >= 0.5, n)}")
say(f"Overlap fraction: median {ovl.median():.3f}, mean {ovl.mean():.3f}")
say(f"Criterion used for the files below: overlap {CUT_TXT}")

# ---- comparison with UKB Phase 7 (per-locus table shipped in ref/)
p7_path = os.path.join(args.ref_dir or "", "ALL_TRACTS_locus_summary_500kb_merged_FA_1e4.csv")
both = None
if args.mod == "FA" and os.path.exists(p7_path):
    p7 = pd.read_csv(p7_path)
    m = p7.merge(loc, on="locus_id", suffixes=("_p7", "_abcd"), validate="one_to_one")
    if len(m) != n or not (m["snp_p7"] == m["snp_abcd"]).all():
        say("WARNING: loci do not line up one-to-one with the UKB Phase 7 table (locus_id / lead SNP differ); "
            "the 'both cohorts' counts below are unreliable")
    o7, oa = m["overlap_percent_p7"].fillna(0), m["overlap_percent_abcd"].fillna(0)
    say("")
    say(f"UKB Phase 7 (P < 1e-4) spatially replicated, from ref/: >50%: {int((o7 > 0.5).sum())}; >=50%: {int((o7 >= 0.5).sum())}")
    say(f"Replicated in BOTH Phase 7 and ABCD, >50% in both:  {cnt((o7 > 0.5) & (oa > 0.5), n)}")
    say(f"Replicated in BOTH Phase 7 and ABCD, >=50% in both: {cnt((o7 >= 0.5) & (oa >= 0.5), n)}")
    r7, ra = rep_ok(o7), rep_ok(oa)
    lab = f"({CUT_TXT})"
    m["PHENOTYPE"] = np.select([r7 & ra, r7 & ~ra, ~r7 & ra],
                               [f"UKB + ABCD Rep {lab}", f"UKB Rep only {lab}", f"ABCD Rep only {lab}"], "No replication")
    say("")
    say(f"Supp. Fig. 24 classes (overlap {CUT_TXT}):")
    for k, v in m["PHENOTYPE"].value_counts().items():
        say(f"  {k:28s} {v}")
    m.rename(columns={"chr_p7": "CHR", "bp_p7": "POS"})[["CHR", "POS", "PHENOTYPE"]].to_csv(
        os.path.join(args.summary_dir, "loci_fa_phenogram.txt"), sep="\t", index=False)
    both = m

    # ---- Table S5 ABCD columns
    t5 = loc[["locus_id", "snp", "chr", "bp"]].copy()
    t5["ABCD_gen_rep_any"] = loc["gen_rep_any"]
    t5["ABCD_n_tracts_gen_rep"] = loc["n_tracts_gen_rep"]
    t5["ABCD_disc_voxels"] = loc["disc_voxels_sum"]
    t5["ABCD_rep_voxels"] = loc["rep_voxels_sum"]
    t5["ABCD_overlap_voxels"] = loc["overlap_voxels_sum"]
    t5["ABCD_overlap_fraction"] = loc["overlap_percent"]
    t5.to_csv(os.path.join(args.summary_dir, "table_S5_ABCD_columns.csv"), index=False)

    # ---- tract composition of the replicated sets (text: which tracts carry the reproduced loci)
    ids_all = set(m["locus_id"])
    sets = {"UKB_Phase7": set(m.loc[r7, "locus_id"]), "ABCD": set(m.loc[ra, "locus_id"]), "both": set(m.loc[r7 & ra, "locus_id"])}
    ex = p7[["locus_id", "tracts"]].assign(tract=p7["tracts"].str.split(r"\s*,\s*")).explode("tract")
    ex = ex[ex["tract"].notna() & (ex["tract"] != "")]
    by_tract = ex.groupby("tract")["locus_id"].nunique().rename("n_loci_all").to_frame()
    for k, s in sets.items():
        by_tract[f"n_loci_rep_{k}"] = ex[ex["locus_id"].isin(s)].groupby("tract")["locus_id"].nunique()
    by_tract = by_tract.fillna(0).astype(int).sort_values("n_loci_rep_ABCD", ascending=False)
    by_tract.to_csv(os.path.join(args.summary_dir, "replicated_loci_by_tract.csv"))
    say("")
    say("Loci per tract (a locus counts in every tract it spans): all / Phase 7 / ABCD / both")
    for t, r in by_tract.iterrows():
        say(f"  {t:5s} {r['n_loci_all']:4d} {r['n_loci_rep_UKB_Phase7']:4d} {r['n_loci_rep_ABCD']:4d} {r['n_loci_rep_both']:4d}")

    # ---- Fisher exact test per tract: are replicated loci over-represented among the loci of the tract?
    def fisher(a, b, c, d):
        lf = lambda x: math.lgamma(x + 1)
        lp = lambda a, b, c, d: lf(a + b) + lf(c + d) + lf(a + c) + lf(b + d) - lf(a + b + c + d) - lf(a) - lf(b) - lf(c) - lf(d)
        p0, tot = lp(a, b, c, d), 0.0
        for x in range(max(0, a + c - (c + d)), min(a + b, a + c) + 1):
            px = lp(x, a + b - x, a + c - x, d - a + x)
            if px <= p0 + 1e-9:
                tot += math.exp(px)
        return min(tot, 1.0), (a * d / (b * c) if b * c > 0 else float("inf"))
    fr = []
    for k in ("ABCD", "both"):
        s = sets[k]
        for t in sorted(ex["tract"].unique()):
            with_t = set(ex.loc[ex["tract"] == t, "locus_id"])
            a = len(s & with_t); b = len(s) - a; c = len(with_t - s); d = len(ids_all) - len(s) - c
            p, orr = fisher(a, b, c, d)
            fr.append({"replicated_set": k, "tract": t, "n_replicated_loci_in_tract": a, "n_loci_in_tract": len(with_t),
                       "odds_ratio": orr, "p_value": p})
    pd.DataFrame(fr).sort_values(["replicated_set", "p_value"]).to_csv(
        os.path.join(args.summary_dir, "tract_enrichment_fisher.csv"), index=False)

    # ---- colocalized-trait groups x tract among loci replicated in both cohorts
    coloc_path = os.path.join(args.ref_dir, "coloc_top_trait_per_locus.tsv")
    if os.path.exists(coloc_path):
        co = pd.read_csv(coloc_path, sep="\t").sort_values("PP.H4.abf", ascending=False).drop_duplicates("locus_id")
        groups = {"Metabolic / renal": ["t2d", "glomerular_filtration_rate"], "Vascular": ["stroke", "coronary_atherosclerosis"],
                  "Neurodegenerative": ["AD", "PD"], "Psychiatric": ["SCZ", "BD", "insomina"],
                  "Neurodevelopment / cognition": ["CP", "intelligence", "EA", "RT"],
                  "Behavioral": ["drinks", "number_sexual_partners"]}
        t2g = {t: gname for gname, ts in groups.items() for t in ts}
        d = m[["locus_id", "PHENOTYPE"]].merge(co[["locus_id", "trait"]], on="locus_id", how="left")
        d["trait_group"] = d["trait"].map(t2g).fillna("Other")
        tr = tl.loc[rep_ok(tl["overlap_percent"]), ["locus_id", "TRACT"]].drop_duplicates()   # tracts in which the locus is reproduced in ABCD
        d = d[d["PHENOTYPE"].str.startswith("UKB + ABCD")].merge(tr, on="locus_id", how="inner")
        tg = d.groupby(["trait_group", "TRACT"])["locus_id"].nunique().rename("n_loci").reset_index()
        tg["prop_within_group"] = tg["n_loci"] / tg.groupby("trait_group")["n_loci"].transform("sum")
        tg["n_loci_in_group"] = tg["trait_group"].map(d.groupby("trait_group")["locus_id"].nunique())
        tg = tg.sort_values(["trait_group", "prop_within_group"], ascending=[True, False])
        tg.to_csv(os.path.join(args.summary_dir, "trait_group_by_tract.csv"), index=False)
        say("")
        say("Colocalized-trait group x tract, loci replicated in both cohorts (top 3 tracts per group; share of the group's tract-level signals):")
        for gname, dd in tg.groupby("trait_group"):
            top = ", ".join(f"{r.TRACT} {r.prop_within_group:.0%}" for r in dd.head(3).itertuples())
            say(f"  {gname:30s} loci={int(dd['n_loci_in_group'].iloc[0]):3d}  {top}")
else:
    say("")
    say("No UKB Phase 7 per-locus table in ref/ for this measure: 'both cohorts' counts, the ideogram file and the tract tables were skipped.")

# ---- change relative to the previous (no-site) ABCD run, if its table is in ref/
old_path = os.path.join(args.ref_dir or "", f"ALL_TRACTS_locus_summary_500kb_merged_{args.mod}_{args.tag}_ABCD_nosite.csv")
if os.path.exists(old_path):
    old = pd.read_csv(old_path)
    oo = old["overlap_percent"].fillna(0)
    say("")
    say("Previous run (no site covariate), from ref/:")
    say(f"  locus-level {int(old['gen_rep_any'].sum())}; overlap >0: {int((oo > 0).sum())}; >50%: {int((oo > 0.5).sum())}; >=50%: {int((oo >= 0.5).sum())}")
    j = old.merge(loc, on="locus_id", suffixes=("_old", "_new"))
    a, b = j["overlap_percent_old"].fillna(0) > 0.5, j["overlap_percent_new"].fillna(0) > 0.5
    say(f"  >50% in both runs: {int((a & b).sum())}; only before: {int((a & ~b).sum())}; only now: {int((~a & b).sum())}; "
        f"correlation of overlap fractions: {j['overlap_percent_old'].fillna(0).corr(j['overlap_percent_new'].fillna(0)):.3f}")

with open(os.path.join(args.summary_dir, "abcd_replication_summary.txt"), "w") as fh:
    fh.write("\n".join(L) + "\n")
print(f"\nWrote tables to {args.summary_dir}")
