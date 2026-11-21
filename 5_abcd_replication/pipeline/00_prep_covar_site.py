#!/usr/bin/env python3
"""Add ABCD site (integer 1-22, from the earliest visit) to the existing ABCD covariate file.

Input : AgeGender_ABCD_withPCs.txt  (FID IID sex age age2 agesex age2sex PC1..PC40)
Output: AgeGender_ABCD_withPCs_site.txt  (same + site)
10,481 of 10,496 subjects were matched.
Run: module load anaconda; python 00_prep_covar_site.py
"""
import pandas as pd

cov_fp  = "/path/to/UKB_PSC/covar/AgeGender_ABCD_withPCs.txt"
site_fp = "/path/to/ABCD/phenotype/cognition_confounders_substance_v013.csv"
out_fp  = "/path/to/UKB_PSC/covar/AgeGender_ABCD_withPCs_site.txt"

SITE_ID_COL, SITE_COL, SESSION_COL = "participant_id", "site", "session_id"

# How AgeGender_ABCD_withPCs.txt was made:
#   age/sex from /path/to/UKB_PSC/covar/AgeGender_ABCD.csv
#     sex = 1 (M) / 0 (F); age z-standardized; age2 = age^2; agesex = age*sex; age2sex = age2*sex
#   PC1-PC40 from /path/to/ABCD/pcs/ABCD_allchr_unimputed_sep10_2019.eigenvec
#   inner join on IID.


def id_key(s):
    """'sub-003RTV85' (site file) and 'NDAR_INV003RTV85' (covariates) -> '003RTV85'."""
    return (s.astype(str).str.strip().str.upper()
             .str.replace(r"^SUB-", "", regex=True)
             .str.replace(r"^NDAR_?INV", "", regex=True))


cov = pd.read_csv(cov_fp, sep=r"\s+")
cov["IID"] = cov["IID"].astype(str).str.strip()
print(f"Existing covariate subjects: {cov.shape[0]} ({cov.shape[1]} columns)")

site = pd.read_csv(site_fp, low_memory=False)[[SITE_ID_COL, SITE_COL, SESSION_COL]].copy()
site[SITE_COL] = site[SITE_COL].astype(str).str.strip()
site = site[~site[SITE_COL].isin(["", "nan", "NaN", "None"])]
site["key"] = id_key(site[SITE_ID_COL])

# earliest visit per subject: ses-00A < ses-02A < ses-04A ...
ses = site[SESSION_COL].astype(str).str.strip().str.upper().str.extract(r"^SES-(\d+)([A-Z]*)$")
site["_ses_num"] = pd.to_numeric(ses[0], errors="coerce")
site["_ses_letter"] = ses[1].fillna("")
site = (site.sort_values(["key", "_ses_num", "_ses_letter"], na_position="last")
            .drop_duplicates("key", keep="first"))
print("Sessions used:", site[SESSION_COL].value_counts().to_dict())

# site01 -> 1 ... site22 -> 22
site["site"] = pd.to_numeric(site[SITE_COL].str.extract(r"(\d+)")[0], errors="coerce")
site = site.dropna(subset=["site"])
site["site"] = site["site"].astype(int)

cov["key"] = id_key(cov["IID"])
merged = cov.merge(site[["key", "site"]], on="key", how="inner").drop(columns="key")
print(f"Merged subjects: {merged.shape[0]} of {cov.shape[0]}")
print(merged["site"].value_counts().sort_index().to_string())

merged.to_csv(out_fp, sep=" ", index=False)
print(f"Wrote {out_fp}")
