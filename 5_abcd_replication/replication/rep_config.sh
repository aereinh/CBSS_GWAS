#!/bin/bash
# =====================================================================
# ABCD replication of the UKB discovery loci (FA, site-adjusted).
# Voxel-level reconstruction, locus-level replication, spatial overlap and
# locus-level aggregation.  Paths of the ABCD GWAS come from ../pipeline/config.sh.
# =====================================================================
REP_SCRIPTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
export REP_SCRIPTS
export PIPE_DIR="${PIPE_DIR:-$(cd "${REP_SCRIPTS}/../pipeline" && pwd)}"
source "${PIPE_DIR}/config.sh"

export MOD="${MOD:-FA}"
set_modality "${MOD}"          # OUT_DIR, LDR_DIR, SEL_LDRS, SUMSTATS_DIR, UKB_BASES_PREFIX, RUN_NAME ...

# ---- replication settings (as in the manuscript)
export SIG_THRESH=0.01         # nominal location-level threshold in ABCD
export TAG=1e2                 # file tag for that threshold
export WINDOW=250000           # +-250 kb around each discovery lead variant
export MERGE_BP=500000         # loci of different tracts within 500 kb = one locus
export OVERLAP_CUT=0.5         # spatial replication: overlap fraction > 0.5 (strict)

# ---- discovery inputs (UKB Phases 1-6 discovery)
export DISC_DIR=/path/to/UKB_PSC_GWAS/FA
export CLUSTER_TRACT_CSV=${DISC_DIR}/cluster_tract_summary.csv
export DISC_TRACT_DIR=${DISC_DIR}/post_screen/disc/tract        # *_candidate_sig_snps.txt, *_clump_r2_0.1_step3.clumped
export DISC_VOX_DIR=/path/to/UKB_PSC_GWAS/shared/Voxel/FA/phase1to6   # *_sig_vgwas_screen.txt
export DISC_PREFIX=UKB_FA_phase16

# ---- HEIG used for the voxel-level reconstruction: the same copy as the rest of the ABCD pipeline
# (HEIG_DIR in ../pipeline/config.sh).
export VGWAS_HEIG=${HEIG_DIR}/heig.py

# ---- outputs
export REP_DIR=${OUT_DIR}/replication/tract          # per-cluster voxel GWAS + per-tract replication files
export SUMMARY_DIR=${OUT_DIR}/replication/summary    # locus-level tables and counts
export REP_LOG_DIR=${LOG_ROOT}/${RUN_NAME}/replication
export REF_DIR=${REP_SCRIPTS}/ref                    # reference tables shipped with these scripts

TRACTS=(ACR ALIC BCC CGH CP EC FX GCC IFO PCR PLIC PTR RLIC SCC SCR SLF SS)
