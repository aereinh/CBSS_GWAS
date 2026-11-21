#!/bin/bash
# What is in place for the ABCD replication, and how far each step got.
REP_SCRIPTS="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"; source "${REP_SCRIPTS}/rep_config.sh"
shopt -s nullglob
# count files matching a glob (an unmatched glob expands to nothing, so ls would list the current folder)
nf() { local f=( $1 ); echo ${#f[@]}; }
ok() { [[ -e "$1" ]] && echo "  [ OK ] $2" || echo "  [MISS] $2: $1"; }
echo "-- inputs"
ok "${VGWAS_HEIG}" "heig.py for voxel GWAS"; ok "${CLUSTER_TRACT_CSV}" "cluster -> tract table"; ok "${SEL_LDRS}" "selLDRs_val.txt"
n=$(nf "${SUMSTATS_DIR}/ABCD_${MOD}_*.snpinfo");            echo "  [$n/430] ABCD LDR sumstats (${SUMSTATS_DIR})"
n=$(nf "${LDR_DIR}/ABCD_${MOD}_*_ldr_cov_top*");            echo "  [$n/430] ABCD LDR covariance files (${LDR_DIR})"
n=$(nf "${DISC_TRACT_DIR}/${DISC_PREFIX}_*_candidate_sig_snps.txt"); echo "  [$n/17] discovery candidate SNP lists"
n=$(nf "${DISC_TRACT_DIR}/${DISC_PREFIX}_Tract_*_clump_r2_0.1_step3.clumped"); echo "  [$n/17] discovery clump files"
n=$(nf "${DISC_VOX_DIR}/${DISC_PREFIX}_*_sig_vgwas_screen.txt"); echo "  [$n] discovery footprint files (${DISC_VOX_DIR})"
ok "${REF_DIR}/ALL_TRACTS_locus_summary_500kb_merged_FA_1e4.csv" "UKB Phase 7 per-locus table (ref/)"
echo "-- step 1 (voxel GWAS)"
n=$(nf "${REP_DIR}/ABCD_${MOD}_*_vgwas_${TAG}.done"); m=$(nf "${REP_DIR}/ABCD_${MOD}_*_vgwas_${TAG}.txt")
echo "  [$n/430] clusters finished; $m with associations below ${SIG_THRESH}; $(grep -l 'no tract' "${REP_DIR}"/*.done /dev/null 2>/dev/null | wc -l) without a tract"
echo "-- steps 2-3 (per tract)"
n=$(nf "${REP_DIR}/${DISC_PREFIX}_Tract_*_clump_replication_by_ABCD_${TAG}.txt"); echo "  [$n/17] locus-level replication files"
n=$(nf "${REP_DIR}/${DISC_PREFIX}_Tract_*_clump_replication_by_ABCD_spatial_overlap_all_loci_${TAG}.txt"); echo "  [$n/17] spatial overlap files"
echo "-- step 4 (aggregate)"
ok "${SUMMARY_DIR}/abcd_replication_summary.txt" "summary"; [[ -e "${SUMMARY_DIR}/abcd_replication_summary.txt" ]] && sed -n '1,40p' "${SUMMARY_DIR}/abcd_replication_summary.txt"
