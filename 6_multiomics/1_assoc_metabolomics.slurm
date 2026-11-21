#!/bin/bash
#SBATCH --job-name=met
#SBATCH --output=/path/to/UKB_PSC/logs/multiomics/metabolomics_%a.out
#SBATCH --error=/path/to/UKB_PSC/logs/multiomics/metabolomics_%a.err
#SBATCH --mem=20G
#SBATCH --cpus-per-task=1
#SBATCH --time=00:30:00
#SBATCH --array=1-430%50

export TMPDIR=/path/to/tmp

IDX=$SLURM_ARRAY_TASK_ID
N_LDRS=$(awk -F',' -v i=${IDX} 'NR == i+1 {print $2}' /path/to/UKB_PSC_GWAS/FA/selLDRs.txt)
N_LDRS_VAL=$(awk -F',' -v i=${IDX} 'NR == i+1 {print $2}' /path/to/UKB_PSC_GWAS/FA_val/selLDRs_val.txt)

module load anaconda
cd /path/to/heig
conda activate heig

cd /path/to/tabular

# DISCOVERY
python main.py \
--ldr-assoc \
--out /path/to/UKB_PSC/multiomics/metabolomics/disc/ldr-assoc/ldr_assoc_${IDX} \
--covar /path/to/pheno_covar/2age_sex_inter_40gpc_assess_cat.txt \
--cat-covar-list sex,assess_center \
--ldrs /path/to/UKB_PSC_GWAS/FA/ldrs/UKB_FA_phase16_${IDX}_ldr_top* \
--tabular-mt /path/to/UKB_PSC/multiomics/data/metabolomics_mt_disc.mt \
--spark-conf /path/to/spark_config.json

python main.py \
--sumstats \
--out /path/to/UKB_PSC/multiomics/metabolomics/disc/ldr-sumstats/ldr_sumstats_${IDX} \
--ldr-raw-sumstats /path/to/UKB_PSC/multiomics/metabolomics/disc/ldr-assoc/ldr_assoc_${IDX}.parquet

python main.py \
--voxel-assoc \
--out /path/to/UKB_PSC/multiomics/metabolomics/disc/voxel-assoc/voxel_assoc_${IDX} \
--n-ldrs ${N_LDRS} \
--ldr-sumstats /path/to/UKB_PSC/multiomics/metabolomics/disc/ldr-sumstats/ldr_sumstats_${IDX} \
--bases /path/to/UKB_PSC_GWAS/FA/fpca/UKB_FA_phase16_${IDX}_bases_top* \
--ldr-cov /path/to/UKB_PSC_GWAS/FA/ldrs/UKB_FA_phase16_${IDX}_ldr_cov_top* \
--sig-thresh 0.05


# VALIDATION
python main.py \
--ldr-assoc \
--out /path/to/UKB_PSC/multiomics/metabolomics/val/ldr-assoc/ldr_assoc_${IDX} \
--covar /path/to/pheno_covar/2age_sex_inter_40gpc_assess_cat.txt \
--cat-covar-list sex,assess_center \
--ldrs /path/to/UKB_PSC_GWAS/FA_val/ldrs/UKB_FA_phase7_${IDX}_ldr_top* \
--tabular-mt /path/to/UKB_PSC/multiomics/data/metabolomics_mt_val.mt \
--spark-conf /path/to/spark_config.json

python main.py \
--sumstats \
--out /path/to/UKB_PSC/multiomics/metabolomics/val/ldr-sumstats/ldr_sumstats_${IDX} \
--ldr-raw-sumstats /path/to/UKB_PSC/multiomics/metabolomics/val/ldr-assoc/ldr_assoc_${IDX}.parquet

python main.py \
--voxel-assoc \
--out /path/to/UKB_PSC/multiomics/metabolomics/val/voxel-assoc/voxel_assoc_${IDX} \
--n-ldrs ${N_LDRS_VAL} \
--ldr-sumstats /path/to/UKB_PSC/multiomics/metabolomics/val/ldr-sumstats/ldr_sumstats_${IDX} \
--bases /path/to/UKB_PSC_GWAS/FA/fpca/UKB_FA_phase16_${IDX}_bases_top* \
--ldr-cov /path/to/UKB_PSC_GWAS/FA_val/ldrs/UKB_FA_phase7_${IDX}_ldr_cov_top* \
--sig-thresh 0.05