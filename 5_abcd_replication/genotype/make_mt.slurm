#!/bin/bash
#SBATCH --job-name=make_mt
#SBATCH --output=/path/to/UKB_PSC/logs/make_mt_abcd.out
#SBATCH --error=/path/to/UKB_PSC/logs/make_mt_abcd.err
#SBATCH --mem=30G
#SBATCH --cpus-per-task=4
#SBATCH --partition=general
#SBATCH --time=10:00:00

module purge
module load anaconda
conda activate heig
cd /path/to/heig

python heig.py \
--make-mt \
--out /path/to/UKB_PSC/output/genotype/genotype_abcd_qc \
--bfile /path/to/ABCD/genotype/ABCD_release_2.0.1_r1 \
--grch37 \
--keep /path/to/UKB_PSC/Modalities/ABCD/abcd_ids.txt \
--spark-conf /path/to/spark_config.json \
--threads 4