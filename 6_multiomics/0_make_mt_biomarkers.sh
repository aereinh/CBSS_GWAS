#!/bin/bash

module load anaconda

cd /path/to/heig

conda activate heig

cd /path/to/tabular

python main.py \
--make-mt \
--tabular-txt /path/to/UKB_PSC/multiomics/data/biomarkers_dat_disc.txt \
--out /path/to/UKB_PSC/multiomics/data/biomarkers_mt_disc \
--spark-conf /path/to/spark_config.json

python main.py \
--make-mt \
--tabular-txt /path/to/UKB_PSC/multiomics/data/biomarkers_dat_val.txt \
--out /path/to/UKB_PSC/multiomics/data/biomarkers_mt_val \
--spark-conf /path/to/spark_config.json