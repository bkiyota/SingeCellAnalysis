#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N bkVIOLIN
#$ -l s_vmem=32G
#$ -l mem_req=32G
#$ -j y
#$ -o qsub_files

# qsub violin_comparisons_plots.sh 3 BP

source $HOME/setenv/miniconda.sh

DAY=$1
CELL_TYPE=$2

RNA_COUNT_PATH=/path/to/raw/test_cite_inputs.h5
METADATA_PATH=/path/to/raw/metadata.csv
ENRICHED_GENES_PATH=/path/to/processed/volcano/volcano_enriched_genes.csv.gz
OUT_DIR=/path/to/processed/volcano

python violin_comparisons_plots.py \
-r $RNA_COUNT_PATH \
-m $METADATA_PATH \
-d $DAY \
-c $CELL_TYPE \
-e $ENRICHED_GENES_PATH \
-o $OUT_DIR
