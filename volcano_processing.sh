#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N bkVOLCANO
#$ -l s_vmem=64G
#$ -l mem_req=64G
#$ -j y
#$ -o qsub_files

# qsub volcano_processing.sh 2 BP

source $HOME/setenv/miniconda.sh

DAY=$1
CELL_TYPE=$2

RNA_COUNT_PATH=/path/to/raw/test_cite_inputs.h5
METADATA_PATH=/path/to/raw/metadata.csv
OUT_DIR=/path/to/processed/volcano

python volcano_processing.py \
-r $RNA_COUNT_PATH \
-m $METADATA_PATH \
-d $DAY \
-c $CELL_TYPE \
-o $OUT_DIR
