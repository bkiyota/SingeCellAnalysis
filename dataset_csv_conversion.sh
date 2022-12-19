#!/bin/bash
#$ -S /bin/bash
#$ -cwd
#$ -N bkCONVERSION
#$ -l s_vmem=64G
#$ -l mem_req=64G
#$ -j y
#$ -o qsub_files

# qsub dataset_csv_conversion.sh

source $HOME/setenv/miniconda.sh

python dataset_csv_conversion.py
