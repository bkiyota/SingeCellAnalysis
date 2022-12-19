#!/bin/bash
#$ -S /bin/sh
#$ -N bkMARKERS
#$ -cwd
#$ -l s_vmem=64G
#$ -l mem_req=64G
#$ -j y
#$ -o qsub_files

# qsub celltype_marker_genes.sh

CELL_TYPE=$1

source $HOME/setenv/miniconda.sh

cmd='which Rscript'

Rscript celltype_marker_genes.R $CELL_TYPE

