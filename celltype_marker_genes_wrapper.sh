#!/bin/bash

for CELL_TYPE in HSC EryP NeuP MasP MkP BP MoP
do
    qsub celltype_marker_genes.sh $CELL_TYPE
done
