#!/bin/bash

for DAY in 2 3 4
do
    for CELL_TYPE in BP EryP HSC MasP MkP MoP NeuP 
    do
        qsub violin_comparisons_plots.sh $DAY $CELL_TYPE
    done
done
