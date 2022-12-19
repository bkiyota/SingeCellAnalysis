#!/bin/bash

for DAY in 2 3 4
do
    for CELL_TYPE in BP EryP HSC MasP MkP MoP NeuP 
    do
        qsub volcano_processing.sh $DAY $CELL_TYPE
    done
done
