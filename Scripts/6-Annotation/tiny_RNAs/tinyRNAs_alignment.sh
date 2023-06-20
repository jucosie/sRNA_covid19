#!/bin/bash
##################################################################
# Script to align the 13-14 nt sequences to miRBase and          #
# plot the alignment to see whether they're from 5' mature miRNA #
##################################################################

results_dir="/path/to/annotation_results/tiny_RNAs"

cd $results_dir

scram_linux profile -r miRBase_mature_human_final.fa -1 sequences_tocheck_tiny.fa -t fa --minLen 13 -o tinyRNA -l 13,14
scram_plot.py profile -a tinyRNA -l 13,14 -png #If you don't need to visualize every plot comment this line, as it will run interactively a plot for every alignment

