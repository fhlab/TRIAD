#!/bin/bash

# wrapper script to generate +9 bp baseline library
# split the gene (1002 bp long) into 20 bp long sections
# 2 bp section needs 52 minutes
# expect one 20 bp section to need ~5 hours to align with needleall
# 20 bp x 10 at oncex 10 hours per batch x 10 batches = 50 hours just go get alignments

# 1. generate short fa files for all sections and align the results
n_pos=10

# one batch
for i in {0..90..$n_pos}: # last step 
    generate_baseline_i9.py -r full_fragement.fa -s 200 -e 97 -o S6.base_i9_step${n_pos}_${i}.fa --position_start $i --position_step $n_pos
    needleall -gapopen 15 -gapextend 0.5 -asequence ~/TRIAD/baseline/full_fragment.fa  -supper1 \
    -bsequence S6.base_i9_step${n_pos}_${i}.fa -supper2 -aformat3 fasta -outfile S6.base_i9_step${n_pos}_${i}.aln
