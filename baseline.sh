#!/bin/bash
../baseline.py S6short.fa 
needleall -gapopen 15 -gapextend 0.5 -asequence S6short.fa -supper1 -bsequence S6.d_baseline.fa -supper2 -aformat3 fasta -outfile S6.d_baseline.aln
../PTE_composition.py -f . -b N -s 3 -e 3 -o d_baseline