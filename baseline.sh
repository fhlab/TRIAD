#!/bin/bash
# ../baseline.py S6-6.fa

FILES=./baseline/S6.base_i6.fa

for fa_file in $FILES
do
    echo $"Processing $fa_file"
    basename=${fa_file%.fa}

    needleall -gapopen 15 -gapextend 0.5 -asequence ./baseline/S6-6.fa -supper1 \
    -bsequence $fa_file -supper2 -aformat3 fasta -outfile ${basename}.aln

#../PTE_composition.py -f . -b N -s 6 -e 6 -o d_baseline
done