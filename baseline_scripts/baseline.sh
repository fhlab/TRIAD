#!/bin/bash
# ../baseline.py S6-6.fa

FILES=./baseline/S6.base_i6.fa

for fa_file in $FILES
do
    echo $"Processing $fa_file"
    basename=${fa_file%.fa}

    needleall -gapopen 15 -gapextend 0.5 -asequence ./baseline/S6-6.fa -supper1 -bsequence $fa_file -supper2 -aformat3 fasta -outfile ${basename}.aln

#../PTE_composition.py -f . -b N -s 6 -e 6 -o d_baseline
done

# TIMING: +9 bp, for 50 positions (×20 for full gene)
# ~/TRIAD/baseline/generate_baseline_i9.py ~/TRIAD/baseline/full_fragment.fa 0  357.80s user 9.00s system 95% cpu 6:23.08 total

# TIMING: +9 bp, for 2 positions (×50 for full gene)
# variant generation:
# ~/TRIAD/baseline/generate_baseline_i9.py ~/TRIAD/baseline/full_fragment.fa 0  19.92s user 0.50s system 94% cpu 21.515 total 
# needle alignment:
# needleall -gapopen 15 -gapextend 0.5 -asequence ~/TRIAD/baseline/full_fragment.fa  -supper1 -bsequence $fa_file -supper2 -aformat3 fasta -outfile ${basename}.aln
# needleall -gapopen 15 -gapextend 0.5 -asequence  -supper1 -bsequence $fa_file  2786.62s user 69.09s system 91% cpu 52:08.29 total

# est. 50 hours for full gene
