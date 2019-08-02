#!/usr/bin/python3

"""
Because these insertions are longer, do not embed them within the full-length gene
sequence. Instead, leave 40 bp flanking sequence on each side.
Hence, use reference file full_fragment.fa (located in TRIAD/baseline).


Correspondence between two references:
position of A in first codon ATG:
S6-6.fa : ref[6] - there are 6 additional bases before
full_fragment.fa: ref[200] - 200 nt come before
stop codon: 97 bases after

Parallelization:
split the gene into 10 100-bp regions and generate insertions within those
Run as:

n_pos=20
for i in range(0, 1005, n_pos): # last step 
    generate_baseline_i9.py -r full_fragement.fa -s 200 -e 97 --position_step $i --position_step $n_pos


1. Generate fasta reads for the fragement
2. Align with needleall
5. Run PTE_composition.py to extract mutations
"""

import sys
import argparse

from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.Data import CodonTable


# Demand Python 3.
if sys.version_info[0] < 3:
    print("Python 3 is required, but you are using Python %i.%i.%i") % (
           sys.version_info[0], sys.version_info[1], sys.version_info[2])
    sys.exit(1)

def main():
    """
    Read reference, generate all possible sub/del reads and write them to file
    """
    parser = argparse.ArgumentParser(description='')
    parser.add_argument('-r', '--reference', required=False)
    parser.add_argument('-s', '--start_offset', help='Number of nt before starting ATG (integer)', required=False, default=3, type=int)
    parser.add_argument('-e', '--end_trail', help='Number of nt after end of gene (integer)', required=False, default=0, type=int)
    parser.add_argument('--position_start', required=True, type=int)
    parser.add_argument('--position_step', required=True, type=int)
    parser.add_argument('-o', '--output', help='Filename')
    args = parser.parse_args()
    
    
    reference = SeqIO.read(args.reference,'fasta', alphabet = IUPAC.ambiguous_dna)
    sequence = str(reference.seq.upper())
    codons = list(CodonTable.unambiguous_dna_by_name["Standard"].forward_table.keys())
    codons.sort()
  
    generation_start =  args.start_offset - 1 + args.position_start
    generation_end = min(len(sequence) - 9 - args.end_trail, generation_start + args.position_step)
    
    with open(args.output, "w") as output:
        for i in range(generation_start, generation_end):
            for triplet_1 in codons:
                for triplet_2 in codons:
                    for triplet_3 in codons:
                        newseq = sequence[i-50:i] + triplet_1 + triplet_2 + triplet_3 + sequence[i:i+50]
                        record = SeqRecord(Seq(newseq, reference.seq.alphabet),
                                           id="insertion9", description="")
                        SeqIO.write(record, output, "fasta")


main()
