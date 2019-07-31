#!/usr/bin/python3

"""
1. Create a reference sequence
2. Generate all possible 3/6/9 bp deletions on them as reads
3. Print in format compatible with needle
(4. Process with needle in case alignment reshuffles bases)
(5. count_one_lane.py)
"""

import sys

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
    reference = SeqIO.read(sys.argv[1],'fasta', alphabet = IUPAC.ambiguous_dna)
    sequence = str(reference.seq.upper())
    deletion = [9]
    insertion = [6]
    codons = list(CodonTable.unambiguous_dna_by_name["Standard"].forward_table.keys())
    codons.sort()

    start_offset = 6
    end_offset = 3

    # for length in deletion:
    #     with open(reference.name + ".base_d" + str(length) + ".fa", "w") as output:
    #         for i in range(start_offset - 1, len(sequence) - length - end_offset):
    #             newseq = sequence[:i] + sequence[i + length:]
    #             record = SeqRecord(Seq(newseq, reference.seq.alphabet),
    #                                id="d" + str(length) + '_' + str(i), description="")
    #             SeqIO.write(record, output, "fasta")

    with open(reference.name + ".base_i6.fa", "w") as output:
        for i in range(start_offset - 1, len(sequence) - 6 - end_offset):
            for triplet_1 in codons:
                for triplet_2 in codons:
                    newseq = sequence[:i] + triplet_1 + triplet_2 + sequence[i:]
                    record = SeqRecord(Seq(newseq, reference.seq.alphabet),
                                       id="insertion6", description="")
                    SeqIO.write(record, output, "fasta")


main()
