#!/usr/bin/env python3

import sys
import pickle
import re
import os
import csv
import argparse
import pprint

import pandas as pd
from collections import defaultdict

from Bio import AlignIO, SeqIO
from Bio.Alphabet import IUPAC
from Bio.Seq import MutableSeq, translate

from ind import trim_read, findEnds, endMatch, findGap, gapAlign
from output import print_coloured_diff, printErrors

# Demand Python 3.
if sys.version_info[0] < 3:
    print("Python 3 is required, but you are using Python %i.%i.%i") % (
        sys.version_info[0], sys.version_info[1], sys.version_info[2])
    sys.exit(1)

# Identifying mutations from fasta alignment

def read_is_wt(read,ref):
    """
    Check whether read sequence is exactly equal to wt.
    :param read:
    :param ref:
    :return:
    """
    trimmed_read = re.search(r'^-+([AGCTN][ACGTN-]+[ACGTN])-+$', str(read))
    if trimmed_read is None:
        return True
    else:
        return re.search(str(trimmed_read.group(1)), str(ref)) is not None


def indel_len(sequence, start):
    l = 0
    while sequence[start + l] == '-':
        l += 1
    return l


def find_DNA_hgvs(read, ref, refname, verbose=False, start_offset=3, end_trail=3):
    """@ read, ref: MutableSeq objects
    :return errors - tuple (position, expected triplet, actual triplet, ) / none if broken read

    The assumption is that the reference includes an offset of 3 nt either side of the gene of interest. The starting triplet is
    reported as 'amino acid 0'. If the offset is less or more, it needs to be set explicitly. end_trail specifies the
    number of nt after end of gene and is ignored
    """
    if read is None:
        if verbose:
            print('no read provided')
        return

    # No gap realignment at this point
    prefix = str(refname) + ':c.'

    # quality control that there are no mutations at ends of reads
    ends = findEnds(read, ref, start_offset)
    if not endMatch(read, ref, ends):
        if verbose:
            print('ends do not match')
        return

    # scan read & reference letter by letter, counting position in reference
    # reads have been trimmed so that reference starts @ 3 (0,1,2 is the extra triplet)
    # in the general case, reference starts @ offset in 0-count
    # This is equal to the number of nt before ATG
    # ref_index denotes HGVS DNA position labeling, i is used for accessing sequence
    dna_errors = []
    ref_index = ends.get('start') - start_offset + 1  # if the read starts at 3, this becomes nt 1 (1-based as is HGVS)
    i = ends.get('start')
    max_i = len(ref) - end_trail

    while i < ends.get('end'):
        if i > max_i:  # the trailing nt are ignored when reading mutations
            break
        # check for differences
        if read[i] == ref[i]:
            ref_index += 1
            i += 1

        elif read[i] == '-':
            # start of a deletion, format depends on length
            l = indel_len(read, i)
            if ref_index > 0:
                if l == 1:  # format is POSdel
                    dna_errors.append(str(ref_index) + 'del')
                else:
                    # format is FIRST_LASTdel
                    dna_errors.append(str(ref_index) + '_' + str(ref_index + l - 1) + 'del')
            i += l
            ref_index += l

        elif ref[i] == '-':
            # start of an insertion, format is FLANK_FLANKinsSEQ
            l = indel_len(ref, i)
            if ref_index > 0:
                dna_errors.append(str(ref_index -1) + '_' + str(ref_index) + 'ins' + str(read[i:i+l]) )
            i += l

        else:
            # substitution: need to include ref. sequence in format 8A>G
            if ref_index > 0:
                dna_errors.append(str(ref_index) + str(ref[i]) + '>' + str(read[i]))
            i += 1
            ref_index += 1

    # format the result including name of sequence
    if len(dna_errors) == 1:
        dna_hgvs = prefix + dna_errors[0]
    else:
        dna_hgvs = prefix + '[' + (';').join(dna_errors) + ']'

    return dna_hgvs


def find_DNA_diff(read, ref, verbose=False, start_offset=3, end_trail=3):
    """
    @ read, ref: MutableSeq objects
    :return errors - tuple (position, expected triplet, actual triplet, ) / none if broken read

    The assumption is that the reference includes 3 nt either side of the gene of interest. The starting triplet is
    reported as 'amino acid 0'.
    As for HGVS, the starting offset and number of trailing nt are variable
    Letter by letter report mutations in NGS read, all counts 1- based in result (code in 0-count).
    - substitution: 78C = nt 78 in reference is changed to C
    - deletions: 78d6 = 6 nt deleted starting with 78: 1-77, d6, 84-end
    - insertion: 78iATC = after nt 78 inserted seq. ATC
    """

    if read is None:
        if verbose:
            print('no read provided')
        return

    # No gap realignment at this point

    # quality control that there are no mutations at ends of reads
    ends = findEnds(read, ref, start_offset)
    if not endMatch(read, ref, ends):
        if verbose:
            print('ends do not match')
        return

    # scan read & reference letter by letter, counting position in reference
    # reads have been trimmed so that reference starts @ offset=3 by default (0,1,2 is the extra triplet)
    dna_errors = []
    ref_index = ends.get('start') - start_offset + 1
    i = ends.get('start')
    max_i = len(ref) - end_trail

    while i < ends.get('end'):
        if i > max_i:
            break
        # check for differences
        if read[i] == ref[i]:
            ref_index += 1
            i += 1

        elif read[i] == '-':
            # start of a deletion
            l = indel_len(read, i)
            # now we know the length of a deletion, check for frameshifts
            if l % 3 == 0:
                if ref_index > 0:
                    dna_errors += [(str(ref_index), 'd', str(l))]  # deletion length l starting at ref_index in 0-count
                i += l
                ref_index += l
            else:
                dna_errors += [(str(ref_index), 'f')]
                break

        elif ref[i] == '-':
            # start of an insertion
            l = indel_len(ref, i)
            # check for frameshifts
            if l % 3 == 0:
                if ref_index > 0:
                    dna_errors += [(str(ref_index), 'i', str(read[i:i+l]) )]
                i += l
            else:
                dna_errors += [(str(ref_index), 'f')]
                break

        else:
            # substitution
            if ref_index > 0:
                dna_errors += [(str(ref_index + 1), 's', str(read[i]) )]
            i += 1
            ref_index += 1

    return tuple(dna_errors)


def find_protein_diff(read, ref, verbose=False, start_offset=3, end_trail=3):

    # quality control
    if read is None:
        return
    ends = findEnds(read, ref, start_offset)
    if not endMatch(read, ref, ends):
        return

    newread = read
    newref = ref

    # scan reference triplet by triplet
    # move letters when encountering an indel
    prot_errors = []
    i = ends.get('aligned')
    ref_index = int((ends.get('aligned') - start_offset)/3) + 1  # reference amino acid index
    max_i = len(ref) - end_trail

    while i <= ends.get('end'):
        if i > max_i:
            break

        if newread is None:
            break
        ref_codon = newref[i:i+3]
        read_codon = newread[i:i+3]

        if '-' in read_codon:  # found a deletion
            # Check if this is the last acid, and it's incomplete, ignore it.
            if re.search('[ATGC]', str(newread[i + 3:])) is None:
                break

            if '-' in ref_codon:  # something very broken
                prot_errors.append((ref_index,'f'))
                return tuple(prot_errors)
            elif read_codon == '---':  # single codon deletion
                if ref_index > 0:
                    prot_errors += [(ref_index, 'd')]
                i += 3
                ref_index += 1

            else:  # check it's not a frame shift
                l = indel_len(newread, i)
                if l % 3 != 0:
                    prot_errors.append((ref_index, 'f'))
                    return tuple(prot_errors)
                # realign gap and repeat loop at same position to compare the codons
                gap = findGap(newread[i - 1:])
                gap = (gap[0] + i - 1, gap[1] + i - 1)
                newread = gapAlign(newread, gap, start_offset)
                continue

        elif '-' in ref_codon:  # found an insertion
            l = indel_len(newref, i)
            if l % 3 != 0:
                prot_errors.append((ref_index, 'f'))
                return tuple(prot_errors)
            gap = findGap(newref[i-1:])
            if gap[0] == 1:  # insertion after codon
                insertion = newread[gap[0] + i - 1:gap[1] + i - 1]
                if '-' in insertion:
                    prot_errors.append((ref_index, 'f'))
                    return tuple(prot_errors)
                if ref_index > 0:
                    prot_errors.append((ref_index, 'i', str(translate(insertion)) ))
                i += l
                ref_index += 1
            else:  # realign gap and repeat loop at same position to compare the codons
                gap = (gap[0] + i - 1, gap[1] + i - 1)
                newref = gapAlign(newref, gap, start_offset)
                continue

        elif translate(read_codon) != translate(ref_codon):  # must be a substitution
            if ref_index > 0:
                prot_errors.append((ref_index, 's', str(translate(read_codon))))
            if str(translate(read_codon)) == '*':
                return tuple(prot_errors)
            i += 3
            ref_index += 1

        else:
            i += 3
            ref_index += 1

    if verbose:
        print(prot_errors)

    return tuple(prot_errors)


# Raw processing of all alignments, get composition

def count_one_fraction(alignment, refname, debug, start_offset, end_trail):
    """
    Don't bother with expected/allowed mutations, just find everything and filter later
    Final format: {DNA error: [(protein error), fraction,
    1. Read reference file
    2. Scan over reference sequence to generate all possible mutations
    3. For each ref & read in multiple alignment:
        - verify the read is good quality
        - call the mutation
        - add to count table
    4. Print counts
    """
    # use a regular dictionary
    # when a protein mutation is first encountered, create an entry
    one_lane_counts = {}

    # reading & looping over read/reference sequence in multiple sequence alignment
    # use AlignIO parser and keep sequence only, allowing it to change (important for gap shifts) 
    for pair in AlignIO.parse(alignment, "fasta", alphabet=IUPAC.ambiguous_dna, seq_count=2):
        # both read and ref are MutableSeq
        ref = pair[0].seq.tomutable()
        read = pair[1].seq.tomutable()
        read = MutableSeq(str(read).replace('N', '.'), read.alphabet)
        readname = pair[1].id

        # trim sequencing read to reference
        ref, read = trim_read(ref, read)

        if read_is_wt(read,ref):
            if debug:
                trimmed_read = re.search(r'^-+([AGCTN][ACGTN-]+[ACGTN])-+$', str(read))
                print()
                print(trimmed_read.group(1))
                printErrors("WT", read, ref, True)
            continue

        dna_errors, dna_hgvs, prot_erros = None, None, None

        try:
            dna_errors = find_DNA_diff(read, ref, debug, start_offset, end_trail)  # errors = a tuple
            dna_hgvs = find_DNA_hgvs(read, ref, refname, debug, start_offset, end_trail)  # string according to HGVS format (ish)
            prot_errors = find_protein_diff(read, ref, debug, start_offset, end_trail)
            # print()
            # print(readname)
            # print(dna_hgvs, prot_errors)
            # printErrors(dna_errors, read, ref, True)

        except:
            if not dna_errors:
                print(dna_errors)
            print_coloured_diff(readname, read, ref, debug)
            raise

        try:
            one_lane_counts[prot_errors]['total'] += 1
            one_lane_counts[prot_errors]['dna'][dna_errors] += 1
            one_lane_counts[prot_errors]['dna_hgvs'][dna_hgvs] += 1
        except KeyError:
            one_lane_counts[prot_errors] = {'dna': defaultdict(int), 'dna_hgvs': defaultdict(int), 'total': 1}
            one_lane_counts[prot_errors]['dna'][dna_errors] += 1
            one_lane_counts[prot_errors]['dna_hgvs'][dna_hgvs] += 1


    # count the mutations
    n = 0
    threshold = 10
    for error in one_lane_counts.keys():
        if one_lane_counts[error]['total'] > threshold:
            n += 1

    print('Fount {0} total protein mutations, of which {1} have more than {2} counts'
          .format(len(one_lane_counts), n, threshold))

    return one_lane_counts


def count_multiple_fractions(folder, baseline, debug, start_offset, end_trail):
    """
    Process all reference.fraction.aln files in given  folder
    :param folder: contains all *.aln
    :param baseline: string containing name of baseline
    :return:
    """
    all_references = {}

    print(os.listdir(folder))
    for f in os.listdir(folder):
        if f.endswith('.aln'):
            aln_path = os.path.join(folder, f)
            refname, fraction, suffix = f.rsplit(".", 2)
            print('Counting alignment {0} in background {1} and activity fraction {2}'
                  .format(f, refname, fraction))

            if refname not in all_references.keys():
                all_references[refname] = {}

            if fraction == baseline:
                fraction = 'baseline'
            all_references[refname][fraction] = count_one_fraction(aln_path, refname, debug, start_offset, end_trail)

    return all_references





# Start composition statistics


if __name__ == "__main__":
    """
    Arguments:
    1.   Alignment of sequencing reads to a reference sequence in FASTA format. The
         reference comes before read sequence. (>Ref, seq, >Read, seq).
    (2.) Optional: DEBUG print a coloured representation of mismatches
    """
    parser = argparse.ArgumentParser(description='Finds all in-frame mutations in a gene')
    parser.add_argument('-f', '--folder', help='Folder containing multiple sequence alignments', required=True)
    parser.add_argument('-r', '--reference', required=False)
    parser.add_argument('-d', '--debug', help='Turn on debugging', required=False, action="store_true")  # Visual
    parser.add_argument('-b', '--baseline', help='Name of baseline fraction', required=False)
    parser.add_argument('-s', '--start_offset', help='Number of nt before starting ATG (integer)', required=False,
                        default=3, type=int)
    parser.add_argument('-e', '--end_trail', help='Number of nt after end of gene (integer)', required=False, default=0,
                        type=int)
    parser.add_argument('-o', '--output', help='Filename')
    args = parser.parse_args()


    # On the first run, analyse all *.aln files in the target folder and create a dictionary of errors
    # Structure: all_ref[background][fraction][protein mutation] = all data about the mutations
    all_ref = count_multiple_fractions(args.folder, args.baseline, args.debug, args.start_offset, args.end_trail)
    export_hgvs(all_ref, 'new')
    with open(args.output +'.p', 'wb') as f:
        pickle.dump(all_ref, f)