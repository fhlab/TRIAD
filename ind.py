"""
All functions for finding substitutions, insertions & deletions
"""

import re
import glob
import os
import csv
import pickle
import numpy as np

from Bio.Emboss.Applications import NeedleallCommandline

from collections import defaultdict

import Bio.AlignIO
import Bio.SeqIO
import Bio.Data.CodonTable
import Bio.Alphabet.IUPAC
from Bio.SeqRecord import SeqRecord
from Bio.Seq import MutableSeq, Seq


"""
STRING UTILS
"""


def get_char(target, index):
    """
    Get a char from a string, or "-" if `index` is out of range.
    """
    if index < 0 or index >= len(target):
        return "-"

    return target[index]


def get_chars(target, index, n):
    """
    Get a sequence of `n` chars from `target` at `index`, safely, via `getChar`
    """

    out = ""
    for i in range(0, n):
        out += get_char(target, index + i)

    return out


def set_char(target_str, char, i):
    """
    Set the character at index `i` in `targetStr` to `char`, returning the new string.
    """

    return target_str[:i] + char + target_str[i + 1:]


"""
INDIVIDUAL READS
Check quality, match ends, find positions of indels
"""


def findEnds(read, ref, start_offset):
    """
    Figure out where the first symbol of interest in the read is. This start point should be:
    - Not a -.
    - Not (yet) triplet aligned
    Similarly find the last non-dash symbol in the read.
    """
    ends = {"start": 0, "end": len(ref)}
    for i in range(0, len(ref)):
        if read[i] != "-":
            ends["start"] = i
            ends["aligned"] = i
            break

    for i in range(0, len(ref)):
        if read[-1 - i] != "-":
            ends["end"] = len(ref) - i
            break

    while ends.get("aligned") % 3 != start_offset % 3:
        ends["aligned"] += 1

    return ends


def endMatch(read, ref, ends, MATCH_N_END=3):
    """
    Aligner errors arise when mutations are at the ends of the read rather than in the middle.
    Trimming ends only shifts the problem. Instead require that read ends match the reference.
    Return False if either end doesn't match
    How many nucleotides need to match at the end of a read for a valid alignment:
    - matching 2 should correct alignment errors, 3 avoids problems with InDel repositioning
    - 3 simplifies the first codon: it's either complete or it's OK to move 1 or 2 bases over to the next triplet
    """
    startRead = read[ends.get("start"):ends.get("start") + MATCH_N_END]
    startRef = ref[ends.get("start"):ends.get("start") + MATCH_N_END]
    endRead = read[ends.get("end") - MATCH_N_END:ends.get("end")]
    endRef = ref[ends.get("end") - MATCH_N_END:ends.get("end")]
    start = startRead == startRef
    end = endRead == endRef
    return start and end


def hasMultipleGaps(read):
    """
    Return true iff read contains multiple gaps.
    """
    # This regex looks for 3 distinct islands of letters. That implies at least two gaps.
    return re.search('[ATGC]+-+[ATGC]+-+[ATGC]+', str(read)) != None


def findGap(read):
    """
    Get the start and end index (as a 2-tuple) of the gap in readMatch, or None if there isn't one.
    """
    # Find the gap, captured in group 1. This gives us the indexes in the string where the gap starts and
    # ends. (we define "gap" as "dashes in between letters". We previously used hasMultipleGaps to reject
    # any strings which have multiple gaps (which would confuse this regex)).
    match = re.search('[ATGC]+(-+)[ATGC]+', str(read))
    if match == None:
        # No gap exists.
        return None

    return (match.start(1), match.end(1))


def gapAlign(read, gap, start_offset):
    """
    Perform the "letter-stealing" operation:
    - Find the gap (if there is one).
    - While the gap startpoint is not aligned to a triplet boundary, move letters from the end of the gap
      to the start of the gap.

    @param read Aligned read to process.
    @param gap The gap, as returned by findGap.
    """

    if gap is None:
        return read

    movingGap = (gap[0], gap[1])
    newread = read

    # Shift letters from the end to the start...
    while movingGap[0] % 3 != start_offset % 3:
        assert (newread[movingGap[0]] == "-")

        # This means we ran out of symbols to steal before we managed to align.
        # This indicates a gap at the end, and one we can't fix. Reject!
        if read[movingGap[1]] == "-":
            return None

        newread[movingGap[0]] = newread[movingGap[1]]
        newread[movingGap[1]] = "-"

        # Shift the gap to the right...
        movingGap = (movingGap[0] + 1, movingGap[1] + 1)

    return newread


def trim_read(ref, read):
    """
    In Sanger sequencing the read will likely extend beyond the gene reference. Trim both to reference length.
    :param ref: MutableSeq for the Sanger read
    :param read: MutableSeq with the reference sequence
    :return: trimmed ref & read
    """
    start, end = 0, len(ref)

    # trim start
    for i in range(len(ref)):
        if ref[i] != '-':
            start = i
            break
    for i in range(len(ref), 1, -1):
        if ref[i-1] != '-':
            end = i
            break

    return ref[start:end], read[start:end]
