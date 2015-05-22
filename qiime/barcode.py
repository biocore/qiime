#!/usr/bin/env python
import numpy

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]  # remember to add yourself
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

"""
provides functions, used to assign a possibly error-fraught DNA
barcode to a list of original barcodes.  correct_barcode uses edit distance
of the DNA, not the bit encoding of the DNA.  correct_barcode_bitwise
uses bit encoding to determine closest match
"""
DEFAULT_GOLAY_NT_TO_BITS = {"A": "11", "C": "00", "T": "10", "G": "01"}
DEFAULT_HAMMING_NT_TO_BITS = {"A": "11", "C": "10", "T": "00", "G": "01"}


def correct_barcode(query_seq, seq_possibilities):
    """ finds closest (by nt seq edit distance) match to query_seq

    assumes:
    all sequences are same length
    no sequence appears twice in seq_possibilities

    returns (best_hit, min_dist)
    * best_hit is closest sequence from seq_possibilities, or None if a tie
    * min_dist is the edit distance between the query_seq and the best hit
    cw1 = AAACCCGGGTTT (12 nucleotides)
    cw2 = AAACCCGGGTTA (edit distance of 1 from cw1)
    cw3 = AAACCCGGGTTG (edit distance of 1 from cw1)
    cw4 = AAACCCGGGTTC (edit distance of 1 from cw1)
    cw5 = AAACCCGGGTGG (edit distance of 2 from cw1)
    """
    dists = [_edit_dist(query_seq, seq) for seq in seq_possibilities]
    min_dist = min(dists)
    number_mins = dists.count(min_dist)
    if number_mins > 1:
        return None, min_dist
    else:
        best_hit = seq_possibilities[dists.index(min_dist)]
        return best_hit, min_dist


def _edit_dist(s1, s2):
    """ computes edit (hamming) between to strings of equal len

    designed for strings of nucleotides, not bits"""
    dist = 0
    for i in range(len(s1)):
        if s1[i] != s2[i]:
            dist += 1
    return dist


def correct_barcode_bitwise(query_seq, seq_possibilities,
                            nt_to_bits=DEFAULT_GOLAY_NT_TO_BITS):
    """ finds closest (by bit distance) match to query_seq

    assumes:
    all sequences are same length
    no sequence appears twice in seq_possibilities

    returns (best_hit, min_dist)
    * best_hit is closest sequence from seq_possibilities, or None if a tie
    * min_dist is the bitwise hamming distance between the query_seq and the
    best hit

    assuming default nt_to_bits of
    { "A":"11",  "C":"00", "T":"10", "G":"01"}, examples:
    cw1 = AAACCCGGGTTT (12 nucleotides)
    cw2 = AAACCCGGGTTA (hamming distance of 1 from cw1)
    cw3 = AAACCCGGGTTG (hamming distance of 1 from cw1)
    cw4 = AAACCCGGGTTC (hamming distance of 2 from cw1)
    cw4 = AAACCCGGGTCC (hamming distance of 4 from cw1)
    """
    if nt_to_bits is None:
        nt_to_bits = DEFAULT_NT_TO_BITS
    dists = []
    query_seq_bits = seq_to_bits(query_seq, nt_to_bits)
    for seq in seq_possibilities:
        possible_seq_bits = seq_to_bits(seq, nt_to_bits)
        dists.append(hamming_dist(query_seq_bits, possible_seq_bits))
    min_dist = min(dists)
    number_mins = dists.count(min_dist)
    if number_mins > 1:
        return None, min_dist
    else:
        best_hit = seq_possibilities[dists.index(min_dist)]
        return best_hit, min_dist


def hamming_dist(v1, v2):
    """ two binary arrays, equal len.  for bits, not nucleotide strings"""
    edits = (v1 != v2)
    return edits.sum()


def seq_to_bits(seq, nt_to_bits):
    """ e.g.: "AAG" -> array([1,1,1,1,0,1]), for a specific nt_to_bits
    """
    bitstring = ""
    for nt in seq:
        bitstring += nt_to_bits[nt]
    bits = numpy.array(map(int, bitstring))
    return bits
