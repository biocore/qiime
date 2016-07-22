#!/usr/bin/env python

from random import sample
from sys import stdout
from string import lowercase
from os.path import split, exists, splitext
from os import mkdir, remove
from collections import defaultdict
from itertools import izip

import numpy as np
from numpy import (nonzero, array, fromstring, repeat, bitwise_or,
    uint8, zeros, arange, finfo, mean, std)

from skbio.alignment import Alignment
from skbio.sequence import DNA
from skbio.parse.sequences import parse_fasta


__author__ = "Dan Knights"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Justin Kuczynski", "Dan Knights"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

"""Contains code for filtering alignments before building trees from them
"""


def mask_to_positions(maskstring):
    """Converts lanemask binary string to array of valid indices."""
    return nonzero(array(map(int, maskstring)))[0]


def get_masked_string(s, p):
    """Extracts valid positions in string s using index array p."""
    return (fromstring(s, dtype=uint8))[p].tostring()


def apply_lane_mask_and_gap_filter(fastalines, mask,
                                   allowed_gap_frac=1.-1e-6,
                                   entropy_threshold=None):
    """Applies a mask and gap filter to fasta file, yielding filtered seqs."""

    # load the alignment
    aln = Alignment.from_fasta_records(parse_fasta(fastalines), DNA)

    # build the entropy mask
    if mask is None and entropy_threshold is None:
        if allowed_gap_frac < 1:
            aln = aln.omit_gap_positions(allowed_gap_frac)
    elif mask is not None and entropy_threshold is not None:
        raise ValueError('only mask or entropy threshold can be provided.')
    elif mask is not None:
        # a pre-computed mask (e.g., Lane mask) was provided, so apply that
        # and then remove highly gapped positions (gap positions have to be
        # removed after the mask-based filtering so that the positions in the
        # mask correspond with the positions in the alignment at the time of
        # filtering)
        entropy_mask = mask_to_positions(mask)
        aln = aln.subalignment(positions_to_keep=entropy_mask)
        if allowed_gap_frac < 1:
            aln = aln.omit_gap_positions(allowed_gap_frac)
    elif entropy_threshold is not None:
        # a mask is being computed on the fly to filter the entropy_threshold
        # most entropic positions. if highly gapped positions are being omitted
        # those are filtered first, so the entropy scores for those positions
        # aren't included when determining the entropy threshold (since the
        # positions that are mostly gaps will be counted as a lot of low
        # entropy positions)
        if not (0 <= entropy_threshold <= 1):
            raise ValueError('entropy_threshold needs to be between 0 and 1'
                             ' (inclusive)')
        if allowed_gap_frac < 1:
            aln = aln.omit_gap_positions(allowed_gap_frac)
        entropy_mask = generate_lane_mask(aln, entropy_threshold)
        entropy_mask = mask_to_positions(entropy_mask)
        aln = aln.subalignment(positions_to_keep=entropy_mask)
    else:
        # it shouldn't be possible to get here
        raise ValueError("Can't resolve parameters for "
                         "apply_lane_mask_and_gap_filter.")

    if aln.sequence_length() == 0:
        raise ValueError("Positional filtering resulted in removal of all "
                         "alignment positions.")
    for seq in aln:
        yield ">%s\n" % seq.id
        yield "%s\n" % seq

def remove_outliers(seqs, num_stds, fraction_seqs_for_stats=.95):
    """ remove sequences very different from the majority consensus

    given aligned sequences, will:
     1. calculate a majority consensus (most common symbol at each position
        of the alignment);
     2. compute the mean/std edit distance of each seq to the consensus;
     3. discard sequences whose edit dist is greater than the cutoff, which is
        defined as being `num_stds` greater than the mean.

    """
    # load the alignment and compute the consensus sequence
    aln = Alignment.from_fasta_records(parse_fasta(seqs), DNA)
    consensus_seq = aln.majority_consensus()
    # compute the hamming distance between all sequences in the alignment
    # and the consensus sequence
    dists_to_consensus = [s.distance(consensus_seq) for s in aln]
    # compute the average and standard deviation distance from the consensus
    average_distance = mean(dists_to_consensus)
    std_distance = std(dists_to_consensus)
    # compute the distance cutoff
    dist_cutoff = average_distance + num_stds * std_distance
    # for all sequences, determine if they're distance to the consensus
    # is less then or equal to the cutoff distance. if so, add the sequence's
    # identifier to the list of sequence identifiers to keep
    seqs_to_keep = []
    for seq_id, dist_to_consensus in izip(aln.ids(), dists_to_consensus):
        if dist_to_consensus <= dist_cutoff:
            seqs_to_keep.append(seq_id)
    # filter the alignment to only keep the sequences identified in the step
    # above
    filtered_aln = aln.subalignment(seqs_to_keep=seqs_to_keep)
    # and return the filtered alignment
    return filtered_aln


def generate_lane_mask(aln, entropy_threshold):
    """ Generates lane mask dynamically by calculating base frequencies

    aln: skbio.Alignment object
    entropy_threshold:  float value that designates the percentage of entropic
     positions to be removed, i.e., 0.10 means the 10% most entropic positions
     are removed.

    """
    if entropy_threshold == 0.0:
        # special case: keep all positions if filtering the 0% most entropic
        # positions
        result = ["1"] * aln.sequence_length()
    else:
        entropies = aln.position_entropies(nan_on_non_standard_chars=False)
        entropy_cutoff = np.percentile(entropies,
                                       (1 - entropy_threshold) *100.,
                                       interpolation='nearest')
        result = []
        for entropy in entropies:
            if entropy >= entropy_cutoff:
                result.append("0")
            else:
                result.append("1")

    return "".join(result)
