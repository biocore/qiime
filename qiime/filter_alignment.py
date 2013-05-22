#!/usr/bin/env python

from random import sample
from sys import stdout
from string import lowercase
from os.path import split, exists, splitext
from os import mkdir, remove
from collections import defaultdict

from cogent.util.unit_test import TestCase, main
from numpy import nonzero, array, fromstring, repeat, bitwise_or, uint8, zeros,\
 arange
import numpy
from cogent import LoadSeqs, DNA
from cogent.core.alignment import DenseAlignment, eps
from cogent.parse.fasta import MinimalFastaParser
from cogent.core.sequence import ModelDnaSequence
from cogent.core.profile import Profile

from qiime.util import get_tmp_filename


__author__ = "Dan Knights"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Justin Kuczynski", "Dan Knights"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Dan Knights"
__email__ = "danknights@gmail.com"
__status__ = "Development"

"""Contains code for filtering alignments before building trees from them
"""


def mask_to_positions(maskstring):
    """Converts lanemask binary string to array of valid indices."""
    return nonzero(array(map(int, maskstring)))[0]

def get_masked_string(s, p):
    """Extracts valid positions in string s using index array p."""
    return (fromstring(s,dtype=uint8))[p].tostring()

def find_gaps(s, gapcode=45):
    """Returns index array indicating locations of gaps ('-') in string s
    """
    return nonzero(fromstring(s,dtype=uint8) == gapcode)

def apply_lane_mask(fastalines, lane_mask, verbose=False):
    """ Applies lanemask to fasta-formatted data, yielding filtered seqs.
    """
    return apply_lane_mask_and_gap_filter(fastalines, lane_mask,\
     allowed_gap_frac=1, verbose=False)

def apply_gap_filter(fastalines, allowed_gap_frac=1-eps, verbose=False):
    """ Applies gap filter to fasta-formatted data, yielding filtered seqs.
    """
    return apply_lane_mask_and_gap_filter(fastalines, None, \
     allowed_gap_frac=allowed_gap_frac, verbose=False)

def apply_lane_mask_and_gap_filter(fastalines, lane_mask,\
    allowed_gap_frac=1-eps, verbose=False, entropy_threshold=None):
    """Applies lanemask and gap filter to fasta file, yielding filtered seqs.
    """

    if entropy_threshold:
        if entropy_threshold < 0 or entropy_threshold > 1:
            raise ValueError,('Entropy threshold parameter (-e) needs to be '+\
             'between 0 and 1')
    
    if lane_mask:
        # convert lane_mask to a numpy index array
        p = mask_to_positions(lane_mask)
        
        # special case: lanemask is all zeros
        if sum(p) == 0:
            for line in fastalines:
                if line.startswith(">"):
                    yield line + '\n'
                else:
                    yield '\n'
            return

    # random temporary file for first-pass results
    tmpfilename = get_tmp_filename()
    try:
        tmpfile = open(tmpfilename,'w')
    except IOError:
        raise IOError, "Can't open temporary file for writing: %s" %\
          tmpfilename

    # the number of gaps seen in each position (length may be unknown here)
    gapcounts = None

    # First pass: apply filter, and track gaps
    if verbose and lane_mask:
        print "Applying lanemask..."
    seq_count = 0
    for k, v in MinimalFastaParser(fastalines):
        seq_count += 1
        # print progress in verbose mode
        if verbose and (seq_count % 100) == 0: status(seq_count)

        # apply lanemask if there is one
        if lane_mask:
            masked = get_masked_string(v,p)
        else:
            masked = v.replace('.', '-')

        # initialize gapcount array to proper length
        if gapcounts == None:
            gapcounts = zeros(len(masked))

        # increment gap counts if requested
        if allowed_gap_frac < 1:
            gapcounts[find_gaps(masked)] += 1
                     
        # write masked sequence to temporary file
        tmpfile.write('>%s\n%s\n' % (k, masked))
    if verbose: print; print
    tmpfile.close()
    tmpfile = open(tmpfilename,'U')
    
    # random temporary file for second-pass results
    tmpfilename_gaps = get_tmp_filename()
    try:
        tmpfile_gaps = open(tmpfilename_gaps,'w')
    except IOError:
        raise IOError, "Can't open temporary file for writing: %s" %\
          tmpfilename_gaps


    # if we're not removing gaps, we're done; yield the temp file contents
    if allowed_gap_frac == 1:
        for line in tmpfile:
            yield line
            
    # else we are removing gaps; do second pass
    else:

        # convert gapcounts to true/false mask
        gapcounts = (gapcounts / float(seq_count) ) <= allowed_gap_frac
        
        # Second pass: remove all-gap positions
        if verbose: print "Remove all-gap positions..."
        seq_count = 0
        for k, v in MinimalFastaParser(tmpfile):
            seq_count += 1
            # print progress in verbose mode
            if verbose and (seq_count % 100) == 0: status(seq_count)
            
            masked = get_masked_string(v.replace('.','-'),gapcounts)
            tmpfile_gaps.write('>%s\n' % (k))
            tmpfile_gaps.write('%s\n' % (masked))
        if verbose: print
    
    tmpfile_gaps.close()
    tmpfile_gaps = open(tmpfilename_gaps, "U")
        
    # If no dynamic entropy calculations, return current values
    if not entropy_threshold:
        for line in tmpfile_gaps:
            yield line
    # Otherwise, filter out positions of post-gap filtered sequences
    else:
        if verbose:
            print "Generating lanemask..."
        lane_mask = generate_lane_mask(tmpfile_gaps, entropy_threshold)
        tmpfile_gaps.close()
        tmpfile_gaps = open(tmpfilename_gaps, "U")
        if verbose:
            print "Applying lanemask..."
        
        p = mask_to_positions(lane_mask)
        
        seq_count = 0
        for k, v in MinimalFastaParser(tmpfile_gaps):
            seq_count += 1
            # print progress in verbose mode
            if verbose and (seq_count % 100) == 0: status(seq_count)

            masked = get_masked_string(v,p)
            
            yield(">%s\n%s\n" % (k, masked))
            
        if verbose:
            print
        
    # delete temporary files
    tmpfile.close()
    tmpfile_gaps.close()
    remove(tmpfilename)
    remove(tmpfilename_gaps)

def remove_outliers(seqs, num_sigmas, fraction_seqs_for_stats=.95):
    """ remove sequences very different from the majority consensus
    
    given aligned seqs, will calculate a majority consensus (most common
    symbol at each position of the alignment), and average edit distance
    of each seq to that consensus.  any seq whose edit dist is > cutoff 
    (roughly seq_dist > num_sigmas * (average edit dist) ) is removed
    when calculating mean and stddev edit distance, only the best
    fraction_seqs_for_stats are used
    
    seqs must be compatible with DenseAlignment: 
    aln = DenseAlignment(data=seqs, MolType=DNA) is called
    """
    aln = DenseAlignment(data=seqs, MolType=DNA)
    cons = DenseAlignment(data=aln.majorityConsensus(), MolType=DNA)
    diff_mtx = cons.SeqData[:,0] != aln.SeqData
    
    # consider only a fraction of seqs for mean, std
    seq_diffs = diff_mtx.sum(1)
    num_to_consider = round(len(seq_diffs)*fraction_seqs_for_stats)
    seq_diffs_considered_sorted = \
        seq_diffs[seq_diffs.argsort()[:num_to_consider]]
    diff_cutoff = seq_diffs_considered_sorted.mean() + \
      num_sigmas*seq_diffs_considered_sorted.std()
    # mean + e.g.: 4 sigma
    seq_idxs_to_keep = numpy.arange(len(seq_diffs))[seq_diffs <= diff_cutoff]
    
    filtered_aln = aln.getSubAlignment(seq_idxs_to_keep)
    return filtered_aln
        

def status(message,dest=stdout,overwrite=True, max_len=100):
    """Writes a status message over the current line of stdout
    """
    message = str(message)
    message_len = max(len(message),max_len)
    if overwrite:
        dest.write('\b' * (message_len+2))
    dest.write(message[0:message_len])
    if not overwrite:
        dest.write('\n')
    dest.flush()
    
def generate_lane_mask(infile, entropy_threshold):
    """ Generates lane mask dynamically by calculating base frequencies
    
    infile: open file object for aligned fasta file
    entropy_threshold:  float value that designates the percentage of entropic
     positions to be removed, i.e., 0.10 means the 10% most entropic positions
     are removed.
     
    """
    
    
    base_freqs = freqs_from_aln_array(infile)
    uncertainty = base_freqs.columnUncertainty()
    uncertainty_sorted = list(uncertainty)
    uncertainty_sorted.sort()
    
    cutoff_index = int(round((len(uncertainty_sorted)-1) *\
     (1 - entropy_threshold)))
    
    max_uncertainty = uncertainty_sorted[cutoff_index]
    
    # This correction is for small datasets with a small possible number of
    # uncertainty values.
    highest_certainty = min(uncertainty_sorted)
    
    lane_mask = ""
    
    for base in uncertainty:
        if base >= max_uncertainty and base != highest_certainty:
            lane_mask += "0"
        else:
            lane_mask += "1"
    
    
    return lane_mask
    

   
def freqs_from_aln_array(seqs):
    """Returns per-position freqs from arbitrary size alignment.

    Warning: fails if all seqs aren't the same length.
    written by Rob Knight
    
    seqs = list of lines from aligned fasta file
    """
    result = None
    for label, seq in MinimalFastaParser(seqs):
        # Currently cogent does not support . characters for gaps, converting
        # to - characters for compatability.
        seq = ModelDnaSequence(seq.replace('.','-'))
        if result is None:
            result = zeros((len(seq.Alphabet), len(seq)),dtype=int)
            indices = arange(len(seq), dtype=int)
        result[seq._data,indices] += 1
    return Profile(result, seq.Alphabet)
