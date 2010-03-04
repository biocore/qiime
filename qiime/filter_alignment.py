#!/usr/bin/env python
from numpy import nonzero, array, fromstring, repeat, bitwise_or, uint8, zeros
from random import sample
from cogent.parse.fasta import MinimalFastaParser
from cogent.core.alignment import eps
from sys import stdout
from string import lowercase
from os.path import split, exists, splitext
from os import mkdir, remove

__author__ = "Dan Knights"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso", "Justin Kuczynski", "Dan Knights"]
__license__ = "GPL"
__version__ = "0.92"
__maintainer__ = "Dan Knights"
__email__ = "danknights@gmail.com"
__status__ = "Release"

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
    allowed_gap_frac=1-eps, verbose=False):
    """Applies lanemask and gap filter to fasta file, yielding filtered seqs.
    """
    
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
    tmpfilename = "/tmp/"+"".join(sample(lowercase, 20)) + ".tmp"
    try:
        tmpfile = open(tmpfilename,'w')
    except IOError:
        raise IOError, "Can't open temporary file for writing: %s" %\
          tmpfilename

    # the number of gaps seen in each position (length may be unknown here)
    gapcounts = None

    # First pass: apply filter, and track gaps
    if verbose: print "First pass: applying lanemask..."
    seq_count = 0
    for k, v in MinimalFastaParser(fastalines):
        seq_count += 1
        # print progress in verbose mode
        if verbose and (seq_count % 100) == 0: status(seq_count)

        # apply lanemask if there is one
        if lane_mask:
            masked = get_masked_string(v,p)
        else:
            masked = v

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
    tmpfile = open(tmpfilename,'r')


    # if we're not removing gaps, we're done; yield the temp file contents
    if allowed_gap_frac == 1:
        for line in tmpfile:
            yield line
            
    # else we are removing gaps; do second pass
    else:

        # convert gapcounts to true/false mask
        gapcounts = (gapcounts / float(seq_count) ) <= allowed_gap_frac
        
        # Second pass: remove all-gap positions
        if verbose: print "Second pass: remove all-gap positions..."
        seq_count = 0
        for k, v in MinimalFastaParser(tmpfile):
            seq_count += 1
            # print progress in verbose mode
            if verbose and (seq_count % 100) == 0: status(seq_count)
            
            masked = get_masked_string(v,gapcounts)
            yield '>%s\n' % (k)
            yield '%s\n' % (masked)
        if verbose: print 

    # delete temporary file
    tmpfile.close()
    remove(tmpfilename)


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
