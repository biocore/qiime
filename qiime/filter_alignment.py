#!/usr/bin/env python
from numpy import nonzero, array, fromstring, repeat, bitwise_or, uint8, zeros
from random import sample
from cogent.parse.fasta import MinimalFastaParser
from cogent.core.alignment import eps
from optparse import OptionParser
from sys import stdout
from string import lowercase
from os.path import split, exists, splitext
from os import mkdir, remove
from qiime.util import qiime_config

__author__ = "Dan Knights"
__copyright__ = "Copyright 2009, the PyCogent Project"
__credits__ = ["Greg Caporaso", "Justin Kuczynski", "Dan Knights"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Dan Knights"
__email__ = "danknights@gmail.com"
__status__ = "Prototype"

"""Contains code for filtering alignments before building trees from them
"""

usage_string = """usage: %prog [options] {-i INPUT_FASTA_FILE}

Description:
 %prog is intended to be applied as a pre-filter for building trees from 
 alignments. It removes positions which are highly variable based on a 
 lane mask file (a string of 1s and 0s which is the length of the alignment
 and where 1s indicate positions to keep; e.g., lanemask_in_1s_and_0s.txt from
 greengenes), and positions which are 100% gap characters. Uses a minimal
 fasta parser, safe to use on large sequence collections.

Example usage:

 # Get detailed usage information
 python filter_alignment.py -h

 # filter 1.fasta using the lanemask in lm.txt, but filtering no gaps (b/c
 # --allowed_gap_frac=1.0, meaning positions can be up to 100% gap); output
 # written to ./1_pfiltered.fasta
 python filter_alignment.py -i 1.fasta -g 1.0 -m lm.txt
 
 # filter 1.fasta using the lanemask in lm.txt and filter positions which are
 # 100% gap (default -g behavior); output written to ./1_pfiltered.fasta
 python filter_alignment.py -i 1.fasta -o ./ -m lm.txt
 
 # filter 1.fasta positions which are 100% gap (default -g behavior) but no lane mask
 # filtering (because no lane mask file provided with -l); output written to 
 # ./1_pfiltered.fasta
 python filter_alignment.py -i 1.fasta -o ./"""


def mask_to_positions(maskstring):
    """Converts lanemask binary string to array of valid indices."""
    return nonzero(array(map(int, maskstring)))[0]

def get_masked_string(s, p):
    """Extracts valid positions in string s using index array p."""
    return (fromstring(s,dtype=uint8))[p].tostring()

def find_gaps(s, gapcode=45):
    """Returns index array indicating locations of gaps ('-') in string s"""
    return nonzero(fromstring(s,dtype=uint8) == gapcode)

def apply_lane_mask(fastalines, lane_mask, verbose=False):
    """ Applies lanemask to fasta-formatted data, yielding filtered seqs."""
    return apply_lane_mask_and_gap_filter(fastalines, lane_mask, allowed_gap_frac=1, verbose=False)

def apply_gap_filter(fastalines, allowed_gap_frac=1-eps, verbose=False):
    """ Applies gap filter to fasta-formatted data, yielding filtered seqs."""
    return apply_lane_mask_and_gap_filter(fastalines, None, allowed_gap_frac=allowed_gap_frac, verbose=False)

def apply_lane_mask_and_gap_filter(fastalines, lane_mask, allowed_gap_frac=1-eps, verbose=False):
    """Applies lanemask and gap filter to fasta file with minimal parser, yielding filtered seqs.
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
        raise IOError, "Can't open temporary file for writing: %s" % tmpfilename

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

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = usage_string
    version = 'Version: %prog 0.1'
    parser = OptionParser(usage=usage, version=version)

    parser.add_option('-i','--input_fasta_file',action='store',\
           type='string',help='the input directory '+\
           '[REQUIRED]')
    parser.add_option('-o','--output_dir',action='store',\
           type='string',help='the output directory '+\
           '[default: %default]',default='.')
    parser.add_option('-m','--lane_mask_fp',action='store',\
            type='string',
            default=qiime_config['template_alignment_lanemask_fp'],
            help='path to lanemask file [default: %default]')
    parser.add_option('-s','--suppress_lane_mask_filter',action='store_true',\
            help='suppress lane mask filtering (necessary to turn off '+\
             'lane-mask-based filtering when a qiime_config default is  '+\
             'provided for --lane_mask_fp) [default: %default]',default=False)
    parser.add_option('-g','--allowed_gap_frac',action='store',\
            type='float',help='gap filter threshold, ' +\
            'filters positions which are gaps in > allowed_gap_frac '+\
            'of the sequences [default: %default]',
            default=1.-eps)
    parser.add_option('--verbose',action='store_true',\
            help='print all information [default: %default]')

    opts,args = parser.parse_args()

    required_options = ['input_fasta_file']
    
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 

    return opts,args
    
if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    
    # build the output filepath and open it any problems can be caught before starting
    # the work
    try:
        mkdir(opts.output_dir)
    except OSError:
        pass
    input_dir, input_filename = split(opts.input_fasta_file)
    input_basename, ext = splitext(input_filename)
    output_fp = '%s/%s_pfiltered.fasta' % (opts.output_dir,input_basename)
    
    try:
        outfile = open(output_fp,'w')
    except IOError:
        raise IOError, "Can't open output_filepath for writing: %s" % output_filepath

    
    # read the lane_mask, if one was provided
    if opts.verbose: print "Reading lane mask..."
    if opts.lane_mask_fp and not opts.suppress_lane_mask_filter:
        lane_mask = open(opts.lane_mask_fp).read().strip()
    else:
        lane_mask = None
    # open the input and output files        
    infile = open(opts.input_fasta_file,'U')

    # apply the lanemask/gap removal
    for result in apply_lane_mask_and_gap_filter(infile, lane_mask, opts.allowed_gap_frac, verbose=opts.verbose):
        outfile.write(result)
    infile.close()
    outfile.close()

