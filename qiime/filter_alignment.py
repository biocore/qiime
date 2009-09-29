#!/usr/bin/env python

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2009, the PyCogent Project"
__credits__ = ["Greg Caporaso", "Justin Kuczynski"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Prototype"

"""Contains code filtering alignments before building trees from them
"""

from optparse import OptionParser
from os.path import split, exists, splitext
from os import mkdir
from cogent.core.alignment import eps, DenseAlignment
from cogent import LoadSeqs

usage_string = """usage: %prog [options] {-i INPUT_FASTA_FILE}

Description:
 %prog is intended to be applied as a pre-filter for building trees from 
 alignments. It removes positions which are highly variable based on a 
 lane mask file (a string of 1s and 0s which is the length of the alignment
 and where 1s indicate positions to keep; e.g., lanemask_in_1s_and_0s.txt from
 greengenes), and positions which are 100% gap characters.

Example usage:

 # Get detailed usage information
 python filter_alignment.py -h

 # filter 1.fasta using the lanemask in lm.txt, but filtering no gaps (b/c
 # --allowed_gap_frac=1.0, meaning positions can be up to 100% gap); output
 # written to ./1_filtered.fasta
 python filter_alignment.py -i 1.fasta -g 1.0 -l lm.txt
 
 # filter 1.fasta using the lanemask in lm.txt and filter positions which are
 # 100% gap (default -g behavior); output written to ./1_filtered.fasta
 python filter_alignment.py -i 1.fasta -o ./ -l lm.txt
 
 # filter 1.fasta positions which are 100% gap (default -g behavior) but no lane mask
 # filtering (because no lane mask file provided with -l); output written to 
 # ./1_filtered.fasta
 python filter_alignment.py -i 1.fasta -o ./ -l lm.txt"""

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = usage_string
    version = 'Version: %prog 0.1'
    parser = OptionParser(usage=usage, version=version)

    # An example string option
    parser.add_option('-i','--input_fasta_file',action='store',\
           type='string',help='the input directory '+\
           '[REQUIRED]')
    parser.add_option('-o','--output_dir',action='store',\
           type='string',help='the output directory '+\
           '[default: %default]',default='.')
    parser.add_option('-l','--lane_mask_fp',action='store',\
            type='string',help='path to lanemask file '+\
            '[default: %default]')
    parser.add_option('-g','--allowed_gap_frac',action='store',\
            type='float',help='gap filter threshold, ' +\
            'filters positions which are gaps in > allowed_gap_frac '+\
            'of the sequences [default: %default]',
            default=1.-eps)

    opts,args = parser.parse_args()

    required_options = ['input_fasta_file']
    
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 

    return opts,args

def apply_lane_mask(aln,lane_mask):
    """ Remove positions corresponding to zeros in lane_mask
    """
    zero_values = {}.fromkeys([0,'0'])
    positions_to_keep = \
     [i for i in range(len(lane_mask)) if lane_mask[i] not in zero_values]
    return aln.takePositions(positions_to_keep)
    
def apply_gap_filter(aln,allowed_gap_frac=1.-eps):
    """ Remove positions that contain gaps in > allowed_gap_frac of sequences
    """
    return aln.omitGapPositions(allowed_gap_frac=allowed_gap_frac)
    
def apply_lane_mask_and_gap_filters(aln,lane_mask,allowed_gap_frac=1.-eps):
    
    # must apply the lane mask first as it relies on the original positions
    if lane_mask:
        result = apply_lane_mask(aln,lane_mask)
    else:
        # if no lane mask, don't do anything to the alignment
        result = aln
    
    # apply the gap filter
    result = apply_gap_filter(result,allowed_gap_frac)
    
    # return the result
    return result
    
if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    
    # create local copies of the options
    lane_mask_fp = opts.lane_mask_fp
    allowed_gap_frac = opts.allowed_gap_frac
    input_fasta_file = opts.input_fasta_file
    output_dir = opts.output_dir
    
    # build the output filepath and open it any problems can be caught before starting
    # the work
    try:
        mkdir(output_dir)
    except OSError:
        pass
    input_dir, input_filename = split(input_fasta_file)
    input_basename, ext = splitext(input_filename)
    output_filepath = '%s/%s_filtered.fasta' % (output_dir,input_basename)
    
    try:
        output_file = open(output_filepath,'w')
    except IOError:
        raise IOError, "Can't open output_filepath for writing: %s" % output_filepath
    
    # load the input alignment
    aln = LoadSeqs(input_fasta_file,aligned=DenseAlignment)
    
    # read the lane_mask, if one was provided
    if lane_mask_fp:
        lane_mask = open(lane_mask_fp).read().strip()
        if len(aln) != len(lane_mask):
            raise ValueError,\
             "Lane mask (%s) and alignment (%s) differ in number of positions." %\
             (lane_mask_fp, input_fasta_file)
    else:
        lane_mask = None
        
    # apply the filters
    filtered_aln = apply_lane_mask_and_gap_filters(aln,lane_mask,allowed_gap_frac)
    
    # write the alignment to the output_file
    output_file.write(filtered_aln.toFasta())
    output_file.close()
    