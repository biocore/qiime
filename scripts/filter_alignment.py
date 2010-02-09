#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Doug Wendel"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso", "Justin Kuczynski", "Dan Knights", "Doug Wendel"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Doug Wendel"
__email__ = "wendel@colorado.edu"
__status__ = "Pre-release"
 

from qiime.filter_alignment import apply_lane_mask_and_gap_filter
from qiime.util import parse_command_line_parameters
from optparse import make_option
from numpy import nonzero, array, fromstring, repeat, bitwise_or, uint8, zeros
from random import sample
from cogent.parse.fasta import MinimalFastaParser
from cogent.core.alignment import eps
from optparse import OptionParser
from sys import stdout
from string import lowercase
from os.path import split, exists, splitext
from os import mkdir, remove
from qiime.util import load_qiime_config


script_description = """%prog is intended to be applied as a pre-filter for building trees from 
alignments. It removes positions which are highly variable based on a 
lane mask file (a string of 1s and 0s which is the length of the alignment
and where 1s indicate positions to keep; e.g., lanemask_in_1s_and_0s.txt from
greengenes), and positions which are 100% gap characters. Uses a minimal
fasta parser, safe to use on large sequence collections."""

script_usage = """
# filter 1.fasta using the lanemask in lm.txt, but filtering no gaps (b/c
# --allowed_gap_frac=1.0, meaning positions can be up to 100% gap); output
# written to ./1_pfiltered.fasta
filter_alignment.py -i 1.fasta -g 1.0 -m lm.txt

# filter 1.fasta using the lanemask in lm.txt and filter positions which are
# 100% gap (default -g behavior); output written to ./1_pfiltered.fasta
filter_alignment.py -i 1.fasta -o ./ -m lm.txt

# filter 1.fasta positions which are 100% gap (default -g behavior) but no lane mask
# filtering (because no lane mask file provided with -l); output written to 
# ./1_pfiltered.fasta
filter_alignment.py -i 1.fasta -o ./"""

required_options = [\
    make_option('-i','--input_fasta_file',action='store',\
         type='string',help='the input directory ')
]

qiime_config = load_qiime_config()
optional_options = [\
    make_option('-o','--output_dir',action='store',\
        type='string',help='the output directory '+\
        '[default: %default]',default='.'),\
   make_option('-m','--lane_mask_fp',action='store',\
        type='string',
        default=qiime_config['template_alignment_lanemask_fp'],
        help='path to lanemask file [default: %default]'),\
   make_option('-s','--suppress_lane_mask_filter',action='store_true',\
        help='suppress lane mask filtering (necessary to turn off '+\
        'lane-mask-based filtering when a qiime_config default is  '+\
        'provided for --lane_mask_fp) [default: %default]',default=False),\
   make_option('-g','--allowed_gap_frac',action='store',\
        type='float',help='gap filter threshold, ' +\
        'filters positions which are gaps in > allowed_gap_frac '+\
        'of the sequences [default: %default]',
        default=1.-eps)
]


def main():
    option_parser, opts, args = parse_command_line_parameters(
        script_description=script_description,
        script_usage=script_usage,
        version=__version__,
        required_options=required_options,
        optional_options=optional_options)
      
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


if __name__ == "__main__":
    main()
