#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Doug Wendel"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso", "Justin Kuczynski", "Dan Knights", \
    "Doug Wendel"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Doug Wendel"
__email__ = "wendel@colorado.edu"
__status__ = "Pre-release"
 

from qiime.filter_alignment import apply_lane_mask_and_gap_filter
from qiime.util import parse_command_line_parameters
from optparse import make_option
from cogent.core.alignment import eps
from os.path import split, exists, splitext
from os import mkdir, remove
from qiime.util import load_qiime_config

script_info={}
script_info['brief_description']="""Filter sequence alignment by removing highly variable regions"""
script_info['script_description']="""This script should be applied to generate a useful tree when aligning against a template alignment (e.g., with PyNAST). This script will remove positions which are gaps in every sequence (common for PyNAST, as typical sequences cover only 200-400 bases, and they are being aligned against the full 16S gene). Additionally, the user can supply a lanemask file, that defines which positions should included when building the tree, and which should be ignored. Typically, this will differentiate between non-conserved positions, which are uninformative for tree building, and conserved positions which are informative for tree building. FILTERING ALIGNMENTS WHICH WERE BUILD WITH PYNAST AGAINST THE GREENGENES CORE SET ALIGNMENT SHOULD BE CONSIDERED AN ESSENTIAL STEP."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Examples:""","""As a simple example of this script, the user can use the following command, which consists of an input FASTA file (i.e. resulting file from align_seqs.py), lanemask template file and the output directory "filtered_alignment/":""","""%prog -i repr_set_seqs_aligned.fna -m lanemask_template -o filtered_alignment/"""))
script_info['script_usage'].append(("","""Alternatively, if the user would like to use a different gap fraction threshold ("-g"), they can use the following command:""","""%prog -i repr_set_seqs_aligned.fna -m lanemask_template -o filtered_alignment/ -g 0.95"""))
script_info['output_description']="""The output of filter_alignment.py consists of a single FASTA file, which ends with "pfiltered.fasta", where the "p" stands for positional filtering of the columns."""
script_info['required_options']= [\
    make_option('-i','--input_fasta_file',action='store',\
         type='string',help='the input directory ')
]

qiime_config = load_qiime_config()

script_info['optional_options']= [\
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
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
      
    # build the output filepath and open it any problems can be caught 
    # before starting the work
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
        raise IOError, "Can't open output_filepath for writing: %s" \
         % output_filepath


    # read the lane_mask, if one was provided
    if opts.verbose: print "Reading lane mask..."
    if opts.lane_mask_fp and not opts.suppress_lane_mask_filter:
        lane_mask = open(opts.lane_mask_fp).read().strip()
    else:
        lane_mask = None
    # open the input and output files        
    infile = open(opts.input_fasta_file,'U')

    # apply the lanemask/gap removal
    for result in apply_lane_mask_and_gap_filter(infile, lane_mask, \
     opts.allowed_gap_frac, verbose=opts.verbose):
        outfile.write(result)
    infile.close()
    outfile.close()


if __name__ == "__main__":
    main()
