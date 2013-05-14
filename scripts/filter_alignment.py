#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Doug Wendel"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Justin Kuczynski", "Dan Knights", \
    "Doug Wendel", "William Walters", "John Chase"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Doug Wendel"
__email__ = "wendel@colorado.edu"
__status__ = "Development"
 

from os.path import split, exists, splitext, getsize
from os import mkdir, remove

from cogent.core.alignment import eps

from qiime.util import load_qiime_config
from qiime.filter_alignment import apply_lane_mask_and_gap_filter, \
  remove_outliers, generate_lane_mask
from qiime.util import parse_command_line_parameters
from qiime.util import make_option

script_info={}
script_info['brief_description']="""Filter sequence alignment by removing highly variable regions"""
script_info['script_description']="""This script should be applied to generate a useful tree when aligning against a template alignment (e.g., with PyNAST). This script will remove positions which are gaps in every sequence (common for PyNAST, as typical sequences cover only 200-400 bases, and they are being aligned against the full 16S gene). Additionally, the user can supply a lanemask file, that defines which positions should included when building the tree, and which should be ignored. Typically, this will differentiate between non-conserved positions, which are uninformative for tree building, and conserved positions which are informative for tree building. FILTERING ALIGNMENTS WHICH WERE BUILD WITH PYNAST AGAINST THE GREENGENES CORE SET ALIGNMENT SHOULD BE CONSIDERED AN ESSENTIAL STEP."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Examples:""","""As a simple example of this script, the user can use the following command, which consists of an input FASTA file (i.e. resulting file from align_seqs.py), lanemask template file and the output directory "filtered_alignment/":""","""%prog -i seqs_rep_set_aligned.fasta -m lanemask_in_1s_and_0s -o filtered_alignment/"""))

script_info['script_usage'].append(("","""Alternatively, if the user would like to use a different gap fraction threshold ("-g"), they can use the following command:""","""%prog -i seqs_rep_set_aligned.fasta -m lanemask_in_1s_and_0s -o filtered_alignment/ -g 0.95"""))

script_info['output_description']="""The output of filter_alignment.py consists of a single FASTA file, which ends with "pfiltered.fasta", where the "p" stands for positional filtering of the columns."""

script_info['required_options']= [\
    make_option('-i','--input_fasta_file',action='store',\
         type='existing_filepath',help='the input directory ')
]

qiime_config = load_qiime_config()

script_info['optional_options']= [\
    make_option('-o','--output_dir',action='store',\
        type='new_dirpath',help='the output directory '+\
        '[default: %default]',default='.'),\
    make_option('-m','--lane_mask_fp',action='store',\
        type='existing_filepath',
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
        default=1.-eps),
    make_option('-r','--remove_outliers',action='store_true',\
        help='remove seqs very dissimilar to the alignment consensus' +\
        ' (see --threshold).  [default: %default]',
        default=False),
    make_option('-t', '--threshold',action='store',\
        type='float',help='with -r, remove seqs whose dissimilarity to the ' +\
        'consensus sequence is approximately > x standard devaitions above '+\
        'the mean of the sequences [default: %default]',
        default=3.0),
    make_option('-e', '--entropy_threshold', action='store',\
        type='float',help='Sets percent threshold for removing base '+\
        'positions with the highest entropy.  For example, if 0.10 were '+\
        'specified, the top 10% most entropic base positions would be '+\
        'filtered.  If this value is used, any lane mask supplied will be '+\
        'ignored.  Entropy filtered occurs after gap filtering.  '+\
        '  [default: %default]', default = None)
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
    
    if getsize(opts.input_fasta_file) == 0:
        raise ValueError, ("An empty fasta file was provided. "
        "Did the alignment complete sucessfully? "
        "Did PyNAST discard all sequences due to too-stringent minimum length "
        "or minimum percent ID settings?")
    
    output_fp = '%s/%s_pfiltered.fasta' % (opts.output_dir,input_basename)

    try:
        outfile = open(output_fp,'w')
    except IOError:
        raise IOError, "Can't open output_filepath for writing: %s" \
         % output_filepath
         


    
    
        
    if opts.lane_mask_fp and not opts.suppress_lane_mask_filter and not\
     opts.entropy_threshold:
        # read the lane_mask, if one was provided
        if opts.verbose: print "Reading lane mask..."
        lane_mask = open(opts.lane_mask_fp).read().strip()
    else:
        lane_mask = None
    # open the input and output files        
    infile = open(opts.input_fasta_file,'U')

    
    if opts.remove_outliers:
        # apply the lanemask/gap removal, then remove outliers

        seq_gen = apply_lane_mask_and_gap_filter(infile, lane_mask,
          opts.allowed_gap_frac, verbose=opts.verbose,
          entropy_threshold=opts.entropy_threshold)
          
        filtered_aln = remove_outliers(seq_gen, opts.threshold)
        for seq in filtered_aln.Seqs:
            outfile.write(seq.toFasta())
            outfile.write('\n')
   
    else:
        # just apply the lanemask/gap removal
        for result in apply_lane_mask_and_gap_filter(infile, lane_mask, \
         opts.allowed_gap_frac, verbose=opts.verbose,
         entropy_threshold=opts.entropy_threshold):
            outfile.write(result)
    infile.close()
    outfile.close()


if __name__ == "__main__":
    main()
