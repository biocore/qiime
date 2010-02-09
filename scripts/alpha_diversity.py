#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Justin Kuczynski", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Pre-release"
 

from qiime.util import parse_command_line_parameters
from optparse import make_option
from qiime.alpha_diversity import (single_file_alpha, multiple_file_alpha,
list_known_metrics)
import os

script_description = """
calculate alpha diversity on each sample in an otu table, using a variety of
alpha divirsity metrics"""

script_usage = """
 see what alpha diversity metrics are available:
 python %prog -s
 
 single analysis: 
 python %prog -i otu_table.txt -m observed_species,chao1,PD_whole_tree -o alphadiv.txt -t repr_set.tre
 creates one file containing the the observed species, chao1 diversity,
 and PD for each sample in otu_table.txt

 or batch example: 
 python %prog -i alpha_rare (a directory) -m observed_species,chao1,PD_whole_tree -t repr_set.tre -o batch_alpha
processes every file in alpha_rare, and for each file, creates a file 
"alpha_" + fname in batch_alpha folder"""

required_options = [\
 # Example required option
 #make_option('-i','--input_dir',help='the input directory'),\

]

optional_options = [\
 # Example optional option
 #make_option('-o','--output_dir',help='the output directory [default: %default]'),\
 make_option('-i', '--input_path',
     help='input path.  directory for batch processing, '+\
      'filename for single file operation'),
     
 make_option('-o', '--output_path',
     help='output path. directory for batch processing, '+\
      'filename for single file operation'),

 make_option('-m', '--metrics',
     help='metrics to use, comma delimited'),
     
 make_option('-s', '--show_metrics', action='store_true', 
     dest="show_metrics",
     help='show available alpha diversity metrics and quit'),

 make_option('-t', '--tree_path', default=None,
     help='path to newick tree file, required for phylogenetic metrics'+\
     ' [default: %default]')
]


def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)

    if opts.show_metrics:
        print("Known metrics are: %s\n" \
              % (', '.join(list_known_metrics()),))
        exit(0)
    almost_required_options = ['input_path','output_path','metrics']
    for option in almost_required_options:
        if getattr(opts,option) == None:
            option_parser.error('Required option --%s omitted.' % option)
    
    if os.path.isdir(opts.input_path):
      multiple_file_alpha(opts.input_path, opts.output_path, opts.metrics, 
        opts.tree_path)
    elif os.path.isfile(opts.input_path):
      try:
          f = open(opts.output_path, 'w')
          f.close()
      except IOError:
          print("ioerror, couldn't create output file")
          exit(1)
      single_file_alpha(opts.input_path, opts.metrics, 
          opts.output_path, opts.tree_path)
    else:
      print("io error, input path not valid. does it exist?")
      exit(1)

if __name__ == "__main__":
    main()