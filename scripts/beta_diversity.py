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
from qiime.beta_diversity import (single_file_beta, multiple_file_beta,
list_known_metrics)
import os

script_description = """Description:
calculate beta diversity (pairwise sample dissimilarity) on one or many
otu tables"""

script_usage = """
 Show available beta diversity metrics:
 python %prog -s
 
 
 Calculate two beta diversity measures on otu_table.txt:
 python %prog -i otu_table.txt -m bray_curtis,unweighted_unifrac -o outdir -t repr_set.tre

this creates two files: outdir/bray_curtis_otu_table.txt etc.

 Batch example: 
python %prog -i beta_rare_dir (a directory) -m bray_curtis,unweighted_unifrac -o outdir -t repr_set.tre
processes every file in beta_rare_dir, and creates a file "metric_" + infilename
in results folder"""

required_options = [\
 # Example required option
 #make_option('-i','--input_dir',help='the input directory'),\
]

optional_options = [\
 # Example optional option
 #make_option('-o','--output_dir',help='the output directory [default: %default]'),\
 make_option('-i', '--input_path',
     help='input path: otu table, or dir of otu tables for batch mode'),
     
 make_option('-o', '--output_dir',
     help="output directory, will be created if doesn't exist"),

 make_option('-m', '--metrics',
     help='metrics to use, comma delimited if >1 metric, '+\
         'no spaces'),
     
 make_option('-s', '--show_metrics', action='store_true', 
     help='show available beta diversity metrics and quit'),

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
    
    almost_required_options = ['input_path','output_dir','metrics']
    for option in almost_required_options:
        if getattr(opts,option) == None:
            option_parser.error('Required option --%s omitted.' % option)
    
    if opts.output_dir.endswith('.txt'):
        stderr.write('output must be a directory, files will be named'+\
          ' automatically.  And we refuse to make .txt directories\n')
        exit(1)
    
    try: 
        os.makedirs(opts.output_dir)
    except OSError:
        pass # hopefully dir already exists 
    if os.path.isdir(opts.input_path):
        multiple_file_beta(opts.input_path, opts.output_dir, opts.metrics, 
            opts.tree_path)
    elif os.path.isfile(opts.input_path):
        single_file_beta(opts.input_path, opts.metrics, opts.tree_path, 
          opts.output_dir)
    else:
      print("io error, input path not valid.  Does it exist?")
      exit(1)

if __name__ == "__main__":
    main()