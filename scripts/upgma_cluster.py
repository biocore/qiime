#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak+1@gmail.com"
__status__ = "Pre-release"
 

from qiime.util import parse_command_line_parameters
from optparse import make_option
from qiime.hierarchical_cluster import single_file_upgma, multiple_file_upgma
import os

script_description = """relate samples with UPGMA (resulting in a tree), using a distance matrix."""

script_usage = """
 single_file_example:
python %prog -i beta_unweighted_unifrac.txt -o sample_cluster.tre
creates file sample_cluster.tre, newick format of upgma clustering based on
distance matrix in beta_unweighted_unifrac.txt

or batch example: 
python %prog -i TEST/rare_unifrac -o TEST/rare_unifrac_upgma
processes every file in rare_unifrac, and creates a file "upgma_" + 
"base_fname.tre"
in rare_unifrac_upgma folder
-o is mandatory here, created if doesn't exist"""

required_options = [\
 # Example required option
 # make_option('-i','--input_dir',help='the input directory'),\
make_option('-i', '--input_path',
     help='input path.  directory for batch processing, '+\
      'filename for single file operation'),\
      

make_option('-o', '--output_path',
            help='output path. directory for batch processing, '+\
             'filename for single file operation'),\
]

optional_options = [\
 # Example optional option
 #make_option('-o','--output_dir',help='the output directory [default: %default]'),\

]




def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)
    if os.path.isdir(opts.input_path):
        multiple_file_upgma(opts.input_path,opts.output_path)
    elif os.path.isfile(opts.input_path):
        single_file_upgma(opts.input_path, opts.output_path)
    else:
        print("io error, check input file path")
        exit(1)

if __name__ == "__main__":
    main()