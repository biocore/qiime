#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Justin Kuczynski","Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Pre-release"
 

from qiime.util import parse_command_line_parameters
from optparse import make_option
from cogent import LoadTree
from qiime.make_bootstrapped_tree import write_pdf_bootstrap_tree
from qiime.parse import parse_bootstrap_support


script_description = """
Takes a tree and bootstrap support file, then writes a pdf, colored by bootstrap \
support.
"""

script_usage = """
Example usage:

make_bootstrapped_tree.py -m master_tree.tre -s jackknife_support.txt -o jackknife_samples.pdf
"""

required_options = [\
 # Example required option
 #make_option('-i','--input_dir',help='the input directory'),\
 make_option('-m', '--master_tree', help='This is the path to the master tree'),
 make_option('-s', '--support', help='This is the path to the bootstrap \
support file'),
 make_option('-o', '--output_file', help="This is the filename where the output \
should be written.  If the suffix does not contain .pdf, then it will be \
attached")
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

    tree = LoadTree(opts.master_tree)
    support_file = open(opts.support)
    bootstraps = parse_bootstrap_support(support_file)
    support_file.close()
    write_pdf_bootstrap_tree(tree, opts.output_file, bootstraps)

if __name__ == "__main__":
    main()