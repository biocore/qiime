#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski","Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Release"

from matplotlib import use
use('Agg',warn=False)
import matplotlib
from qiime.util import make_option
from qiime.parse import parse_newick
from qiime.util import parse_command_line_parameters
from qiime.make_bootstrapped_tree import write_pdf_bootstrap_tree
from qiime.parse import parse_bootstrap_support

script_info={}
script_info['brief_description']="""Make bootstrapped tree"""
script_info['script_description']="""This script takes a tree and bootstrap support file and creates a pdf, colored by bootstrap support."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""In this example, the user supplies a tree file and a text file containing the jackknife support information, which results in a pdf file:""","""make_bootstrapped_tree.py -m master_tree.tre -s jackknife_support.txt -o jackknife_samples.pdf"""))
script_info['output_description']="""The result of this script is a pdf file."""
script_info['required_options']=[\
 make_option('-m', '--master_tree', help='This is the path to the master tree'),
 make_option('-s', '--support', help='This is the path to the bootstrap \
support file'),
 make_option('-o', '--output_file', help="This is the filename where the output \
should be written.  If the suffix does not contain .pdf, then it will be \
attached")
]
script_info['optional_options']=[]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    tree = parse_newick(open(opts.master_tree,'U'))
    support_file = open(opts.support)
    bootstraps = parse_bootstrap_support(support_file)
    support_file.close()
    write_pdf_bootstrap_tree(tree, opts.output_file, bootstraps)

if __name__ == "__main__":
    main()
