#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Pre-release"
 

from qiime.util import parse_command_line_parameters
from qiime.tree_compare import load_tree_files, bootstrap_support, write_bootstrap_support_files
from optparse import make_option
import os.path
from os.path import exists

script_description = """Compares jackknife/bootstrapped trees with a master, \
and outputs support for nodes.  The primary function in this module is \
bootstrap_support() """

script_usage = """usage: %prog [options] {-o OUTPUT_DIR -m MASTER_TREE -s SUPPORT_DIR}

Example usage:
python %prog -m sample_cluster.tre -s rare_unifrac_upgma/ -o unifrac_jackknife/
this makes the folder unifrac_jackknife.  In that is the master tree,
with internal nodes named uniquely, a separate bootstrap/jackknife support file,
and a jackknife_named_nodes.tre tree, for use with e.g.: figtree

output jackknife support values are in the range [0,1]

master tree must have the same tips as support trees.  if your support trees
omit some tips (e.g.: samples with few sequences),
make a new master tree with those tips omitted"""

required_options = [\
make_option('-m', '--master_tree', help='master tree filepath'),\
make_option('-s', '--support_dir', help='path to dir containing support trees'),\
make_option('-o', '--output_dir', help='output directory, writes three files here '+\
"makes dir if it doesn't exist")
]

optional_options = [\
 # Example optional option
 #make_option('-o','--output_dir',help='the output directory [default: %default]'),\
]

def main():
    option_parser, options, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)

    if not os.path.exists(options.output_dir):
        os.makedirs(options.output_dir)

    master_tree, support_trees = load_tree_files(options.master_tree,
        options.support_dir)
    # get support of each node in master
    new_master, bootstraps = bootstrap_support(master_tree, support_trees)

    write_bootstrap_support_files(new_master, bootstraps, options.output_dir,
    len(support_trees))

if __name__ == "__main__":
    main()