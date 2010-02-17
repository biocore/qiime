#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Julia Goodrich"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Julia Goodrich"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Julia Goodrich"
__email__ = "julia.goodrich@colorado.edu"
__status__ = "Pre-release"
 
from os import getcwd
from qiime.util import parse_command_line_parameters
from optparse import make_option
from qiime.summarize_otu_by_cat import summarize_by_cat

script_description = """This script generates the otu table for a specific \
category"""

script_usage = """ Create otu table for a user specified category. The script\
uses the OTU file otus.txt (-c) and the user mapping file meta.txt. The user\
must also specify a metadata category equivalent to one of the column names \
in the mapping file. If the user wants the counts to be normalized by sample \
use th normalize flag (-n) the default is False meaning it is only the raw \
counts. The output is a file called <meta category>_otu_table.txt, it will be \
put int the current working directory unless specified by the user (-o).

python ~/code/Qiime/trunk/qiime/make_otu_network.py -c otus.txt -i \
meta.txt -m time -o /Users/bob/qiime_run/ -n"""

required_options = [\
make_option('-i', '--input_map', dest='map_file',action='store',type='string',\
                help='name of input map file [REQUIRED]'),
make_option('-c', '--otu_file', dest='counts_file',\
			action='store',type='string',
            help='name of otu table file [REQUIRED]'),
make_option('-m', '--meta_category', dest='category',\
				action='store',type='string',
               help='name of category for OTU table [REQUIRED]')
]

optional_options = [\
make_option('-o', '--dir-prefix', dest='dir_path',action='store',type='string',\
               help='directory prefix for all analyses [default: cwd]'),
make_option('-n', '--normalize_flag', dest='normalize',
     help='if True will normalize counts [default: %default]',default=False,
                      action = 'store_true')
]



def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)

    if not opts.counts_file:
        parser.error("An otu table file must be specified")

    if not opts.map_file:
        parser.error("A Map file must be specified")

    dir_path = opts.dir_path
    category = opts.category
    norm = opts.normalize

    if dir_path == "./" or dir_path is None:
        dir_path = getcwd()

    map_lines = open(opts.map_file,'U').readlines()
    otu_sample_lines = open(opts.counts_file,'U').readlines()

    summarize_by_cat(map_lines,otu_sample_lines,category,dir_path,norm)

if __name__ == "__main__":
    main()
