#!/usr/bin/env python
# File created on 09 Feb 2010
#summarize_otu_by_cat.py
from __future__ import division

__author__ = "Julia Goodrich"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Julia Goodrich"]
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Julia Goodrich"
__email__ = "julia.goodrich@colorado.edu"
__status__ = "Development"
 
from os import getcwd
from qiime.util import parse_command_line_parameters
from optparse import make_option
from qiime.summarize_otu_by_cat import summarize_by_cat


script_info={}
script_info['brief_description']="""Create a summarized OTU table for a specific metadata category"""
script_info['script_description']="""This script generates an otu table where the SampleIDs are replaced by a specific category from the user-generated mapping file. The script uses the OTU file otus.txt (-c) and the user mapping file meta.txt. The user must also specify a metadata category equivalent to one of the column names in the mapping file. If the user wants the counts to be normalized by sample use th normalize flag (-n) the default is False meaning it is only the raw counts. The output is a file called <meta category>_otu_table.txt, it will be put int the current working directory unless specified by the user (-o)."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""Create an otu table for a user specified category. This script uses an OTU table (otu_table.txt) and a user-generated mapping file (mapping_file.txt). The user must also specify a metadata category equivalent to one of the column names in their mapping file (i.e. time). If the user wants the counts to be normalized by sample, they can use the normalize flag (-n), however; the default value for this flag is False, which means it will use the raw counts. The resulting files will be it will be written in the current working directory, unless specified by the user (-o).""","""summarize_otu_by_cat.py -c otu_table.txt -i mapping_file.txt -m time -o qiime_run/ -n"""))
script_info['output_description']="""The output is an otu table called <meta category>_otu_table.txt, """
script_info['required_options']=[\
make_option('-i', '--input_map', dest='map_file',action='store',type='string',\
                help='name of input map file [REQUIRED]'),
make_option('-c', '--otu_file', dest='counts_file',\
            action='store',type='string',
            help='name of otu table file [REQUIRED]'),
make_option('-m', '--meta_category', dest='category',\
                action='store',type='string',
               help='name of category for OTU table [REQUIRED]')
]

script_info['optional_options']=[\
make_option('-o', '--dir-prefix', dest='dir_path',action='store',type='string',\
               help='directory prefix for all analyses [default: cwd]'),
make_option('-n', '--normalize_flag', dest='normalize',
     help='if True will normalize counts [default: %default]',default=False,
                      action = 'store_true')
]

script_info["version"] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

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
