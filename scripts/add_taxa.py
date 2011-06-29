#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight","Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Daniel McDonald"
__email__ = "wasade@gmail.com"
__status__ = "Release"
 

from qiime.util import parse_command_line_parameters
from qiime.util import make_option
from sys import stdout
from qiime.add_taxa import rewrite_otu_table_with_taxonomy

script_info={}
script_info['brief_description']="""Add taxa to OTU table"""
script_info['script_description']="""This script adds taxa to an OTU table."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""Add taxa to otu file from otus.txt from file taxa.txt:""","""add_taxa.py -i otus.txt -t taxa.txt"""))
script_info['output_description']="""The result of this script is written to stdout, or optionally to the specified file."""
script_info['required_options']=[\
    make_option('-i','--otu_file',action='store',\
        type='string',dest='otu_fp',help='Path to read otu file'),
    make_option('-t','--taxonomy_file',action='store',\
        type='string',dest='taxon_fp',help='Path to read taxonomy file')
]

script_info['optional_options']=[\
    make_option('-o','--output_file',action='store',\
        type='string',dest='out_fp',help='Path to write '+\
        'output file [default: stdout]'),
    make_option('-m','--id_map_file',action='store',\
        type='string',dest='id_map_fp',help='Path to read '+\
        'seq id to otu map file [default: %default]')
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if opts.out_fp:
        outfile = open(output_fname, 'w')
    else:
        outfile = stdout
    taxon_lines = open(opts.taxon_fp, 'U')
    if opts.id_map_fp:
        id_map_lines = open(opts.id_map_fp, 'U')
    else:
        id_map_lines = None

    otu_lines = open(opts.otu_fp, 'U')
    rewrite_otu_table_with_taxonomy(taxon_lines, otu_lines, id_map_lines,
        outfile)

if __name__ == "__main__":
    main()
