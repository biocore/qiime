#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Rob Knight", "Kyle Bittinger"]
__license__ = "GPL"
__version__ = "0.92-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Pre-release"
 

from qiime.util import parse_command_line_parameters
from optparse import make_option
from qiime.sra_spreadsheet_to_map_files import write_map_files

script_info={}
script_info['brief_description']="""Create mapping file from SRA submission spreadsheet"""
script_info['script_description']="""This script reads an SRA submission spreadsheet and generates QIIME mapping files."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Simple example""","""Take an SRA submission spreadsheet input_spreadsheet.txt and write out map files as a series of files input_spreadsheet_[STUDY].txt.map.""","""sra_spreadsheet_to_map_files.py -i input_spreadsheet.txt"""))
script_info['output_description']="""Produces one map file per (STUDY, RUN_PREFIX) combination. Note that the output will include extra stuff not actually needed by QIIME. The intention is just to pull out the info needed for split_libraries.py and downstream analyses. Currently, this does not combine this with the data in the per-sample mapping file."""

script_info['required_options'] = [
    make_option('-i', '--input_file',
        help='the input SRA submission spreadsheet'),
]
script_info['optional_options'] = []
script_info['version'] = __version__



def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    write_map_files(opts.input_file)

if __name__ == "__main__":
    main()
