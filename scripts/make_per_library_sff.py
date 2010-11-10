#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Rob Knight", "Daniel McDonald"]
__license__ = "GPL"
__status__ = "1.2.0-dev"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Development"
 

from optparse import make_option

from qiime.make_per_library_sff import make_per_library_sffs
from qiime.util import parse_command_line_parameters


#make_per_library_sff.py
script_info={}
script_info['brief_description']="""Make per-library sff files from id lists"""
script_info['script_description']="""This script generates per-library sff files using the id lists."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""Make per-library sff files using input.sff and a directory of libs where each file in the directory contains the id lists for each library:""","""make_per_library_sff.py -i input.sff -l libs"""))
script_info['output_description']="""The result of this script generates sff files for each library."""

script_info['required_options'] = [
    make_option("-i", "--input_sff",
        help="Input sff file (separate multiple files w/ comma)"),
    make_option("-l", "--libdir",
        help="Directory containing per-library id files"),
]

script_info['optional_options'] = [
    make_option("-p", "--sfffile_path",
        help="Path to sfffile binary [default: use sfffile in $PATH]"),
    make_option('--use_sfftools', action='store_true', default=False,
        help=('Use external sfffile program instead of equivalent Python '
              'routines.')),
    make_option('--debug', action='store_true', default=False,
        help="Print debugging output to stdout [default: %default]"),
]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    sff_fps = opts.input_sff.split(',')
    make_per_library_sffs(
        sff_fps,
        opts.libdir,
        opts.use_sfftools,
        opts.sfffile_path,
        opts.debug,
        ) 

if __name__ == "__main__":
    main()
