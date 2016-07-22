#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Kyle Bittinger"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"

from qiime.util import parse_command_line_parameters
from qiime.trim_sff_primers import (
    get_technical_lengths, set_sff_trimpoints, set_sff_trimpoints_with_sfftools,
)
from qiime.util import make_option

script_info = {}
script_info['brief_description'] = """Trim sff primers"""
script_info[
    'script_description'] = """Finds the technical read regions for each library, and resets the left trim."""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Simple example""",
     """Trim a directory of per-sff files in sff_dir (-l sff_dir/) using an input map (-m input_map.txt). This script uses the sff utility binaries which must be in your path.""",
     """trim_sff_primers.py -l sff_dir/ -m input_map.txt"""))
script_info[
    'output_description'] = """This script replaces the original sff files with the trimmed versions."""

script_info['required_options'] = [
    make_option("-l", "--libdir", dest='libdir',
                type='existing_path',
                help="The directory containing per-library sff files"),
    make_option("-m", "--input_map", dest='input_map',
                type='existing_filepath',
                help="Path to the input mapping file describing the libraries"),
]

script_info['optional_options'] = [
    make_option("-p", "--sfffile_path", default='sfffile', type='string',
                help="Path to sfffile binary [default: %default]"),
    make_option("-q", "--sffinfo_path", default='sffinfo', type='string',
                help="Path to sffinfo binary [default: %default]"),
    make_option('--use_sfftools', action='store_true', default=False,
                help=('Use external sffinfo and sfffile programs instead of '
                      'equivalent Python implementation.')),
    make_option('--debug', default=False, action='store_true',
                help="Print command-line output for debugging [default: %default]"),
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    technical_lengths = get_technical_lengths(
        open(opts.input_map, 'U'), opts.debug)

    if opts.use_sfftools:
        set_sff_trimpoints_with_sfftools(
            opts.libdir, technical_lengths, opts.sffinfo_path,
            opts.sfffile_path, opts.debug)
    else:
        set_sff_trimpoints(opts.libdir, technical_lengths)

if __name__ == "__main__":
    main()
