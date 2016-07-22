#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Kyle Bittinger", "Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"


from qiime.util import (parse_command_line_parameters, get_options_lookup,
                        make_option, create_dir)
from os import mkdir
from os.path import isdir, isfile, split
from os.path import splitext, split, join, isfile, isdir
from qiime.process_sff import prep_sffs_in_dir

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = """Convert sff to FASTA and QUAL files"""
script_info['script_description'] = """This script converts a directory of sff files into FASTA, QUAL and flowgram files.
"""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Simple example""",
     """Convert all the sffs in directory \"sffs/\" to fasta and qual.""",
     """process_sff.py -i sffs/"""))
script_info['script_usage'].append(
    ("""""",
     """Convert a single sff to fasta and qual.""",
     """process_sff.py -i sffs/test.sff"""))
script_info['script_usage'].append(
    ("""Flowgram example""",
     """Convert all the sffs in directory \"sffs/\" to fasta and qual, along with a flowgram file.""",
     """process_sff.py -i sffs/ -f"""))
script_info['script_usage'].append(
    ("""""",
     """Convert a single sff to fasta and qual, along with a flowgram file.""",
     """process_sff.py -i sffs/test.sff -f"""))
script_info['script_usage'].append(
    ("""Output example""",
     """Convert all the sffs in directory \"sffs/\" to fasta and qual, along with a flowgram file and write them to another directory.""",
     """process_sff.py -i sffs/ -f -o output_dir"""))
script_info[
    'output_description'] = """This script results in FASTA and QUAL formatted files."""
script_info['required_options'] = [
    make_option(
        '-i',
        '--input_dir',
        type='existing_path',
        help='Input directory of sff files or a single sff filepath'),
]

script_info['optional_options'] = [
    make_option('--no_trim', action='store_true', default=False,
                help='do not trim sequence/qual (requires --use_sfftools option) [default: %default]'),
    make_option('-f', '--make_flowgram', action='store_true', default=False,
                help='generate a flowgram file. [default: %default]'),
    make_option('-t', '--convert_to_FLX', action='store_true', default=False,
                help='convert Titanium reads to FLX length. [default: %default]'),
    make_option('--use_sfftools', action='store_true', default=False,
                help=('use the external programs sfffile and sffinfo for processing, '
                      'instead of the equivalent python implementation')),
    make_option('-o', '--output_dir', default=None, type='new_dirpath',
                help='Input directory of sff files [default: same as input dir]'),
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    output_dir = opts.output_dir

    if output_dir:
        create_dir(output_dir)
    else:
        if isfile(opts.input_dir):
            # if output_dir is empty after the split, then a relative path was
            # passed, and the input file is in the current directory
            output_dir = split(opts.input_dir)[0] or '.'

        else:  # opts.input_dir is a directory
            output_dir = opts.input_dir

    if opts.no_trim and not opts.use_sfftools:
        raise ValueError(
            "When using the --no_trim option you must have the sfftools installed and must also pass the --use_sfftools option")

    prep_sffs_in_dir(
        opts.input_dir,
        output_dir,
        make_flowgram=opts.make_flowgram,
        convert_to_flx=opts.convert_to_FLX,
        use_sfftools=opts.use_sfftools,
        no_trim=opts.no_trim)

if __name__ == "__main__":
    main()
