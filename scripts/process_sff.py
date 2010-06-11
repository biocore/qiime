#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Rob Knight", "Kyle Bittinger", "Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"
 

from qiime.util import parse_command_line_parameters
from optparse import make_option
from qiime.process_sff import prep_sffs_in_dir

script_info={}
script_info['brief_description']="""Convert sff to FASTA and QUAL files"""
script_info['script_description']="""This script converts a directory of sff files into FASTA and QUAL files.

This script requires that 454's off-instrument apps (sffinfo, sfffile) are in your path."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Simple example""","""Convert all the sffs in directory \"sffs/\" to fasta and qual.""","""process_sff.py -i sffs/"""))
script_info['script_usage'].append(("""Flowgram example""","""Convert all the sffs in directory \"sffs/\" to fasta and qual, along with a flowgram file.""","""process_sff.py -i sffs/ -f"""))
script_info['output_description']="""This script results in FASTA and QUAL formatted files."""
script_info['required_options'] = [\
    make_option('-i', '--input_dir', help='Input directory of sff files'),
]
script_info['optional_options']=[\
    make_option('-f', '--make_flowgram', action='store_true', help='this allows for generating a flowgram file. [default: %default]', default=False),
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    prep_sffs_in_dir(opts.input_dir,opts.make_flowgram)

if __name__ == "__main__":
    main()
