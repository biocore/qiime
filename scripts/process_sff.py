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
 

from qiime.util import parse_command_line_parameters,get_options_lookup
from optparse import make_option
from os import mkdir
from os.path import splitext,split,join,isfile, isdir
from qiime.process_sff import prep_sffs_in_dir

options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""Convert sff to FASTA and QUAL files"""
script_info['script_description']="""This script converts a directory of sff files into FASTA, QUAL and flowgram files.

This script requires that 454's off-instrument apps (sffinfo, sfffile) are in your path."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Simple example""","""Convert all the sffs in directory \"sffs/\" to fasta and qual.""","""process_sff.py -i sffs/"""))
script_info['script_usage'].append(("""Flowgram example""","""Convert all the sffs in directory \"sffs/\" to fasta and qual, along with a flowgram file.""","""process_sff.py -i sffs/ -f"""))
script_info['script_usage'].append(("""Output example""","""Convert all the sffs in directory \"sffs/\" to fasta and qual, along with a flowgram file and write them to another directory.""","""process_sff.py -i sffs/ -f -o output_dir"""))
script_info['output_description']="""This script results in FASTA and QUAL formatted files."""
script_info['required_options'] = [\
    make_option('-i', '--input_dir', help='Input directory of sff files'),
]
script_info['optional_options']=[\
    make_option('-f', '--make_flowgram', action='store_true', help='this allows for generating a flowgram file. [default: %default]', default=False),
    make_option('-t', '--convert_to_FLX', action='store_true', help='this converts Titanium length read to FLX length. [default: %default]', default=False),
    options_lookup['output_dir'],
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    if opts.output_dir:
        #try to make the output directory
        try:
            mkdir(opts.output_dir)
        except OSError:
            pass
            
        #create the output pathname depending on whether a file is supplied
        #or a directory
        if isdir(opts.input_dir):
            output_pathname=opts.output_dir
        elif isfile(opts.input_dir):
            output_pathname=join(opts.output_dir, \
                                    split(splitext(opts.input_dir)[0])[-1])
    else:
        #create the output pathname depending on whether a file is supplied
        #or a directory
        if isdir(opts.input_dir):
            output_pathname=opts.input_dir
        elif isfile(opts.input_dir):
            output_pathname=splitext(opts.input_dir)[0]
            
    prep_sffs_in_dir(opts.input_dir,opts.make_flowgram,output_pathname,\
                        opts.convert_to_FLX)

if __name__ == "__main__":
    main()
