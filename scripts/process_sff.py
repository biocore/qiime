#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Kyle Bittinger", "Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Release"
 

from qiime.util import parse_command_line_parameters,get_options_lookup
from qiime.util import make_option
from os import mkdir
from os.path import splitext,split,join,isfile, isdir
from qiime.process_sff import prep_sffs_in_dir

options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""Convert sff to FASTA and QUAL files"""
script_info['script_description']="""This script converts a directory of sff files into FASTA, QUAL and flowgram files.
"""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Simple example""","""Convert all the sffs in directory \"sffs/\" to fasta and qual.""","""process_sff.py -i sffs/"""))
script_info['script_usage'].append(("""Flowgram example""","""Convert all the sffs in directory \"sffs/\" to fasta and qual, along with a flowgram file.""","""process_sff.py -i sffs/ -f"""))
script_info['script_usage'].append(("""Output example""","""Convert all the sffs in directory \"sffs/\" to fasta and qual, along with a flowgram file and write them to another directory.""","""process_sff.py -i sffs/ -f -o output_dir"""))
script_info['output_description']="""This script results in FASTA and QUAL formatted files."""
script_info['required_options'] = [\
    make_option('-i', '--input_dir', help='Input directory of sff files'),
]

script_info['optional_options'] = [
    make_option('-f', '--make_flowgram', action='store_true', default=False,
        help='generate a flowgram file. [default: %default]'),
    make_option('-t', '--convert_to_FLX', action='store_true', default=False,
        help='convert Titanium reads to FLX length. [default: %default]'),
    make_option('--use_sfftools', action='store_true', default=False,
        help=('use the external programs sfffile and sffinfo for processing, '
              'instead of the equivalent python implementation')),
    make_option('-o', '--output_dir',default=None,
     help='Input directory of sff files [default: same as input dir]'),
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
    else:
        opts.output_dir = opts.input_dir
            
    prep_sffs_in_dir(
        opts.input_dir,
        opts.output_dir,
        make_flowgram=opts.make_flowgram,
        convert_to_flx=opts.convert_to_FLX,
        use_sfftools=opts.use_sfftools,
        )

if __name__ == "__main__":
    main()
