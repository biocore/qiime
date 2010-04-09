#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Rob Knight", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Release"
 

from qiime.util import parse_command_line_parameters
from optparse import make_option
from os import walk, system
from os.path import splitext, join

#make_per_library_sff.py
script_info={}
script_info['brief_description']="""Make per-library sff files from id lists"""
script_info['script_description']="""This script generates per-library sff files using the id lists."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""Make per-library sff files using input.sff and a directory of libs where each file in the directory contains the id lists for each library:""","""make_per_library_sff.py -i input.sff -l libs"""))
script_info['output_description']="""The result of this script generates sff files for each library."""

script_info['required_options']=[\
    make_option("-i","--input_sff",dest='in_sff',\
        help="The path to an input sff file (or files: separate w/ comma, no spaces)"),
    make_option("-l","--libdir",dest='libdir',\
        help=""" The directory containing per-library id files""")
]

script_info['optional_options']=[\
    make_option("-p","--sfffile_path",dest='sfffile_path',\
        help=""" Path to sfffile binary [default: %default]""", 
        default='sfffile'),
    make_option('--debug', dest='debug', default=False,
        action='store_true',
        help="Print command-line for debugging [default: %default]")
]
script_info['version'] = __version__

cmd = "%s -i %s -o %s %s"

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if opts.debug:
        print "Making debug output"
    input_sff_names = opts.in_sff.replace(',',' ')
    for dirpath, dirnames, fnames in walk(opts.libdir):
        for fname in fnames:
            if fname.startswith('.'):
                continue
            basename, ext = splitext(fname)
            cmd_str = cmd % (opts.sfffile_path, join(dirpath,fname),
                join(dirpath, basename+'.sff'), input_sff_names)
            if opts.debug:
                print cmd_str
            system(cmd_str)


if __name__ == "__main__":
    main()
