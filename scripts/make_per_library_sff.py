#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Pre-release"
 

from qiime.util import parse_command_line_parameters
from optparse import make_option
from os import walk, system
from os.path import splitext, join

script_description = """Makes per-library sff files from id lists."""

script_usage = """python make_per_library_sff.py {-i input_sff -l lib_dir} [options]

Example: Make per-library sff files using input.sff and a directory of libs where each file in the directory contains the id lists for each library.

python make_per_library_sff.py -i input.sff -l libs
"""

required_options = [\
    make_option("-i","--input_sff",dest='in_sff',\
        help="The path to an input sff file (or files: separate w/ comma, no spaces)"),
    make_option("-l","--libdir",dest='libdir',\
        help=""" The directory containing per-library id files""")
]

optional_options = [\
    make_option("-p","--sfffile_path",dest='sfffile_path',\
        help=""" Path to sfffile binary [default: %default]""", 
        default='sfffile'),
    make_option('--debug', dest='debug', default=False,
        action='store_true',
        help="Print command-line for debugging [default: %default]")
]

cmd = "%s -i %s -o %s %s"

def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)

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
