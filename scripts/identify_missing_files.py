#!/usr/bin/env python
# File created on 01 Jun 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from qiime.util import parse_command_line_parameters, make_option
from os.path import exists

script_info = {}
script_info[
    'brief_description'] = "This script checks for the existence of expected files in parallel runs."
script_info[
    'script_description'] = "This script checks for the existence of expected files in parallel runs, and is useful for checking the status of a parallel run or for finding out what poller.py is waiting on in a possibly failed run."
script_info['script_usage'] = [("Example",
                                "Check for the existence of files listed in expected_out_files.txt from a "
                                "PyNAST alignment run, and print a warning for any that are missing.",
                                "%prog -e ALIGN_BQ7_/expected_out_files.txt")]
script_info['output_description'] = """
This script does not create any output files.
"""
script_info['required_options'] = [
    make_option('-e', '--expected_out_fp',
                type="existing_filepath",
                help='the list of expected output files')
]
script_info['optional_options'] = []
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    filepaths = [l.strip() for l in open(opts.expected_out_fp, 'U')]
    all_exist = True
    for fp in filepaths:
        if not exists(fp):
            print "Filepath doesn't exist: %s" % fp
            all_exist = False
    if all_exist:
        print "All filepaths exist."


if __name__ == "__main__":
    main()
