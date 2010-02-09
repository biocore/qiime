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
from qiime.sra_spreadsheet_to_map_files import write_map_files
script_description = """This script reads the SRA submission spreadsheet, makes QIIME map files.

Produces one map file per (STUDY, RUN_PREFIX) combination. Note that the 
output will include extra stuff not actually needed by QIIME. Intention is 
just to pull out the info needed for split_libaries and downstream analysis. 
Does not currently combine this with the data in the per-sample mapping file, 
but this is planned for the future."""

script_usage = """Take an SRA submission spreadsheet input_spreadsheet.txt and write out map files as a series of files input_spreadsheet_[STUDY].txt.map:

sra_spreadsheet_to_map_files.py -i input_spreadsheet.txt
"""

required_options = [\
    make_option('-i','--input_file',help='the input SRA submission spreadsheet'),\
]

optional_options = [\
 # Example optional option
 #make_option('-o','--output_dir',help='the output directory [default: %default]'),\
]

def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)

    write_map_files(opts.input_file)

if __name__ == "__main__":
    main()
