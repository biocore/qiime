#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "William Walters"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["William Walters"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"
__status__ = "Pre-release"
 

from qiime.util import parse_command_line_parameters
from optparse import make_option
from qiime.check_id_map import check_mapping_file
from string import letters, digits

script_description = """ Parse mapping file, checking for a number of undesirable characteristics.

Specifically, we check that:
    - the filename does not contain spaces (warn + rewrite if it does)
    - there is an overall run description at the start (warn + add if missing)
    - there is a SampleID field and, if barcoded, a BarcodeSequence field
      (warn + correct reasonable variants of both)
    - there is a LinkerPrimerSequence field
    - The BarcodeSequence and LinkerPrimerSequences have valid IUPAC DNA chars
    - there are not duplicate header fields (error)
    - there are not duplicate near-unique but not exactly unique values 
      within each column (warning)

Overall strategy:
    - maintain list of errors and warnings (initially empty).
    - for each check in the following classes:
      - check_name
      - check_run_description
      - check_col_headers
      - check_cols
    - run the check f(data) -> msg
    - collect the messages, assigning to error or warning
    - return errors and warnings

Returns both errors and warnings as lists of formatted strings: should not
raise exceptions itself under normal circumstances (e.g. if file is
mis-formatted).

SampleID column is required - must contain unique values. 
BarcodeSequence field required (if is barcoded) and must contain unique values

It is somewhat inefficient to read in the whole table, but on the other hand
if reading in the mapping file is the bottleneck the downstream analysis is
likely to prove somewhat challenging as well..."""

script_usage = """
Check the test_mapping.txt mapping file for problems, supplying the required
mapping file and output directory (in this case mapping_info).\n\n 
check_id_map.py -m test_mapping.txt -o mapping_info/ """

required_options = [\
 make_option('-m', '--map', dest='map_fname',
        help='Mapping file filepath'),\
 make_option('-o', '--output_dir',
        help='Required output directory for mapping file with corrected '+\
        'characters and log file (by default, invalid characters will be '+\
        'converted to underscores)')
]

optional_options = [\
 make_option('-c', '--char_replace', dest='char_replace',
        help='Changes the default character used to replace invalid '+\
        'characters found in the mapping file.  Must be a valid character ('+\
        'alphanumeric or underscore).  NOT IMPLEMENTED CURRENTLY '+\
        '[default: %default]', default="_"),\
 make_option('-b', '--not_barcoded',
        action='store_true', default=False,
        help='Use -b if barcodes are not present. [default: %default]')
]




def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)
      
    infile_name = opts.map_fname
    has_barcodes = not opts.not_barcoded
    output_dir = opts.output_dir
    char_replace = opts.char_replace
    verbose = opts.verbose
    
    valid_replacement_chars=digits+letters+"_"
    if char_replace not in valid_replacement_chars:
        option_parser.error('-c option requires alphanumeric or underscore character')
    if len(char_replace) != 1:
        option_parser.error('-c parameter must be a single character.')
    
    check_mapping_file(infile_name, output_dir, has_barcodes, char_replace,\
     verbose)


if __name__ == "__main__":
    main()
