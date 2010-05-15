#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "William Walters"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["William Walters"]
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"
__status__ = "Development"
 

from qiime.util import parse_command_line_parameters, get_options_lookup
from optparse import make_option
from qiime.check_id_map import check_mapping_file
from string import letters, digits

#check_id_map.py
options_lookup = get_options_lookup()
script_info={}
script_info['brief_description']="""Checks user's metadata mapping file for required data, valid format"""
script_info['script_description']="""Specifically, we check that:

    - The filename does not contain spaces (warn + rewrite if it does)
    - There are headers for SampleID, LinkerPrimerSequence, and BarcodeSequence if barcodes are used (returns errors if these are absent or misspelled)
    - The BarcodeSequence and LinkerPrimerSequences fields have valid IUPAC DNA characters
    - There are not duplicate header fields (error)
    - There are not duplicate near-unique but not exactly unique values within each column (warning)
    - The headers do not contain invalid characters (alphanumeric and underscore only)
    - The data fields do not contain invalid characters (alphanumeric, underscore, and +-%. characters)
    - There are no duplicates when the primer and barcodes are appended
    
    Errors and warnings are saved to a log file.  Errors are generally caused 
    by problems with the headers, and should be resolved before attempting to 
    correct any warnings.  Warnings can arise from invalid characters, 
    near-duplicate metadata, duplicate sample descriptions/barcodes, or missing
    data fields. Warnings will contain a reference to the cell (row,column) 
    that the warning arose from.
    
    In addition to the log file, a "corrected_mapping" file will be created.
    Invalid characters will be replaced by underscores in this corrected mapping
    file if there were any such characters in the input metadata mapping file.
    If there were no invalid characters to replace, the corrected mapping file 
    will contain comments saying as much.
    
    check_id_map.py should not raise exceptions itself under normal 
    circumstances, except for situations such as having a misformatted input 
    metadata mapping file.
"""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""Check the test_mapping.txt mapping file for problems, supplying the required mapping file and output directory (in this case mapping_info)""","""check_id_map.py -m test_mapping.txt -o mapping_info/"""))
script_info['output_description']="""A log file and corrected_mapping.txt file will be written to the mapping_info directory."""
script_info['required_options']= [\
    make_option('-m', '--map', dest='map_fname',
        help='Metadata mapping file filepath'),
    make_option('-o', '--output_dir',
        help='Required output directory for log file and corrected mapping '+\
        'file (by default, invalid characters will be '+\
        'converted to underscores)')
]
script_info['optional_options']= [\
    make_option('-c', '--char_replace', dest='char_replace',
        help='Changes the default character used to replace invalid '+\
        'characters found in the mapping file.  Must be a valid character ('+\
        'alphanumeric or underscore).  NOT IMPLEMENTED CURRENTLY '+\
        '[default: %default]', default="_"),
    make_option('-b', '--not_barcoded',
        action='store_true', default=False,
        help='Use -b if barcodes are not present. [default: %default]'),
    make_option('-B', '--variable_len_barcodes',
        action='store_true', default=False,
        help='Use -B if variable length barcodes are present to suppress '+\
        'warnings about barcodes of unequal length. [default: %default]'),
    make_option('-p', '--disable_primer_check',
        action='store_true', default=False,
        help='Use -p to disable checks for primers. [default: %default]')
]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
      
    infile_name = opts.map_fname
    has_barcodes = not opts.not_barcoded
    output_dir = opts.output_dir
    char_replace = opts.char_replace
    var_len_barcodes = opts.variable_len_barcodes
    verbose = opts.verbose
    disable_primer_check = opts.disable_primer_check
    
    valid_replacement_chars=digits+letters+"_"
    if char_replace not in valid_replacement_chars:
        option_parser.error('-c option requires alphanumeric or '+\
        'underscore character')
    if len(char_replace) != 1:
        option_parser.error('-c parameter must be a single character.')
    
    check_mapping_file(infile_name, output_dir, has_barcodes, char_replace,\
     verbose, var_len_barcodes, disable_primer_check)


if __name__ == "__main__":
    main()
