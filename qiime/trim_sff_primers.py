#!/usr/bin/env python
#file trim_sff_primers.py: resets trim values in sff file based on primers.
from qiime.parse import parse_mapping_file

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"
"""Finds the technical read regions for each library, and resets the left trim.

Replaces the sff files with the trimmed versions.
"""
sfffile_cmd = "%s -t %s -o %s %s"

sffinfo_cmd = '%s %s %s'

def get_technical_lengths(input_map, debug=False):
    """Returns per-sample info on technical lengths.
    
    Note: KEY_SEQ, BARCODE and PRIMER fields are required. LINKER optional.
    """
    if debug:
        print "Making debug output"
    body, header, comments = parse_mapping_file(input_map)
    if debug:
        print "HEADER:", header
    key_index = header.index('KEY_SEQ')
    bc_index = header.index('BARCODE')
    if 'LINKER' in header:
        linker_index = header.index('LINKER')
    else:
        linker_index = None
    primer_index = header.index('PRIMER')
    technical_lengths = {}
    for fields in body:
        curr_tech_len = len(fields[key_index]) + len(fields[bc_index]) + \
            len(fields[primer_index])
        if linker_index is not None:
            curr_tech_len += len(fields[linker_index]) 
        technical_lengths[fields[0]] = curr_tech_len
    if debug:
        print "Technical lengths:"
        print technical_lengths
    return technical_lengths
