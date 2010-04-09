#!/usr/bin/env python
#file sra_spreadsheet_to_map_files.py
from qiime.parse import parse_mapping_file
from os.path import split, join
from collections import defaultdict
"""This script reads the SRA submission spreadsheet, makes QIIME map files.

Produces one map file per (STUDY, RUN_PREFIX) combination. Note that the 
output will include extra stuff not actually needed by QIIME. Intention is 
just to pull out the info needed for split_libaries and downstream analysis. 
Does not currently combine this with the data in the per-sample mapping file, 
but this is planned for the future.
"""
__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Rob Knight","Greg Caporaso"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

def strip_quotes(s):
    """strips whitespace and terminal quotes from a string."""
    return s.strip().strip('"')

def collect_study_groups(fields, study_index=0, run_prefix_index=1):
    """Returns dict of {study_group:[fields]}"""
    result = defaultdict(list)
    for f in fields:
        result[(f[study_index], f[run_prefix_index])].append(f)
    return result

def remap_lines(col_names, lines):
    """Remaps fields in lines, substituting 'None' for blank fields."""
    sample_id_index = col_names.index('POOL_MEMBER_NAME')
    barcode_sequence_index = col_names.index('BARCODE')
    primer_sequence_index = col_names.index('PRIMER')
    if 'LINKER' in col_names:
        linker_sequence_index = col_names.index('LINKER')
    else:
        linker_sequence_index = None

    result = [['#SampleID', 'BarcodeSequence', 'LinkerPrimerSequence'] + 
        [i.lstrip('#') for i in col_names] + ['Description']]
    for line in lines:
        curr = [line[sample_id_index], line[barcode_sequence_index]]
        p = line[primer_sequence_index]
        if linker_sequence_index is not None:
            p = line[linker_sequence_index] + p
        curr.append(p)
        curr.extend(line + ['None'])
        result.append(curr)
    return result

def get_study_groups(infile):
    """Return study groups for each study covered in infile."""
    map_lines, map_header, map_comments = parse_mapping_file(infile)
    col_names = map(strip_quotes, map_header)
    study_groups = collect_study_groups(map_lines, 
        col_names.index('STUDY_REF'), col_names.index('RUN_PREFIX'))
    return col_names, study_groups

def write_map_files(in_path):
    """Writes map files for each study covered in infile."""
    infile = open(in_path, 'U')
    col_names, study_groups = get_study_groups(infile)
    basedir, filename = split(in_path)
    for name, lines in study_groups.items():
        outfile = open(join(basedir, '_'.join(name))+'.map', 'w')
        for line in remap_lines(col_names, lines):
            outfile.write('\t'.join(line) + '\n')


