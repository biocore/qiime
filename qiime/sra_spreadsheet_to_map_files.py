#!/usr/bin/env python
#file sra_spreadsheet_to_map_files.py
from qiime.parse import parse_map
from os.path import split, join
from collections import defaultdict
"""This script reads the SRA submission spreadsheet, makes QIIME map files.

Produces one map file per EXPERIMENT -- note that this will include extra
stuff not actually needed by QIIME. Intention is just to pull out the info
needed for split_libaries and downstream analysis. Does not currently combine
this with the data in the per-sample mapping file, but this is planned for
the future.
"""

def strip_quotes(s):
    """strips whitespace and terminal quotes from a string."""
    return s.strip().strip('"')

def collect_study_groups(fields, col_index=0):
    """Returns dict of {study_group:[fields]}"""
    result = defaultdict(list)
    for f in fields:
        result[f[col_index]].append(f)
    return result

def remap_lines(col_names, lines):
    """Remaps fields in lines, substituting 'None' for blank fields."""
    sample_id_index = col_names.index('POOL_MEMBER_NAME')
    barcode_sequence_index = col_names.index('BARCODE')
    result = [['#SampleID', 'BarcodeSequence'] + 
        [i.lstrip('#') for i in col_names] + ['Description']]
    for line in lines:
        result.append([line[sample_id_index], line[barcode_sequence_index]] +
            line + ['None'])
    return result

def write_map_files(infile):
    """Writes map files for each study covered in infile."""
    map_lines = parse_map(infile)
    map_header, map_lines = map_lines[0], map_lines[1:]
    col_names = map(strip_quotes, map_header)
    study_groups = collect_study_groups(map_lines, 
        col_names.index('STUDY_REF'))
    basedir, filename = split(argv[1])
    for name, lines in study_groups.items():
        outfile = open(join(basedir, name)+'.map', 'w')
        for line in remap_lines(col_names, lines):
            outfile.write('\t'.join(line) + '\n')

if __name__ == '__main__':
    from sys import argv
    write_map_files(open(argv[1], 'U'))

