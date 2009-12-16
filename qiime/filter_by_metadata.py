#!/usr/bin/env python
#filter_by_metadata: reads otu table and map, returns only allowed states

from qiime.parse import parse_otus, parse_map
from string import strip
from sys import argv,stdout
from numpy import array
from optparse import OptionParser

__author__ = "Rob Knight"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Rob Knight"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Prototype"

def parse_states(state_string):
    """From string in format 'col1:good1,good2;col2:good1' return dict."""
    result = {}
    cols = map(strip, state_string.split(';'))
    for c in cols:
        colname, vals = map(strip, c.split(':'))
        vals = map(strip, vals.split(','))
        result[colname] = set(vals)
    return result

def get_sample_ids(map_data, states):
    """Takes col states in {col:[good_vals] format.

    TODO: If val starts with !, exclude rather than include.
    Combines cols with and, states with or.
    """
    name_to_col = dict([(s,map_data[0].index(s)) for s in states])
    good_ids = []
    for row in map_data:
        for s, vals in states.items():
            curr_state = row[name_to_col[s]]
            if curr_state in vals or '*' in vals:
                good_ids.append(row[0])
    return good_ids

def find_good_cols(header_line, good_ids):
    """Takes OTU tab-delim header line of sample ids, returns indices"""
    fields = map(strip, header_line.split('\t'))
    result = [fields.index(i) for i in good_ids if i in fields]
    return [0] + result #always include first field with id

def filtered_line(line, good_fields, min_count=0, outfile=stdout):
    fields = map(strip, line.split('\t'))
    result = [fields[i] for i in good_fields]
    if min_count:
        if array(map(int, result[1:]), dtype=int).sum() < min_count:
            return
    outfile.write('\t'.join([fields[i] for i in good_fields])+'\n')

def filter_map(map_data, good_sample_ids):
    """Filters map according to several criteria.

    - keep only sample ids in good_sample_ids
    - drop cols that are different in every sample (except id)
    - drop cols that are the same in every sample
    """
    d = array(map_data) #note, will contain row/col headers
    first_col = list(d[:,0])    #assume contains sample ids
    good_row_indices = [0] + [first_col.index(i) for i in good_sample_ids 
        if i in first_col]
    d = d[good_row_indices]
    cols = d.T
    good_col_indices = [0] + [i+1 for (i, col) in enumerate(cols[1:-1]) 
        if 2 <= len(set(col[1:])) < (len(col) - 1)] + [len(cols)-1]
    d = d[:,good_col_indices]
    return map(list, d)



def make_cmd_parser():
    """Returns command-line options"""
    parser = OptionParser()
    parser.add_option('-o', '--otu', dest='otu_fname',
        help='name of otu file')
    parser.add_option('-m', '--map', dest='map_fname',
        help='name of map file')
    parser.add_option('-s', '--states', dest='valid_states',
        help="string containing valid states, e.g. 'STUDY_NAME:DOG'")
    parser.add_option('-n', '--num_seqs_per_otu', dest='num_seqs_per_otu',
        type=int, default=1, help='minimum counts across samples to keep OTU')
    options, args = parser.parse_args()
    return options, args

if __name__ == '__main__':
    from sys import exit
    options, args = make_cmd_parser()
    map_file_name, otu_file_name, valid_states_str, num_seqs_per_otu = \
        options.map_fname, options.otu_fname, options.valid_states, \
        options.num_seqs_per_otu

    map_file = open(map_file_name, 'U')
    map_data, header = parse_map(map_file, return_header=True)
    map_file.close()

    valid_states = parse_states(valid_states_str)

    sample_ids = get_sample_ids(map_data, valid_states)

    # write out the filtered mapping file
    map_outfile = open(map_file_name+'.filtered.xls', 'w')
    map_outfile.write('\n'.join(
        map('\t'.join, filter_map(map_data, sample_ids))))

    # write out the filtered OTU file
    otu_file = open(otu_file_name, 'U')
    otu_outfile = open(otu_file_name+'.filtered.xls', 'w')
    for line in otu_file:
        if line.startswith('#OTU ID'):
            fields = map(strip, line.split('\t'))
            cols = find_good_cols(line, sample_ids)
            filtered_line(line, cols, min_count=None, outfile=otu_outfile)
        elif line.startswith('#'):
            otu_outfile.write(line)
        else:
            filtered_line(line, cols, min_count=num_seqs_per_otu, 
                outfile=otu_outfile)
        
        

