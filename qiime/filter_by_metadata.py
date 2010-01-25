#!/usr/bin/env python
#filter_by_metadata: reads otu table and map, returns only allowed states

from qiime.parse import parse_otus, parse_map
from string import strip
from sys import argv, stdout, stderr
from numpy import array
from StringIO import StringIO
from optparse import OptionParser

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Rob Knight"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "0.9"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Pre-release"

def parse_states(state_string):
    """From string in format 'col1:good1,good2;col2:good1' return dict."""
    result = {}
    state_string = state_string.strip()
    if state_string:
        cols = map(strip, state_string.split(';'))
        for c in cols:
            colname, vals = map(strip, c.split(':'))
            vals = map(strip, vals.split(','))
            result[colname] = set(vals)
    return result

def get_sample_ids(map_data, states):
    """Takes col states in {col:[vals]} format.

    If val starts with !, exclude rather than include.
    
    Combines cols with and, states with or.

    For example, Study:Dog,Hand will return rows where Study is Dog or Hand;
    Study:Dog,Hand;BodySite:Palm,Stool will return rows where Study is Dog
    or Hand _and_ BodySite is Palm or Stool; Study:*,!Dog;BodySite:*,!Stool
    will return all rows except the ones where the Study is Dog or the BodySite
    is Stool.
    """
    name_to_col = dict([(s,map_data[0].index(s)) for s in states])
    good_ids = []
    for row in map_data[1:]:    #remember to exclude header
        include = True
        for s, vals in states.items():
            curr_state = row[name_to_col[s]]
            include = include and (curr_state in vals or '*' in vals) \
                and not '!'+curr_state in vals
        if include:        
            good_ids.append(row[0])
    return good_ids

def find_good_cols(header_line, good_ids):
    """Takes OTU tab-delim header line of sample ids, returns indices.
    
    Note: always includes index 0 (the OTU id); will also include the last
    index if it is taxonomy.
    """
    fields = map(strip, header_line.split('\t'))
    result = [fields.index(i) for i in good_ids if i in fields]
    #include consensus lineage if present
    if fields[-1] == 'Consensus Lineage':
        result.append(-1)
    return [0] + result #always include first field with id

def filter_line(line, good_fields, min_count=None, outfile=stdout):
    """Returns a filtered line of the OTU table, keeping only good samples.

    line: line of the OTU table
    good_fields: list of field indices to include
    min_count: minimum count of OTUs to accept (if None, doesn't try to check
        the counts, so can use to get header info).
    outfile: file handle to write result into (we are doing this rather than
        keeping in memory because the OTU tables are often very large).
    """
    #find slice of the data that will be numeric
    if good_fields[-1] == -1:   #includes taxomony:
        num_slice = slice(1,-1)
    else:
        num_slice = slice(1,None)
    fields = map(strip, line.split('\t'))
    result = [fields[i] for i in good_fields]
    if min_count is not None:
        if array(map(float, result[num_slice]), dtype=float).sum() < min_count:
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

def filter_otus_and_map(map_infile, otu_infile, map_outfile, otu_outfile, 
    valid_states_str, num_seqs_per_otu):
    """Filters OTU and map files according to specified criteria."""
    map_data, header = parse_map(map_infile, return_header=True)
    map_infile.close()
    valid_states = parse_states(valid_states_str)
    sample_ids = get_sample_ids(map_data, valid_states)

    # write out the filtered mapping file
    map_outfile.write('\n'.join(
        map('\t'.join, filter_map(map_data, sample_ids))))
    if not isinstance(map_outfile, StringIO):
        map_outfile.close()

    # write out the filtered OTU file
    for line in otu_infile:
        if line.startswith('#OTU ID'):
            fields = map(strip, line.split('\t'))
            cols = find_good_cols(line, sample_ids)
            filter_line(line, cols, min_count=None, outfile=otu_outfile)
        elif line.startswith('#'):
            otu_outfile.write(line)
        else:
            filter_line(line, cols, min_count=num_seqs_per_otu, 
                outfile=otu_outfile)
    if not isinstance(otu_outfile, StringIO):
        otu_outfile.close()

def make_cmd_parser():
    """Returns command-line options"""
    parser = OptionParser()
    parser.add_option('-i', '--otu', dest='otu_fname',
        help='name of otu file')
    parser.add_option('-m', '--map', dest='map_fname',
        help='name of map file')
    parser.add_option('-o', '--otu_outfile', dest='otu_out_fname', default=None,
        help='name of otu output file, default is otu_filename.filtered.xls')
    parser.add_option('-p', '--map_outfile', dest='map_out_fname', default=None,
        help='name of map output file, default is map_filename.filtered.xls')
    parser.add_option('-s', '--states', dest='valid_states',
        help="string containing valid states, e.g. 'STUDY_NAME:DOG'")
    parser.add_option('-n', '--num_seqs_per_otu', dest='num_seqs_per_otu',
        type=int, default=1, help='minimum counts across samples to keep OTU,' +
        ' default is only to keep OTUs that are present in the samples.')
    options, args = parser.parse_args()
    return options, args

if __name__ == '__main__':
    from sys import exit
    options, args = make_cmd_parser()
    map_file_name, otu_file_name, valid_states_str, num_seqs_per_otu = \
        options.map_fname, options.otu_fname, options.valid_states, \
        options.num_seqs_per_otu

    map_infile = open(map_file_name, 'U')
    otu_infile = open(otu_file_name, 'U')

    map_out_fname = options.map_out_fname
    if map_out_fname is None:
        map_out_fname = map_file_name + '.filtered.xls'

    otu_out_fname = options.otu_out_fname
    if otu_out_fname is None:
        otu_out_fname = map_file_name + '.filtered.xls'


    # write out the filtered mapping file
    map_outfile = open(map_out_fname, 'w')
    otu_outfile = open(otu_out_fname, 'w')

    filter_otus_and_map(map_infile, otu_infile, map_outfile, otu_outfile,
        valid_states_str, num_seqs_per_otu)

