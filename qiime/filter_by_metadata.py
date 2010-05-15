#!/usr/bin/env python
#filter_by_metadata: reads otu table and map, returns only allowed states

from string import strip
from sys import argv, stdout, stderr
from numpy import array
from StringIO import StringIO
from qiime.parse import (parse_mapping_file,
    parse_metadata_state_descriptions)

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Rob Knight", "Antonio Gonzalez Pena"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

"""This filter allows for the removal of sequences and OTUs that either do or don't match specified 
metadata, for instance, isolating samples from a specific set of studies or body sites. This script 
identifies samples matching the specified metadata criteria, and outputs a filtered mapping file 
and OTU table containing only the specified samples."""

def get_sample_ids(map_data, map_header, states):
    """Takes col states in {col:[vals]} format.

    If val starts with !, exclude rather than include.
    
    Combines cols with and, states with or.

    For example, Study:Dog,Hand will return rows where Study is Dog or Hand;
    Study:Dog,Hand;BodySite:Palm,Stool will return rows where Study is Dog
    or Hand _and_ BodySite is Palm or Stool; Study:*,!Dog;BodySite:*,!Stool
    will return all rows except the ones where the Study is Dog or the BodySite
    is Stool.
    """
    
    name_to_col = dict([(s,map_header.index(s)) for s in states])
    good_ids = []
    for row in map_data:    #remember to exclude header
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

def filter_map(map_data, map_header, good_sample_ids):
    """Filters map according to several criteria.

    - keep only sample ids in good_sample_ids
    - drop cols that are different in every sample (except id)
    - drop cols that are the same in every sample
    """
    d = array(map_data) #note, will contain row headers
    first_col = list(d[:,0])    #assume contains sample ids
    good_row_indices = [first_col.index(i) for i in good_sample_ids 
        if i in first_col]
    d = d[good_row_indices]
    cols = d.T
    
    # Determine the good column indices
    # always include the sample id col
    good_col_indices = [0]
    # for all cols between the first and last, 
    # include the col if there are more than
    # 2 states, and less than 1 state per sample ID 
    good_col_indices +=\
     [i+1 for (i, col) in enumerate(cols[1:-1]) 
      if 2 <= len(set(col)) < (len(col))]
    # always include the description col
    good_col_indices += [len(cols)-1]
    h = [map_header[i] for i in good_col_indices]
    d = d[:,good_col_indices]
    
    return h, map(list, d)

def filter_otus_and_map(map_infile, otu_infile, map_outfile, otu_outfile, 
    valid_states_str, num_seqs_per_otu):
    """Filters OTU and map files according to specified criteria."""
    map_data, map_header, map_comments = parse_mapping_file(map_infile)
    map_infile.close()
    valid_states = parse_metadata_state_descriptions(valid_states_str)
    sample_ids = get_sample_ids(map_data, map_header, valid_states)

    # write out the filtered mapping file
    out_headers, out_data = filter_map(map_data, map_header, sample_ids)
    header_line = '#' + '\t'.join(out_headers)
    map_outfile.write('\n'.join([header_line] + map('\t'.join, out_data)))
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

