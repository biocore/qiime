#!/usr/bin/env python
# File created on 30 Nov 2009.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from qiime.parse import parse_mapping_file
from collections import defaultdict

def merge_mapping_files(mapping_files,no_data_value='no_data'):
    """ Merge list of mapping files into a single mapping file 
    
        mapping_files: open file objects containing mapping data
        no_data_value: value to be used in cases where there is no
         mapping field associated with a sample ID (default: 'no_data')
    """
    mapping_data = defaultdict(dict)
    all_headers = set([])
    
    # iterate over mapping files, parsing each
    for mapping_file in mapping_files:
        current_data, current_headers, current_comments = \
           parse_mapping_file(mapping_file,strip_quotes=False)
        all_headers.update(set(current_headers))

        for entry in current_data:
            current_values = {k:v for k,v in zip(current_headers, entry)}
            sample_id = current_values['SampleID']

            if sample_id in mapping_data:
                # if the sample id has already been seen, confirm that
                # there is no conflicting values across the different 
                # mapping files (e.g., pH=5.0 and pH=6.0)- if there is, 
                # raise a ValueError
                previous_data = mapping_data[sample_id]
                
                for key in current_values:
                    if key not in previous_data:
                        continue

                    if current_values[key] != previous_data[key]:
                        raise ValueError("Different values provided for %s for"
                                      "sample %s in different mapping files."\
                                      % (key, sample_id))

            mapping_data[sample_id].update(current_values)
    
    # remove and place the fields whose order is important
    ordered_beginning = []
    for e in ['SampleID','BarcodeSequence','LinkerPrimerSequence']:
        if e in all_headers:
            all_headers.remove(e)
            ordered_beginning.append(e)
            
    ordered_end = []
    for e in ['Description']:
        if e in all_headers:
            all_headers.remove(e)
            ordered_end.append(e)
    
    ordered_headers = ordered_beginning + list(all_headers) + ordered_end
    
    # generate the mapping file lines containing all fields
    result = ['#' + '\t'.join(ordered_headers)]
    for sample_id, data in mapping_data.items():
        values = [data.get(k, no_data_value) for k in ordered_headers]
        result.append('\t'.join(values))
    
    return result

def write_mapping_file(mapping_data,output_fp):
    """ Write list of mapping_data to output_fp """
    f = open(output_fp,'w')
    f.write('\n'.join(mapping_data))
    f.write('\n')
    f.close()
