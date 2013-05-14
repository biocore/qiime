#!/usr/bin/env python
# File created on 30 Nov 2009.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.7.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from qiime.parse import parse_mapping_file
    
def merge_mapping_files(mapping_files,no_data_value='no_data'):
    """ Merge list of mapping files into a single mapping file 
    
        mapping_files: open file objects containing mapping data
        no_data_value: value to be used in cases where there is no
         mapping field associated with a sample ID (default: 'no_data')
    """
    mapping_data = {}
    all_headers = []
    result = []
    
    # iterate over mapping files, parsing each
    for mapping_file in mapping_files:
        current_data, current_headers, current_comments = \
           parse_mapping_file(mapping_file,strip_quotes=False)
        all_headers += current_headers
        for entry in current_data:
            sample_id = entry[0]
            current_values = {}
            for header,value in zip(current_headers[1:],entry[1:]):
                current_values[header] = value
            if sample_id in mapping_data:
                # if the sample id has already been seen, confirm that
                # there is no conflicting values across the different 
                # mapping files (e.g., pH=5.0 and pH=6.0)- if there is, 
                # raise a ValueError
                previous_data = mapping_data[sample_id]
                for header,value in current_values.items():
                    if header in previous_data and value != previous_data[header]:
                        raise ValueError,\
                         "Different values provided for %s for sample %s in different mapping files."\
                          % (header,sample_id)
                mapping_data[sample_id].update(current_values)
            else:
                mapping_data[sample_id] = current_values
    all_headers = {}.fromkeys(all_headers)
    
    # remove and place the fields whose order is important
    ordered_beginning = []
    for e in ['SampleID','BarcodeSequence','LinkerPrimerSequence']:
        try:
            del all_headers[e]
            ordered_beginning.append(e)
        except KeyError:
            pass
            
    ordered_end = []
    for e in ['Description']:
        try:
            del all_headers[e]
            ordered_end.append(e)
        except KeyError:
            pass
    ordered_headers = ordered_beginning  + list(all_headers) + ordered_end
    
    
    # generate the mapping file lines containing all fields
    result.append('#' + '\t'.join(ordered_headers))
    for sample_id, data in mapping_data.items():
        result.append('\t'.join([sample_id] + \
          [data.get(h,no_data_value) for h in ordered_headers[1:]]))
    return result

def write_mapping_file(mapping_data,output_fp):
    """ Write list of mapping_data to output_fp """
    f = open(output_fp,'w')
    f.write('\n'.join(mapping_data))
    f.close()
