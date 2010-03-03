#!/usr/bin/env python
# File created on 30 Nov 2009.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "0.91"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

from optparse import OptionParser
from qiime.parse import parse_map
    
def merge_mapping_files(mapping_files,no_data_value='no_data'):
    """ Merge list of mapping files into a single mapping file 
    
        mapping_files: open file objects containing mapping data
        no_data_value: value to be used in cases where there is no
         mapping field associated with a sample ID (default: 'no_data')
    """
    mapping_lines = []
    all_headers = {}
    result = []
    
    # iterate over mapping files, parsing each
    for mapping_file in mapping_files:
        # the parse map doesn't function as it says it should, but 
        # when that gets fixed, the following two lines should be 
        # replaced with the comments lines just below them
        data = parse_map(mapping_file,return_header=True,strip_quotes=False)
        current_headers = data[0][0]
        #data, current_headers = \
        # parse_map(mapping_file,return_header=True,strip_quotes=False)
        all_headers.update(dict.fromkeys(current_headers))
        for d in data[0][1:]:
            current_values = {}
            for i,v in enumerate(d):
                current_values[current_headers[i]] = v
            mapping_lines.append(current_values)
    
    # remove and place the fields whose order is important
    del all_headers['#SampleID']
    del all_headers['BarcodeSequence']
    del all_headers['Description']
    all_headers = ['#SampleID','BarcodeSequence'] \
     + list(all_headers) + ['Description']
    
    # generate the mapping file lines containing all fields
    result.append('\t'.join(all_headers))
    for mapping_line in mapping_lines:
        result.append('\t'.join(\
         [mapping_line.get(h,no_data_value) for h in all_headers]))
    
    return result

def write_mapping_file(mapping_data,output_fp):
    """ Write list of mapping_data to output_fp """
    f = open(output_fp,'w')
    f.write('\n'.join(mapping_data))
    f.close()
