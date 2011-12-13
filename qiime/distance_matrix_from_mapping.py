#!/usr/bin/env python
# File created on 27 Sep 2011

from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "1.4.0"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Release"

"""Contains code for performing distance matrix analyses of a column from a mapping file.

command line usage help: python distance_matrix_from_mapping.py -h

This module has the responsibility for creating a distance matrix from sample to sample using
a column from a mapping file:
    - a tab-delimited mapping file
    - the header of the column to relate

The output is a matrix of distances, incl. row/col headers.
    Note that parser expects first field to be blank, i.e. first char of file
    is expected to be a tab.
"""

from qiime.util import FunctionWithParams
from qiime.parse import parse_mapping_file_to_dict
from qiime.format import format_distance_matrix
from numpy import array, reshape
from sys import exit

def distance_matrix(input_path, column):
    """ calculates distance matrix on a single column of a mapping file
    
    inputs:
     input_path (file handler)
     column (str)
    """
    data, comments = parse_mapping_file_to_dict(input_path)
    column_data = []
    column_headers = []
    for i in data:
        if column not in data[i]:
            stderr.write("\n\nNo column: '%s' in the mapping file. Existing columns are: %s\n\n" % (column,data[i].keys()))
            exit(1)
        try:  
            column_data.append(float(data[i][column]))
        except ValueError:
            stderr.write("\n\nall the values in the column '%s' must be numeric but '%s' has '%s'\n\n"\
                % (column,i,data[i][column]))
            exit(1)
            
            
        column_headers.append(i)
    
    data_row = array(column_data)
    data_col = reshape(data_row, (1, len(data_row)))
    dist_mtx = abs(data_row-data_col.T)
    
    return format_distance_matrix(column_headers, dist_mtx)