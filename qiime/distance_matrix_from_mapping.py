#!/usr/bin/env python
# based on beta_diversity.py

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "1.3.0-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Development"

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
from StringIO import StringIO
from sys import exit, stderr
import os.path

import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')
import cogent.maths.distance_transform as distance_transform #avoid hard-coding metrics

from qiime.util import FunctionWithParams
from qiime.parse import parse_mapping_file_to_dict
from qiime.format import format_distance_matrix
from numpy import array, diagflat

def get_nonphylogenetic_metric(name):
    """Gets metric by name from distance_transform.
    
    Metrics should be f(matrix) -> distances.
    """
    # looks for name, inserting possible dist_ to find functions
    # in distance_transform.py named e.g.:
    # binary_dist_chisq / dist_bray_curtis
    try:
        return getattr(distance_transform, 'dist_' + name.lower())
    except AttributeError:
        try:
            return getattr(distance_transform, 
              name.replace('binary', 'binary_dist').lower())
        except AttributeError:
            return getattr(distance_transform, 
              name.lower())

def list_known_nonphylogenetic_metrics():
    """Lists known metrics by name from distance_transform.

    Assumes that functions starting with dist_ or binary_dist are metrics.
    """
    result = []
    for name in dir(distance_transform):
        if name.startswith('dist_'):
            result.append(name[5:])
        elif name.startswith('binary_dist_'):
            result.append('binary_' + name[12:])
    result.sort()
    return result

def list_known_metrics():
    return list_known_nonphylogenetic_metrics()

def distance_matrix(input_path, output_dir, metrics, column):
    """ does distance matrix calc on a single column of a mapping file
    
    uses name in metrics to name outputdistance matrix files
    inputs:
     input_path (str)
     output_dir (str)
     metrics (str, comma delimited if more than 1 metric)
     column (str)
    """
    f = open(input_path,'U')
    data, comments = parse_mapping_file_to_dict(f)
    f.close()
    
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
    column_data = diagflat(array(column_data))
    
    metrics_list = metrics.split(',')
    for metric in metrics_list:
        outfilepath = os.path.join(output_dir, column + '_' + metric + '.txt')
        try:
            metric_f = get_nonphylogenetic_metric(metric)
            is_phylogenetic = False
        except AttributeError:
            stderr.write("you can not use phylogenetic metrics: %s\n" % (metric,))
            exit(1)
        dissims = metric_f(column_data.T)

        f = open(outfilepath,'w')
        f.write(format_distance_matrix(column_headers, dissims))
        f.close()