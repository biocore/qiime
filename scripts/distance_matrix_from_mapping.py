#!/usr/bin/env python
# File created on 23 Jul 2011
# based on beta_diversity.py
from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Development"
 

from qiime.util import parse_command_line_parameters, make_option
from qiime.format import format_distance_matrix
from qiime.distance_matrix_from_mapping import distance_matrix

import os.path

import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')


script_info={}
script_info['brief_description']="""Calculate the pairwise dissimilarity on one column of a mappping file"""
script_info['script_description']="""The input for this script is a mapping file and the name of a column from which a distance matrix will be created. The output of this script is a distance matrix containing a dissimilarity value for each pairwise comparison.

As this is a univariate procedure only one metric is supported: d = c-b."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Pairwise dissimilarity:""","""To calculate the distance matrix (using e.g. euclidean distance) on a column of the mapping file, where the results are output to DOB_euclidean.txt, use the following command:""","""distance_matrix_from_mapping.py -i Fasting_Map.txt -c DOB -m euclidean -o distance_matrix/"""))
script_info['output_description']="""The output of distance_matrix_from_mapping.py is a file containing a distance matrix between rows corresponding to a column in a mapping file."""
script_info['required_options']=[
 make_option('-i', '--input_path',
     help='Mapping filepath.', type='existing_path'),
 make_option('-c', '--column', type='string',\
        help="string containing the name of the column in the mapping file, e.g. 'DOB'"),
]
script_info['optional_options']=[
 make_option('-o', '--output_dir',
     help="Output directory. One will be created if it doesn't exist. [default=%default]",
     type='new_dirpath', default="map_distance_matrix")
]
script_info['option_label']={'input_path':'Mapping filepath',
                             'column':'List of samples for compute',
                             'output_dir': 'Output directory'}
                             
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
  
    if opts.output_dir.endswith('.txt'):
        stderr.write('output must be a directory, files will be named'+\
          ' automatically.  And we refuse to make .txt directories\n')
        exit(1)
    
    try: 
        os.makedirs(opts.output_dir)
    except OSError:
        pass # hopefully dir already exists
    
    dtx_txt = distance_matrix(open(opts.input_path,'U'), opts.column)
    
    outfilepath = os.path.join(opts.output_dir, opts.column + '.txt')
    f = open(outfilepath,'w')
    f.write(dtx_txt)
    f.close()
    

if __name__ == "__main__":
    main()
