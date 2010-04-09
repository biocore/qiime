#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Justin Kuczynski", "Rob Knight", "Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Release"
 
from qiime.util import parse_command_line_parameters, get_options_lookup
from optparse import make_option
from qiime.principal_coordinates import pcoa
from qiime.parse import parse_distmat

options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""Principal Coordinates Analysis (PCoA)"""
script_info['script_description']="""Principal Coordinate Analysis (PCoA) is commonly used to compare groups of samples based on phylogenetic or count-based distance metrics (see section on beta_diversity.py)."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example""","""For this script, the user supplies a distance matrix (i.e., resulting file from beta_diversity.py), along with the output filename (e.g. beta_div_coords.txt), as follows:""","""principal_coordinates.py -i beta_div.txt -o beta_div_coords.txt"""))
script_info['output_description']="""The resulting output file consists of each component (columns) along with the loading for each sample (rows). Pairs of components can then be graphed to view the relationships between samples. The bottom of the output file contains the eigenvalues and % variation explained for each component."""
script_info['required_options']=[\
   make_option('-i', '--input_distance_matrix_fp',\
     help='path to the input distance matrix file (i.e., the output from beta_diversity.py)'),\
   options_lookup['output_fp']\
]
script_info['optional_options']=[]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
      
    infilepath = opts.input_distance_matrix_fp
    f = open(infilepath,'U')
    pcoa_res_string = pcoa(f)
    f.close()

    f = open(opts.output_fp, 'w')
    f.write(pcoa_res_string)
    f.close()


if __name__ == "__main__":
    main()