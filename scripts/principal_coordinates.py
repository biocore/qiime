#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Rob Knight", "Antonio Gonzalez Pena",\
    "Catherine Lozupone"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Development"
 
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from qiime.principal_coordinates import pcoa, multiple_file_pcoa
from qiime.parse import parse_distmat
import os

options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""Principal Coordinates Analysis (PCoA)"""
script_info['script_description']="""Principal Coordinate Analysis (PCoA) is commonly used to compare groups of samples based on phylogenetic or count-based distance metrics (see section on beta_diversity.py)."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""PCoA (Single File)""","""For this script, the user supplies a distance matrix (i.e. resulting file from beta_diversity.py), along with the output filename (e.g.  beta_div_coords.txt), as follows:""","""%prog -i beta_div.txt -o beta_div_coords.txt"""))
script_info['script_usage'].append(("""PCoA (Multiple Files):""","""The script also functions in batch mode if a folder is supplied as input (e.g. from beta_diversity.py run in batch). No other files should be present in the input folder - only the distance matrix files to be analyzed. This script operates on every distance matrix file in the input directory and creates a corresponding principal coordinates results file in the output directory, e.g.:""","""%prog -i beta_div_weighted_unifrac/ -o beta_div_weighted_pcoa_results/"""))
script_info['output_description']="""The resulting output file consists of the principal coordinate (PC) axes (columns) for each sample (rows). Pairs of PCs can then be graphed to view the relationships between samples. The bottom of the output file contains the eigenvalues and % variation explained for each PC."""
script_info['required_options']=[\

make_option('-i', '--input_path',type='existing_path',\
     help='path to the input distance matrix file(s) (i.e., the output from beta_diversity.py). Is a directory for batch processing and a filename for a single file operation.'),\

make_option('-o', '--output_path',type='new_path',
     help='output path. directory for batch processing, '+\
       'filename for single file operation'),\
]

script_info['optional_options']=[]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if os.path.isdir(opts.input_path):
        multiple_file_pcoa(opts.input_path,opts.output_path)
    elif os.path.isfile(opts.input_path):
        f = open(opts.input_path,'U')
        pcoa_res_string = pcoa(f)
        f.close()

        f = open(opts.output_path, 'w')
        f.write(pcoa_res_string)
        f.close()
    else:
        print("io error, check input file path")
        exit(1)
      



if __name__ == "__main__":
    main()
