#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Pre-release"


from qiime.util import parse_command_line_parameters
from optparse import make_option
from qiime.merge_mapping_files import merge_mapping_files, write_mapping_file

script_description = """This script provides a convenient interface to merge mapping which contain data on different samples."""

script_usage = """ Merge two mapping files into a new mapping file (merged_mapping.txt).
  In cases where a mapping field is not provided for some samples, add the value 
  'Data not collected'.
 merge_mapping_files.py -m inseqs1_mapping.txt,inseqs2_mapping.txt  -o merged_mapping.txt -n 'Data not collected'
"""

required_options = [\
 make_option('-m','--mapping_fps',\
         help='the input mapping files in a comma-separated list'),\
 make_option('-o','--output_fp',\
         help='the output mapping file to write'),
]

optional_options = [\
 make_option('-n','--no_data_value',\
         help='value to represent missing data (i.e., when all '+\
         'fields are not defined in all mapping files) [default: %default]',\
         default='no_data'),\
]


def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)
      
    verbose = opts.verbose
    output_fp = opts.output_fp
    mapping_files = [open(fp,'U') for fp in opts.mapping_fps.split(',')]
    no_data_value = opts.no_data_value
    
    mapping_data = merge_mapping_files(mapping_files,\
                                       no_data_value=no_data_value)
    write_mapping_file(mapping_data,output_fp)

if __name__ == "__main__":
    main()