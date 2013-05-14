#!/usr/bin/env python
# File created on 25 May 2010
from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Antonio Gonzalez Pena", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.7.0"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Development"

from biom.table import table_factory
from qiime.util import make_option
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.parse import parse_trflp
from qiime.format import format_biom_table
from os.path import isfile
from sys import stderr

options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""Convert TRFLP text file to an OTU table"""
script_info['script_description']="""The input for this script is a TRLFP text file. The output of this script is an OTU table text file that can be use with QIIME for further analysis """
script_info['script_usage'] = [
  ("""Usage:""","""You need to pass a TRFLP text, the script will remove not wanted chars sample and otus names, and will add zeros as need it""","""%prog -i trflp_in.txt -o otu_table.biom""")
]
script_info['output_description']= ""
script_info['required_options']=[
  make_option('-i', '--input_path',
     type='existing_filepath',
     help='input path: TRFLP text file'),
  make_option('-o', '--output_path',
     type='new_filepath',
     help="output file: OTU table"),
]
script_info['optional_options']=[]
script_info['version'] = __version__


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
 
    if not isfile(opts.input_path):
       raise IOError, \
        "Input path (%s) not valid.  Does it exist?" % opts.input_path
    
    samples, otus, data = parse_trflp(open(opts.input_path,'U'))
    
    output_f = open(opts.output_path, 'w')
    t = table_factory(data,samples,otus)
    output_f.write(format_biom_table(t))
    output_f.close()
 

if __name__ == "__main__":
    main()
