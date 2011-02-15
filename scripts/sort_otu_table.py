#!/usr/bin/env python
# File created on 15 Feb 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.2.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 

from optparse import make_option
from qiime.parse import parse_mapping_file, parse_otu_table
from qiime.format import format_otu_table
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.sort import sort_otu_table

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Script for sorting the sample IDs in an OTU table based on a specified value in a mapping file."
script_info['script_description'] = ""
script_info['script_usage'] = [("",
                                "sort samples by the age field in the mapping file",
                                "sort_otu_table.py -i otu_table.txt -o age_sorted_otu_table.txt -m map.txt -s Age")]
script_info['output_description']= ""
script_info['required_options'] = [
 make_option('-i','--input_otu_table',help='the input otu table'),
 make_option('-o','--output_fp',help='output otu table filepath'),
 make_option('-m','--mapping_fp',help='the mapping file'),
 make_option('-s','--sort_field',help='field to sort by'),
 
]
script_info['optional_options'] = []
script_info['version'] = __version__



def main():
    option_parser, opts, args =\
      parse_command_line_parameters(**script_info)

    mapping_data = parse_mapping_file(open(opts.mapping_fp,'U'))
    otu_table_data = parse_otu_table(open(opts.input_otu_table,'U'))
    
    result = sort_otu_table(otu_table_data,
                           mapping_data,
                           opts.sort_field)

    result_str = format_otu_table(result[0],result[1],result[2],result[3])
    of = open(opts.output_fp,'w')
    of.write(result_str)
    of.close()

if __name__ == "__main__":
    main()