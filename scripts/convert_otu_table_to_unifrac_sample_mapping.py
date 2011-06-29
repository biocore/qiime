#! /usr/bin/env python

__author__ = "Cathy Lozupone"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Catherine Lozupone"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Cathy Lozupone"
__email__ = "lozupone@colorado.edu"
__status__ = "Release"

from qiime.parse import parse_otu_table
from qiime.util import make_option
from qiime.util import parse_command_line_parameters
from qiime.format import format_unifrac_sample_mapping


script_info={}
script_info['brief_description']="""Convert a QIIME OTU table to a UniFrac sample mapping file"""
script_info['script_description']="""This script allows users who have picked OTUs in QIIME to convert it to a sample mapping (environment) file for use with the Unifrac web interface."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""Convert an OTU table (e.g. otu_table.txt) to a unifrac sample mapping (environment) file: ""","""convert_otu_table_to_unifrac_sample_mapping.py -i otu_table.txt -o sample_mapping.txt"""))
script_info['output_description']="""The result of this script is a sample mapping file for the UniFrac web interface."""
script_info['required_options']=[\
    make_option('-i', '--otu_table_fp', dest='otu_table_fp',\
        help='path to the otu table'),
    make_option('-o', '--output_fp', dest='output_fp', \
        help='path to output file')]

script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    otu_table_fp = opts.otu_table_fp
    output_fp = opts.output_fp
    verbose = opts.verbose
    
    otu_table_lines = open(otu_table_fp, 'U')
    sample_ids, otu_ids, otu_table_array, lineages = \
        parse_otu_table(otu_table_lines, float)
    result = format_unifrac_sample_mapping(sample_ids, otu_ids, otu_table_array)
    of = open(output_fp, 'w')
    of.write('\n'.join(result))
    of.close()

if __name__ == "__main__":
    main()

