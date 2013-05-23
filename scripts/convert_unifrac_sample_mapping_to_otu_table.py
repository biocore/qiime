#! /usr/bin/env python

__author__ = "Cathy Lozupone"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Catherine Lozupone", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Cathy Lozupone"
__email__ = "lozupone@colorado.edu"
__status__ = "Development"

from qiime.parse import sample_mapping_to_biom_table
from qiime.util import make_option
from qiime.util import parse_command_line_parameters
from qiime.format import format_biom_table


script_info={}
script_info['brief_description']="""Convert a UniFrac sample mapping file to an OTU table"""
script_info['script_description']="""This script allows users that have already created sample mapping (environment) files for use with the Unifrac web interface to use QIIME. QIIME records this data in an OTU table."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""Convert a UniFrac sample mapping (environment) file into a biom-formatted OTU table: ""","""%prog -i otu_table.sample_mapping.txt -o otu_table.biom"""))

script_info['output_description']="""The result of this script is an OTU table."""

script_info['required_options']=[
    make_option('-i', '--sample_mapping_fp',type='existing_filepath',help='path to the sample mapping file'),
    make_option('-o', '--output_fp',type='new_filepath',help='path to output file')
]

script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    sample_mapping_fp = opts.sample_mapping_fp
    output_fp = opts.output_fp
    verbose = opts.verbose
    
    sample_mapping_file = open(sample_mapping_fp, 'U')
    result = sample_mapping_to_biom_table(sample_mapping_file)
    open(output_fp, 'w').write(format_biom_table(result))

if __name__ == "__main__":
    main()

