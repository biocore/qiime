#!/usr/bin/env python
# File created on 20 Dec 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from qiime.util import parse_command_line_parameters, make_option
from qiime.pycogent_backports.parse_biom import parse_otu_table_to_rich_otu_table, parse_biom_table

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","Convert an \"old style\" OTU table to biom format.","%prog -i otu_table.txt -o otu_table_biom.txt")]
script_info['output_description']= ""
script_info['required_options'] = [
 make_option('-i','--input_fp',type="existing_filepath",help='the input filepath'),
 make_option('-o','--output_fp',type="new_filepath",help='the output filepath'),
]

script_info['optional_options'] = [
  make_option('--dense',action='store_true',default=False,
   help="Write in dense biom format [default: %default]")
 ]
script_info['version'] = __version__

def convert_otu_table_to_biom(otu_table_f,count_map_f,dense):
    """ """
    otu_table = parse_otu_table_to_rich_otu_table(otu_table_f,
                                                  count_map_f=count_map_f,
                                                  dense=dense)
    return otu_table.getBiomFormatJsonString()

def convert_biom_to_otu_table(biom_f):
    """ """
    otu_table = parse_biom_table(biom_f)
    return otu_table.delimitedSelf(header_key='taxonomy', 
                                   header_value='Consensus Lineage',
                                   metadata_formatter=lambda x: ';'.join(x))


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    input_f = open(opts.input_fp,'U')
    output_f = open(opts.output_fp,'w')
    otu_table_to_biom = not opts.input_fp.endswith('.biom')
    dense = opts.dense
    count_map_f = int
    
    if otu_table_to_biom:
        output_f.write(convert_otu_table_to_biom(input_f,
                                                 count_map_f,
                                                 dense))
    else:
        output_f.write(convert_biom_to_otu_table(input_f))
    input_f.close()
    output_f.close()
    

if __name__ == "__main__":
    main()