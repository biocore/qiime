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
from qiime.pycogent_backports.parse_biom import convert_otu_table_to_biom, convert_biom_to_otu_table, parse_biom_table

script_info = {}
script_info['brief_description'] = "Converts between classic otu table and biom formatted OTU tables."
script_info['script_description'] = ""
script_info['script_usage'] = [("","Convert a classic OTU table to biom format.","%prog -i otu_table.txt -o otu_table.biom"),
 ("","Convert a biom format otu table to a classic OTU table","%prog -i otu_table.biom -o otu_table.txt -b")]
script_info['output_description']= ""
script_info['required_options'] = [
 make_option('-i','--input_fp',type="existing_filepath",help='the input filepath'),
 make_option('-o','--output_fp',type="new_filepath",help='the output filepath'),
]

script_info['optional_options'] = [
  make_option('-t','--biom_type',type='choice',choices=['sparse','dense'],
   default='sparse',
   help="Type of biom file to write (dense or sparse) when passed a classic otu table [default: %default]"),
  make_option('-b','--biom_to_classic_otu_table',action='store_true',
   help="Convert biom file to classic otu table file [default: convert "
        "classic otu table file to biom file]",default=False),
  make_option('--sparse_biom_to_dense_biom',action='store_true',
   help="Convert sparse biom file to a dense biom file [default: convert "
        "classic otu table file to biom file]",default=False),
  make_option('--dense_biom_to_sparse_biom',action='store_true',
   help="Convert dense biom file to a sparse biom file [default: convert "
        "classic otu table file to biom file]",default=False),
  make_option('-m','--mapping_fp',type="existing_filepath",
         help='the mapping filepath (will add sample metadata to '+\
              'biom file, if provided) [default: %default]'),
 ]
script_info['version'] = __version__


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    biom_to_classic_otu_table = opts.biom_to_classic_otu_table
    sparse_biom_to_dense_biom = opts.sparse_biom_to_dense_biom
    dense_biom_to_sparse_biom = opts.dense_biom_to_sparse_biom
    
    if sum([biom_to_classic_otu_table,
            sparse_biom_to_dense_biom,
            dense_biom_to_sparse_biom]) > 1:
        option_parser.error("The --biom_to_classic_otu_table, --sparse_biom_to_dense_biom, "
         "and --dense_biom_to_sparse_biom options are mutually exclusive. Pass only one at a time.")
    
    input_f = open(opts.input_fp,'U')
    output_f = open(opts.output_fp,'w')
    
    dense = opts.biom_type == 'dense'
    count_map_f = int
    mapping_fp = opts.mapping_fp
    
    if mapping_fp != None:
        mapping_f = open(mapping_fp,'U')
    else:
        mapping_f = None
    
    if biom_to_classic_otu_table:
        try:
            output_f.write(convert_biom_to_otu_table(input_f))
        except ValueError:
            raise ValueError, "Input does not look like a .biom file. Did you accidentally specify -b?"
    elif sparse_biom_to_dense_biom:
        try:
            otu_table = parse_biom_table(input_f,dense_object=True)
        except ValueError:
            raise ValueError, "Input does not look like a .biom file. Did you accidentally specify -b?"        
        output_f.write(otu_table.getBiomFormatJsonString())
    elif dense_biom_to_sparse_biom:
        try:
            otu_table = parse_biom_table(input_f,dense_object=False)
        except ValueError:
            raise ValueError, "Input does not look like a .biom file. Did you accidentally specify -b?"
        output_f.write(otu_table.getBiomFormatJsonString())
    else:
        try:
            output_f.write(convert_otu_table_to_biom(input_f,
                                                     count_map_f,
                                                     dense,
                                                     mapping_f))
        except ValueError:
            raise ValueError, "Input does not look like a classic OTU table. Do you need to pass -b?"
    input_f.close()
    output_f.close()
    

if __name__ == "__main__":
    main()