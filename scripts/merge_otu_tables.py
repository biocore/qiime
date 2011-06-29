#!/usr/bin/env python
# File created on 30 Aug 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"
 

from qiime.util import make_option
from qiime.util import (parse_command_line_parameters, 
                        get_options_lookup,
                        merge_n_otu_tables)
from qiime.format import format_otu_table

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Merge two or more OTU tables into a single OTU table."
script_info['script_description'] = """This script merges two or more OTU tables into a single OTU table. This is useful, for example, when you've created several reference-based OTU tables for different analyses and need to combine them for a larger analysis. 

Requirements: This process requires that the sample IDs are distinct across the different OTU tables. It is also very important that your OTUs are consistent across across the different OTU tables. For example, you cannot safely merge OTU tables from two independent de novo OTU picking runs. Finally, either all or none of the OTU tables can contain taxonomic information: you can't merge some OTU tables with taxonomic data and some without taxonomic data."""
script_info['script_usage'] = [\
 ("",
  "Merge two OTU tables into a single OTU table",
  "%prog -i otu_table1.txt,otu_table2.txt -o merged_otu_table.txt")]
script_info['output_description']= ""
script_info['required_options'] = [\
 # Example required option
 make_option('-i','--input_fps',help='the otu tables (comma-separated)'),\
 make_option('-o','--output_fp',help='the output otu table filepath'),\
]
script_info['optional_options'] = []
script_info['version'] = __version__

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
       
    input_fs = []
    for input_fp in opts.input_fps.split(','):
        input_fs.append(open(input_fp,'U'))
    
    sample_names, otu_names, data, taxonomy = merge_n_otu_tables(input_fs)
    
    out_f = open(opts.output_fp,'w')
    out_f.write(format_otu_table(sample_names=sample_names,
                                 otu_names=otu_names,
                                 data=data,
                                 taxonomy=taxonomy))
    out_f.close()






if __name__ == "__main__":
    main()