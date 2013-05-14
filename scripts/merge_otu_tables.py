#!/usr/bin/env python
# File created on 30 Aug 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 

from qiime.util import make_option
from qiime.util import (parse_command_line_parameters, 
                        get_options_lookup)
from biom.parse import parse_biom_table
from qiime.format import format_biom_table

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Merge two or more OTU tables into a single OTU table."
script_info['script_description'] = """This script merges two or more OTU tables into a single OTU table. This is useful, for example, when you've created several reference-based OTU tables for different analyses and need to combine them for a larger analysis. 

Requirements: It is also very important that your OTUs are consistent across the different OTU tables. For example, you cannot safely merge OTU tables from two independent de novo OTU picking runs. Finally, either all or none of the OTU tables can contain taxonomic information: you can't merge some OTU tables with taxonomic data and some without taxonomic data."""
script_info['script_usage'] = [\
 ("",
  "Merge two OTU tables into a single OTU table",
  "%prog -i otu_table1.biom,otu_table2.biom -o merged_otu_table.biom")]
script_info['output_description']= ""
script_info['required_options'] = [\
 # Example required option
 make_option('-i','--input_fps',type='existing_filepaths',
             help='the otu tables in biom format (comma-separated)'),\
 make_option('-o','--output_fp',type='new_filepath',
             help='the output otu table filepath'),\
]
script_info['optional_options'] = []
script_info['version'] = __version__

def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    input_fps = opts.input_fps
    
    master = parse_biom_table(open(input_fps[0],'U'))
    for input_fp in input_fps[1:]:
        master = master.merge(parse_biom_table(open(input_fp,'U')))

    out_f = open(opts.output_fp,'w')
    out_f.write(format_biom_table(master))
    out_f.close()

if __name__ == "__main__":
    main()
