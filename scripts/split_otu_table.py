#!/usr/bin/env python
# File created on 15 Jun 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"
 

from os.path import split, splitext, join
from qiime.util import parse_command_line_parameters, make_option, create_dir
from qiime.parse import parse_mapping_file
from qiime.filter_by_metadata import filter_otus_and_map

script_info = {}
script_info['brief_description'] = "Split in a single OTU table into one OTU table per value in a specified field of the mapping file."
script_info['script_description'] = ""
script_info['script_usage'] = [("","Split otu_table.txt into per-study OTU tables, and store the results in ./per_study_otu_tables/","%prog -i otu_table.txt -m mapping.txt -f Study -o per_study_otu_tables")]
script_info['output_description']= ""
script_info['required_options'] = [
 make_option('-i','--otu_table_fp',type="existing_filepath",help='the input otu table'),
 make_option('-m','--mapping_fp',type="existing_filepath",help='the mapping file path'),
 make_option('-f','--mapping_field',help="mapping column to split otu table on"),
 make_option('-o','--output_dir',type="new_dirpath",help='the output directory'),
]
script_info['optional_options'] = []
script_info['version'] = __version__


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    otu_table_fp = opts.otu_table_fp
    mapping_fp = opts.mapping_fp
    mapping_field = opts.mapping_field
    output_dir = opts.output_dir
    
    otu_table_base_name = splitext(split(otu_table_fp)[1])[0]
    
    mapping_data, headers, comments = parse_mapping_file(open(mapping_fp,'U'))
    try:
        field_index = headers.index(mapping_field)
    except ValueError:
        option_parser.error("Field is not in mapping file (search is case "+\
        "and white-space sensitive). \n\tProvided field: "+\
        "%s. \n\tValid fields: %s" % (mapping_field,' '.join(headers)))
    
    mapping_values = set([e[field_index] for e in mapping_data])
    
    create_dir(output_dir)
    
    for v in mapping_values:
        v_fp_str = v.replace(' ','_')
        otu_table_output_fp = join(output_dir,'%s_%s.txt' % (otu_table_base_name, v_fp_str))
        mapping_output_fp = join(output_dir,'mapping_%s.txt' % v_fp_str)
        filter_otus_and_map(open(mapping_fp,'U'), 
                            open(otu_table_fp,'U'), 
                            open(mapping_output_fp,'w'), 
                            open(otu_table_output_fp,'w'),
                            valid_states_str="%s:%s" % (mapping_field,v),
                            num_seqs_per_otu=1)
    


if __name__ == "__main__":
    main()