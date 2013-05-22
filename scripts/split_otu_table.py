#!/usr/bin/env python
# File created on 15 Jun 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 

from os.path import split, splitext, join
from numpy import inf
from qiime.util import parse_command_line_parameters, make_option, create_dir
from qiime.parse import parse_mapping_file
from qiime.split import split_mapping_file_on_field, split_otu_table_on_sample_metadata

script_info = {}
script_info['brief_description'] = "Split in a single OTU table into one OTU table per value in a specified field of the mapping file."
script_info['script_description'] = ""
script_info['script_usage'] = [("","Split otu_table.biom into per-study OTU tables, and store the results in ./per_study_otu_tables/","%prog -i otu_table.biom -m Fasting_Map.txt -f Treatment -o per_study_otu_tables")]
script_info['output_description']= ""
script_info['required_options'] = [
 make_option('-i','--otu_table_fp',type="existing_filepath",help='the input otu table'),
 make_option('-m','--mapping_fp',type="existing_filepath",help='the mapping file path'),
 make_option('-f','--mapping_field',type='string',help="mapping column to split otu table on"),
 make_option('-o','--output_dir',type="new_dirpath",help='the output directory'),
]
script_info['optional_options'] = [
 # this is known issue, see https://github.com/qiime/qiime/issues/417
 # and https://github.com/qiime/qiime/issues/941
 # make_option('-c','--column_rename_ids',type='string',help='Mapping column used as sample id in the output files.' +\
 #                ' Has to be unique in the splited samples. This option can be helpful to create otu tables' +
 #                ' and mapping files for Procustes analysis.', default=None),
 # make_option('--include_repeat_cols',action='store_true', help='By default the new mapping files' +\
 #                ' will not have the columns that have the same information, to include them use this' +\
 #                ' option. This can be helpful to create mapping files for Procrustes analysis.', 
 #                default=False)
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    otu_table_fp = opts.otu_table_fp
    mapping_fp = opts.mapping_fp
    mapping_field = opts.mapping_field
    output_dir = opts.output_dir
    # column_rename_ids = opts.column_rename_ids
    # include_repeat_cols = opts.include_repeat_cols
    
    create_dir(output_dir)
    
    # split mapping file
    mapping_f = open(mapping_fp,'U')
    for fp_str, sub_mapping_s in split_mapping_file_on_field(mapping_f,mapping_field):
        mapping_output_fp = join(output_dir,'mapping_%s.txt' % fp_str)
        open(mapping_output_fp,'w').write(sub_mapping_s)
    
    # split otu table
    otu_table_base_name = splitext(split(otu_table_fp)[1])[0]
    mapping_f = open(mapping_fp,'U')
    otu_table_f = open(otu_table_fp,'U')
    for fp_str, sub_otu_table_s in split_otu_table_on_sample_metadata(otu_table_f,
                                                                      mapping_f,
                                                                      mapping_field):
        otu_table_output_fp = join(output_dir,'%s_%s.biom' % (otu_table_base_name, fp_str))
        open(otu_table_output_fp,'w').write(sub_otu_table_s)


if __name__ == "__main__":
    main()
