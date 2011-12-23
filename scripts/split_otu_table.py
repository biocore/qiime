#!/usr/bin/env python
# File created on 15 Jun 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"
 

from os.path import split, splitext, join
from numpy import inf
from qiime.util import parse_command_line_parameters, make_option, create_dir
from qiime.parse import parse_mapping_file
from qiime.filter import filter_mapping_file, sample_ids_from_metadata_description
from qiime.format import format_mapping_file
from qiime.pycogent_backports.parse_biom import parse_biom_table

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
script_info['optional_options'] = [
 make_option('-c','--column_rename_ids',help='Mapping column used as sample id in the output files.' +\
                ' Has to be unique in the splited samples. This option can be helpful to create otu tables' +
                ' and mapping files for Procustes analysis.', default=None),
 make_option('--include_repeat_cols',action='store_true', help='By default the new mapping files' +\
                ' will not have the columns that have the same information, to include them use this' +\
                ' option. This can be helpful to create mapping files for Procrustes analysis.', 
                default=False),
]
script_info['version'] = __version__


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    otu_table_fp = opts.otu_table_fp
    mapping_fp = opts.mapping_fp
    mapping_field = opts.mapping_field
    output_dir = opts.output_dir
    column_rename_ids = opts.column_rename_ids
    include_repeat_cols = opts.include_repeat_cols
    
    otu_table_base_name = splitext(split(otu_table_fp)[1])[0]
    
    mapping_data, headers, comments = parse_mapping_file(open(mapping_fp,'U'))
    try:
        field_index = headers.index(mapping_field)
    except ValueError:
        option_parser.error("Field is not in mapping file (search is case "+\
        "and white-space sensitive). \n\tProvided field: "+\
        "%s. \n\tValid fields: %s" % (mapping_field,' '.join(headers)))
    if column_rename_ids: 
        try:
            column_rename_ids = headers.index(column_rename_ids)
        except ValueError:
            option_parser.error("Field is not in mapping file (search is case "+\
                 "and white-space sensitive). \n\tProvided field: "+\
                 "%s. \n\tValid fields: %s" % (mapping_field,' '.join(headers)))
    
    mapping_values = set([e[field_index] for e in mapping_data])
    
    create_dir(output_dir)
    
    mapping_data, mapping_headers, _ = parse_mapping_file(open(mapping_fp,'U'))
    otu_table = parse_biom_table(open(otu_table_fp,'U'))
    
    for v in mapping_values:
        v_fp_str = v.replace(' ','_')
        otu_table_output_fp = join(output_dir,'%s_%s.biom' % (otu_table_base_name, v_fp_str))
        mapping_output_fp = join(output_dir,'mapping_%s.txt' % v_fp_str)
        sample_ids_to_keep = sample_ids_from_metadata_description(
            open(mapping_fp,'U'),valid_states_str="%s:%s" % (mapping_field,v))
        
        # parse mapping file each time though the loop as filtering operates on values
        mapping_data, mapping_headers, _ = parse_mapping_file(open(mapping_fp,'U'))
        mapping_headers, mapping_data = filter_mapping_file(
                                         mapping_data, 
                                         mapping_headers,
                                         sample_ids_to_keep,
                                         include_repeat_cols=include_repeat_cols, 
                                         column_rename_ids=column_rename_ids)
        open(mapping_output_fp,'w').write(format_mapping_file(mapping_headers, mapping_data))
        filtered_otu_table = otu_table.filterSamples(
                              lambda values,id_,metadata: id_ in sample_ids_to_keep,
                              invert=True)
        open(otu_table_output_fp,'w').write(filtered_otu_table.getBiomFormatJsonString())
    


if __name__ == "__main__":
    main()