#!/usr/bin/env python
# File created on 21 Dec 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso",
               "Rob Knight",
               "Jesse Stombaugh",
               "Dan Knights",
               "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from itertools import izip
from numpy import inf, isinf
from qiime.pycogent_backports.parse_biom import parse_biom_table
from qiime.util import parse_command_line_parameters, make_option
from qiime.filter import (sample_ids_from_metadata_description, 
                          filter_samples_from_otu_table,
                          filter_mapping_file)

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","","")]
script_info['output_description']= ""
script_info['required_options'] = [
 make_option('-i','--input_fp',type="existing_filepath",
             help='the input otu table filepath in biom format'),
 make_option('-o','--output_fp',type="new_filepath",
             help='the output filepath in biom format'),
]
script_info['optional_options'] = [
 make_option('-m',
             '--mapping_fp',
             type='existing_filepath',
             help='path to the map file [default: %default]'),
 make_option('--output_mapping_fp',
             type='new_filepath',
             help='path to write filtered mapping file [default: filtered mapping file is not written]'),
 make_option('-s', 
             '--valid_states',
             help="string describing valid states (e.g. 'Treatment:Fasting') [default: %default]"),
 make_option('-n', 
             '--min_count',
             type='int',
             default=0,
             help="the minimum total observation count in a sample for that sample to be retained [default: %default]"),
 make_option('-x', 
             '--max_count',
             type='int',
             default=inf,
             help="the maximum total observation count in a sample for that sample to be retained [default: infinity]")

]
script_info['version'] = __version__


def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    input_fp = opts.input_fp
    output_fp = opts.output_fp
    
    mapping_fp = opts.mapping_fp
    output_mapping_fp = opts.output_mapping_fp
    valid_states = opts.valid_states
    min_count = opts.min_count
    max_count = opts.max_count
    
    if not ((mapping_fp and valid_states) or 
            min_count != 0 or 
            not isinf(max_count)):
        option_parser.error("No filtering requested. Must provide either "
                     "mapping_fp and valid states, min counts, or "
                     "max counts (or some combination of those).")
    if output_mapping_fp and not mapping_fp:
        option_parser.error("Must provide input mapping file to generate"
                            " output mapping file.")

    otu_table = parse_biom_table(open(opts.input_fp,'U'))
    output_f = open(opts.output_fp,'w')
    
    if (mapping_fp and valid_states):
        sample_ids_to_keep = sample_ids_from_metadata_description(
                              open(mapping_fp,'U'),valid_states)
    else:
        sample_ids_to_keep = otu_table.SampleIds
    
    filtered_otu_table = filter_samples_from_otu_table(otu_table,
                                                        sample_ids_to_keep,
                                                        min_count,
                                                        max_count)
    output_f.write(filtered_otu_table.getBiomFormatJsonString())
    output_f.close()
    
    # filter mapping file if requested
    if output_mapping_fp:
        filtered_mapping_str = \
         filter_mapping_file(open(mapping_fp,'U'),filtered_otu_table.SampleIds)
        open(output_mapping_fp,'w').write(filtered_mapping_str)
        


if __name__ == "__main__":
    main()