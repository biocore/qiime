#!/usr/bin/env python
# File created on 15 Feb 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Daniel McDonald", "Emily TerAvest",
               "Yoshiki Vazquez Baeza"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from biom import load_table

from qiime.parse import parse_mapping_file
from qiime.util import (parse_command_line_parameters, get_options_lookup,
                        make_option, write_biom_table)
from qiime.sort import (sort_otu_table, sort_otu_table_by_mapping_field,
                        natsort_case_insensitive)

options_lookup = get_options_lookup()

script_info = {}
script_info[
    'brief_description'] = "Script for sorting the sample IDs in an OTU table based on a specified value in a mapping file."
script_info['script_description'] = ""
script_info['script_usage'] = [("Default",
                                "case insensitive natural sorting"
                                " i.e. SAMP0, samp1, SAMP2, samp10, samp12",
                                "%prog -i otu_table.biom -o sorted_otu_table.biom"),
                               ("",
                                "sort samples by the age field in the mapping file",
                                "%prog -i otu_table.biom -o dob_sorted_otu_table.biom -m Fasting_Map.txt -s DOB"),
                               ("",
                                "sort samples based on order in a file where each line starts with a sample id",
                                "%prog -i otu_table.biom -o sorted_otu_table.biom -l sample_id_list.txt")]
script_info['output_description'] = ""
script_info['required_options'] = [
    make_option('-i', '--input_otu_table',
                help='Input OTU table filepath in BIOM format.',
                type='existing_filepath'),
    make_option('-o', '--output_fp',
                help='Output OTU table filepath.',
                type='new_filepath'),
]
script_info['optional_options'] = [
    make_option('-m', '--mapping_fp',
                help='Input metadata mapping filepath. [default: %default]',
                type='existing_filepath'),
    make_option('-s', '--sort_field', type='string',
                help='Category to sort OTU table by. [default: %default]'),
    make_option('-l', '--sorted_sample_ids_fp',
                help='Sorted sample id filepath [default: %default]',
                type='existing_filepath')
]
script_info[
    'option_label'] = {'input_otu_table': 'OTU table filepath in BIOM format',
                       'output_fp': 'Output filepath',
                       'mapping_fp':
                       'QIIME-formatted mapping filepath',
                       'sort_field': 'Category to sort by',
                       'sorted_sample_ids_fp': 'Sorted sample id filepath'}
script_info['version'] = __version__


def sample_ids_from_f(lines):
    result = []
    for line in lines:
        line = line.strip()
        if line and not line.startswith('#'):
            result.append(line.split()[0])
    return result


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    otu_table_data = load_table(opts.input_otu_table)
    sort_field = opts.sort_field
    mapping_fp = opts.mapping_fp
    sorted_sample_ids_fp = opts.sorted_sample_ids_fp

    if sort_field and mapping_fp:
        mapping_data = parse_mapping_file(open(mapping_fp, 'U'))
        result = sort_otu_table_by_mapping_field(otu_table_data, mapping_data,
                                                 sort_field)
    elif sorted_sample_ids_fp:
        sorted_sample_ids = sample_ids_from_f(open(sorted_sample_ids_fp, 'U'))
        result = sort_otu_table(otu_table_data,
                                sorted_sample_ids)
    else:
        result = sort_otu_table(otu_table_data,
            natsort_case_insensitive(otu_table_data.ids()))

    write_biom_table(result, opts.output_fp)

if __name__ == "__main__":
    main()
