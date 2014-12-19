#!/usr/bin/env python
# File created on 09 Feb 2010
# summarize_otu_by_cat.py
from __future__ import division

__author__ = "Julia Goodrich"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Julia Goodrich", "Greg Caporaso", "Justin Kuczynski",
               "Jose Carlos Clemente Litran", "Yoshiki Vazquez Baeza"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "wasade@gmail.com"

from os import getcwd, makedirs

from biom.parse import parse_biom_table
from biom.util import biom_open

from qiime.parse import (parse_mapping_file_to_dict, parse_mapping_file,
                         mapping_file_to_dict)
from qiime.util import (parse_command_line_parameters, make_option,
                        write_biom_table)
from qiime.colors import combine_map_label_cols, get_map

script_info = {}
script_info['brief_description'] = ("Summarize an OTU table by a single column "
                                    "in the mapping file.")
script_info['script_description'] = ("Collapse an OTU table based on values in a"
                                     " single column in the mapping file. For example, if you have 10 samples, "
                                     "five of which are from females and five of which are from males, you could"
                                     " use this script to collapse the ten samples into two corresponding based "
                                     "on their values in a 'Sex' column in your mapping file.")
script_info['script_usage'] = []
script_info['script_usage'].append(("Example:", " Collapsed otu_table.biom "
                                    "on the 'Treatment' column in Fasting_Map.txt and write the resulting OTU "
                                    "table to otu_table_by_treatment.biom", "%prog -i otu_table.biom -m "
                                    "Fasting_Map.txt -c Treatment -o otu_table_by_treatment.biom"))
script_info['script_usage'].append(("", " Combine two categories and collapse "
                                    "otu_table.biom on the 'Sex' and 'Age' columns in map.txt and write the "
                                    "resulting OTU table to otu_table_by_sex_and_age.txt", "%prog -i "
                                    "otu_table.biom -m Fasting_Map.txt -c 'Treatment&&DOB' -o "
                                    "otu_table_by_sex_and_age.txt"))
script_info['script_usage'].append(("Use normalization",
                                    "Normalize, combine two categories and"
                                    " collapse otu_table.biom on the "
                                    "'Treatment' and 'DOB' columns in "
                                    "Fasting_Map.txt and write the resulting "
                                    "OTU table to normalized_otu_table_by_treat"
                                    "ment_and_dob.txt",
                                    "%prog -i otu_table.biom -m Fasting_Map.txt"
                                    " -c 'Treatment&&DOB' -n -o normalized_otu_"
                                    "table_by_treatment_and_dob.txt"))

script_info['output_description'] = """"""
script_info['required_options'] = [
    make_option('-m', '--mapping_fp',
                help='Input metadata mapping filepath [REQUIRED]',
                type='existing_filepath'),
    make_option('-i', '--otu_table_fp',
                help='Input OTU table filepath. [REQUIRED]',
                type='existing_filepath'),
    make_option('-c', '--mapping_category', type='string',
                help='Summarize OTU table using this category. The user can '
                'also combine columns in the mapping file by separating the '
                'categories by "&&" without spaces. [REQUIRED]'),
    make_option('-o', '--output_fp', dest='output_fp',
                help='Output OTU table filepath. [REQUIRED]',
                type='new_filepath'),
]
script_info['optional_options'] = [
    make_option('-n', '--normalize',
                help='Transform the OTU counts to relative abundance.',
                default=False, action='store_true')
]
script_info['option_label'] = {'otu_table_fp': 'OTU table filepath',
                               'output_fp': 'Output filepath',
                               'mapping_fp':
                               'QIIME-formatted mapping filepath',
                               'mapping_category': 'Summarize category',
                               'normalize': 'Normalize counts'}

script_info["version"] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    mapping_fp = opts.mapping_fp
    mapping_category = opts.mapping_category
    otu_table_fp = opts.otu_table_fp
    output_fp = opts.output_fp
    normalize = opts.normalize

    # define a function that returns the bin a sample shouldbe placed into
    bin_function = lambda id_, sample_metadata:\
        sample_metadata[mapping_category]
    # parse the sample metadata and add it to the OTU table (we assume that
    # sample metadata is not already present in the table)
    mapping, headers, comments = parse_mapping_file(open(mapping_fp, 'U'))

    # added in ability to combine metadata columns and summarize based on the
    # new combined category
    if '&&' in mapping_category:
        new_mapping = []
        new_mapping.append(headers)
        for i in range(len(mapping)):
            new_mapping.append(mapping[i])
        # Create an array using multiple columns from mapping file
        combinecolorby = mapping_category.split('&&')
        mapping = combine_map_label_cols(combinecolorby, new_mapping)

    sample_metadata = mapping_file_to_dict(mapping, headers)
    with biom_open(otu_table_fp, 'U') as biom_file:
        table = parse_biom_table(biom_file)
    table.add_metadata(sample_metadata)
    # create a new OTU table where samples are binned based on their return
    # value from bin_function
    result = table.collapse(bin_function, norm=False, min_group_size=1,
                            axis='sample')

    # normalize the result if requested by the user
    if normalize:
        result.norm(axis='sample', inplace=True)

    # write a new BIOM file
    write_biom_table(result, output_fp)

if __name__ == "__main__":
    main()
