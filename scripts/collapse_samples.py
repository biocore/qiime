#!/usr/bin/env python
# File created on 08 Oct 2014
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

import numpy as np
from biom import load_table
from biom.exception import UnknownAxisError

from qiime.util import (parse_command_line_parameters, make_option,
                         write_biom_table)
from qiime.group import (collapse_samples, get_collapse_fns,
                          mapping_lines_from_collapsed_df)

collapse_fns = get_collapse_fns()
collapse_modes = collapse_fns.keys()

script_info = {}
script_info['brief_description'] = "Collapse samples in a BIOM table and mapping file."
script_info['script_description'] = ("Collapse samples in a BIOM table and mapping file. "
    "Values in the BIOM table are collapsed in one of several different ways; see the "
    "available options for --collapse_mode. Values in the mapping file are collapsed "
    "by grouping the values if they differ for the grouped samples, and by providing the single "
    "value if they don't differ for the grouped samples.")
script_info['script_usage'] = [
    ("Collapse samples in biom table and mapping file",
     "Collapse samples by taking the median value for each observation in each group, where group is defined by having the same values for subject in the mapping file.",
     "%prog -b table.biom -m map.txt --output_biom_fp collapsed.biom --output_mapping_fp collapsed_map.txt --collapse_mode median --collapse_fields subject"),
    ("Collapse samples in biom table and mapping file",
     "Collapse samples by taking the median value for each observation in each group, where group is defined by having the same values for both subject and replicate-group in the mapping file.",
     "%prog -b table.biom -m map.txt --output_biom_fp collapsed.biom --output_mapping_fp collapsed_map.txt --collapse_mode median --collapse_fields replicate-group,subject"),
    ("Collapse samples in biom table and mapping file, and normalize the table",
     "Collapse samples by taking the median value for each observation in each group, where group is defined by having the same values for both subject and replicate-group in the mapping file. Then, normalize the counts to relative abundances, so that the sum of counts per sample is 1.0.",
     "%prog -b table.biom -m map.txt --output_biom_fp collapsed-normed.biom --output_mapping_fp collapsed_map.txt --collapse_mode median --collapse_fields replicate-group,subject --normalize")
     ]
script_info['output_description']= "A collapsed mapping file and BIOM table will be generated at the requested paths."
script_info['required_options'] = [
    make_option('-b', '--input_biom_fp', type='existing_filepath',
                help='the biom table containing the samples to be collapsed'),
    make_option('-m', '--mapping_fp', type='existing_filepath',
                help='the sample metdata mapping file'),
    make_option('--output_biom_fp', type="new_filepath",
                help="path where collapsed biom table should be written"),
    make_option('--output_mapping_fp', type="new_filepath",
                help="path where collapsed mapping file should be written"),
    make_option('--collapse_fields',
                help="comma-separated list of fields to collapse on")
]
script_info['optional_options'] = [
    make_option('--collapse_mode', type='choice', choices=collapse_modes,
        help="the mechanism for collapsing counts within groups; "
        "valid options are: %s. " % ', '.join(collapse_modes) + 
        "[default: %default]", default='sum'),
    make_option('--normalize',
        help='Normalize observation counts to relative abundances, so the '
             'counts within each sample sum to 1.0. [default: %default]',
             default=False, action='store_true')]

script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    mapping_fp = opts.mapping_fp
    collapse_fields = opts.collapse_fields.split(',')
    input_biom_fp = opts.input_biom_fp
    collapse_mode = opts.collapse_mode
    output_biom_fp = opts.output_biom_fp
    output_mapping_fp = opts.output_mapping_fp
    normalize = opts.normalize

    collapsed_metadata, collapsed_table = \
        collapse_samples(load_table(input_biom_fp),
                         open(mapping_fp, 'U'),
                         collapse_fields,
                         collapse_mode)

    if normalize:
        collapsed_table.norm(axis='sample', inplace=True)

    write_biom_table(collapsed_table, output_biom_fp)
    output_map_lines = mapping_lines_from_collapsed_df(collapsed_metadata)
    with open(output_mapping_fp, 'w') as output_mapping_f:
        output_mapping_f.write('\n'.join(output_map_lines))

if __name__ == "__main__":
    main()
