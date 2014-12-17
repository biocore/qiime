#!/usr/bin/env python
# File created on 08 Oct 2014
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from functools import partial

import numpy as np
from biom import load_table
from biom.exception import UnknownAxisError

from qiime.util import (parse_command_line_parameters, make_option,
                         write_biom_table)
from qiime.group import (group_by_sample_metadata, collapse_metadata,
                          _collapse_to_median, _collapse_to_first,
                          _collapse_to_random, _partition_f,
                          _mapping_lines_from_collapsed_df)

collapse_fns = {'median': collapse_to_median,
                'first': collapse_to_first,
                'random': collapse_to_random}
collapse_modes = ['sum', 'mean'] + collapse_fns.keys()

script_info = {}
script_info['brief_description'] = "Collapse samples in a BIOM table and mapping file."
script_info['script_description'] = ""
script_info['script_usage'] = [
    ("Collapse samples in biom table and mapping file",
     "Collapse samples by taking the median value for each observation in each group, where group is defined by having the same values for both subject and replicate-group in the mapping file.",
     "%prog -b table.biom -m map.txt --output_biom_fp collapsed.biom --output_mapping_fp collapsed_map.txt --collapse_mode median --collapse_fields replicate-group,subject")]
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
    make_option('--collapse_mode', type='choice', choices=collapse_modes,
                help="the mechanism for collapsing counts from replicates; "
                     "valid options are: %s" % ' '.join(collapse_modes)),
    make_option('--collapse_fields',
                help="comma-separated list of fields to collapse on")
]
script_info['optional_options'] = []

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    mapping_fp = opts.mapping_fp
    collapse_fields = opts.collapse_fields
    input_biom_fp = opts.input_biom_fp
    collapse_mode = opts.collapse_mode
    output_biom_fp = opts.output_biom_fp
    output_mapping_fp = opts.output_mapping_fp

    collapsed_metadata = collapse_metadata(open(mapping_fp, 'U'),
                                           collapse_fields.split(','))

    new_index_to_group, old_index_to_new_index = \
        group_by_sample_metadata(open(mapping_fp, 'U'),
                                 collapse_fields.split(','))
    partition_f = partial(_partition_f, sid_to_group_id=old_index_to_new_index)
    input_table = load_table(input_biom_fp)

    if collapse_mode == 'sum':
        output_table = input_table.collapse(
            partition_f, norm=False, axis='sample')
    elif collapse_mode == 'mean':
        output_table = input_table.collapse(
            partition_f, norm=True, axis='sample')
    else:
        # A KeyError is not possible here, since a collapse_mode that is not
        # in collapse_fns or one of the modes listed above would be caught
        # by the option parser.
        collapse_f = collapse_fns[collapse_mode]
        output_table = input_table.collapse(
            partition_f, collapse_f=collapse_f, norm=False, axis='sample')

    write_biom_table(output_table, output_biom_fp)
    output_mapping_f = open(output_mapping_fp, 'w')
    output_map_lines = _mapping_lines_from_collapsed_df(collapsed_metadata)
    output_mapping_f.write('\n'.join(output_map_lines))
    output_mapping_f.close()

if __name__ == "__main__":
    main()
