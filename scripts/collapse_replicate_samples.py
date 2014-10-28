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
from qiime.group import group_by_sample_metadata, collapse_metadata


collapse_modes = ['sum', 'mean', 'median', 'first', 'random']

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
# Members of the tuple in script_usage are (title, description, example call)
script_info['script_usage'] = [
    ("Collapse samples in biom table and mapping file",
     "Collapse samples by taking the median value for each observation in each group, where group is defined by having the same values for both subject and replicate-group in the mapping file.",
     "%prog -b table.biom -m map.txt --output_biom_fp collapsed.biom --output_mapping_fp collapsed_map.txt --collapse_mode median --collapse_fields replicate-group,subject")]
script_info['output_description']= ""
script_info['required_options'] = [
    # Example required option
    make_option('-b', '--input_biom_fp', type='existing_filepath',
                help='the biom table containing the samples to be collapsed'),
    make_option('-m', '--mapping_fp', type='existing_filepath',
                help='the sample metdata mapping file'),
    make_option('--output_biom_fp', type="new_filepath",
                help="path where collapsed biom table should be written"),
    make_option('--output_mapping_fp', type="new_filepath",
                help="path where collapsed mapping file should be written"),
    make_option('--collapse_mode', type='choice', choices=collapse_modes,
                help="the mechanism for collapsing counts from replicates"),
    make_option('--collapse_fields',
                help="comma-separated list of fields to collapse on")
]
script_info['optional_options'] = [

]
script_info['version'] = __version__

def _partition_f(id_, md, sid_to_group_id):
    return '.'.join(map(str, sid_to_group_id[id_]))

def collapse_to_first(t, axis):
    return np.asarray([e[0] for e in t.iter_data(axis=axis, dense=True)])

def collapse_to_median(t, axis):
    return np.asarray([np.median(e) for e in t.iter_data(axis=axis, dense=True)])

def collapse_to_random(t, axis):
    # this is a little clunky - waiting for an answer here
    # http://stackoverflow.com/q/26050412/3424666
    if axis == 'sample':
        length = t.shape[0]
    elif axis == 'observation':
        length = t.shape[1]
    else:
        raise UnknownAxisError(axis)
    n = np.random.randint(length)
    return np.asarray([e[n] for e in t.iter_data(axis=axis, dense=True)])

def mapping_lines_from_collapsed_df(collapsed_df):
    lines = []
    lines.append('\t'.join(['#SampleID', 'represented-sample-id'] +\
                           list(collapsed_df.columns)[1:]))

    for r in collapsed_df.iterrows():
        new_idx = '.'.join(map(str, r[0]))
        new_values = map(str,[e[0] for e in r[1]])
        lines.append('\t'.join([new_idx] + new_values))
    return lines


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    collapsed_metadata = collapse_metadata(open(opts.mapping_fp, 'U'),
                                           opts.collapse_fields.split(','))

    new_index_to_group, old_index_to_new_index = \
        group_by_sample_metadata(open(opts.mapping_fp, 'U'),
                                 opts.collapse_fields.split(','))
    partition_f = partial(_partition_f, sid_to_group_id=old_index_to_new_index)
    input_table = load_table(opts.input_biom_fp)

    collapse_mode = opts.collapse_mode
    if collapse_mode == 'sum':
        output_table = input_table.collapse(
            partition_f, norm=False, axis='sample')
    elif collapse_mode == 'mean':
        output_table = input_table.collapse(
            partition_f, norm=True, axis='sample')
    elif collapse_mode == 'median':
        output_table = input_table.collapse(
            partition_f, collapse_f=collapse_to_median, norm=False,
            axis='sample')
    elif collapse_mode == 'first':
        output_table = input_table.collapse(
            partition_f, collapse_f=collapse_to_first, norm=False,
            axis='sample')
    elif collapse_mode == 'random':
        output_table = input_table.collapse(
            partition_f, collapse_f=collapse_to_random, norm=False,
            axis='sample')
    else:
        option_parser.error("Unknown collapse_mode (%s). Valid choices are: %s"
                            % (collapse_mode, ', '.join(collapse_modes)))

    write_biom_table(output_table, opts.output_biom_fp)
    output_mapping_f = open(opts.output_mapping_fp, 'w')
    output_map_lines = mapping_lines_from_collapsed_df(collapsed_metadata)
    output_mapping_f.write('\n'.join(output_map_lines))
    output_mapping_f.close()

if __name__ == "__main__":
    main()
