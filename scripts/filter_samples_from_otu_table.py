#!/usr/bin/env python
# File created on 21 Dec 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso",
               "Rob Knight",
               "Jesse Stombaugh",
               "Dan Knights",
               "Daniel McDonald",
               "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from itertools import izip
from numpy import inf, isinf
from biom.parse import parse_biom_table
from biom.util import biom_open

from qiime.util import (parse_command_line_parameters, make_option,
                        write_biom_table)
from qiime.parse import parse_mapping_file
from qiime.filter import (sample_ids_from_metadata_description,
                          filter_samples_from_otu_table,
                          filter_mapping_file)
from qiime.format import format_mapping_file

script_info = {}
script_info[
    'brief_description'] = "Filters samples from an OTU table on the basis of the number of observations in that sample, or on the basis of sample metadata. Mapping file can also be filtered to the resulting set of sample ids."
script_info['script_description'] = ""
script_info[
    'script_usage'] = [("Abundance filtering (low coverage)", "Filter samples with fewer than 150 observations from the otu table.", "%prog -i otu_table.biom -o otu_table_no_low_coverage_samples.biom -n 150"),
                       ("Abundance filtering (high coverage)", "Filter samples with greater than 149 observations from the otu table.",
                        "%prog -i otu_table.biom -o otu_table_no_high_coverage_samples.biom -x 149"),
                       ("Metadata-based filtering (positive)", "Filter samples from the table, keeping samples where the value for 'Treatment' in the mapping file is 'Control'",
                        "%prog -i otu_table.biom -o otu_table_control_only.biom -m map.txt -s 'Treatment:Control'"),
                       ("Metadata-based filtering (negative)", "Filter samples from the table, keeping samples where the value for 'Treatment' in the mapping file is not 'Control'",
                        "%prog -i otu_table.biom -o otu_table_not_control.biom -m map.txt -s 'Treatment:*,!Control'"),
                       ("List-based filtering", "Filter samples where the id is listed in samples_to_keep.txt", "%prog -i otu_table.biom -o otu_table_samples_to_keep.biom --sample_id_fp samples_to_keep.txt")]
script_info['output_description'] = ""
script_info['required_options'] = [
    make_option('-i', '--input_fp', type="existing_filepath",
                help='the input otu table filepath in biom format'),
    make_option('-o', '--output_fp', type="new_filepath",
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
    make_option('--sample_id_fp',
                type='existing_filepath',
                help='path to file listing sample ids to keep [default: %default]'),
    make_option('-s',
                '--valid_states', type='string',
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
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    input_fp = opts.input_fp
    output_fp = opts.output_fp

    mapping_fp = opts.mapping_fp
    output_mapping_fp = opts.output_mapping_fp
    valid_states = opts.valid_states
    min_count = opts.min_count
    max_count = opts.max_count
    sample_id_fp = opts.sample_id_fp

    if not ((mapping_fp and valid_states) or
            min_count != 0 or
            not isinf(max_count) or
            sample_id_fp is not None):
        option_parser.error("No filtering requested. Must provide either "
                            "mapping_fp and valid states, min counts, "
                            "max counts, or sample_id_fp (or some combination "
                            "of those).")
    if output_mapping_fp and not mapping_fp:
        option_parser.error("Must provide input mapping file to generate"
                            " output mapping file.")

    with biom_open(opts.input_fp) as biom_file:
        otu_table = parse_biom_table(biom_file)

    if mapping_fp and valid_states:
        sample_ids_to_keep = sample_ids_from_metadata_description(
            open(mapping_fp, 'U'), valid_states)
    else:
        sample_ids_to_keep = otu_table.sample_ids

    if sample_id_fp is not None:
        sample_id_f_ids = set([l.strip().split()[0]
                              for l in open(sample_id_fp, 'U') if not l.startswith('#')])
        sample_ids_to_keep = set(sample_ids_to_keep) & sample_id_f_ids

    filtered_otu_table = filter_samples_from_otu_table(otu_table,
                                                       sample_ids_to_keep,
                                                       min_count,
                                                       max_count)
    write_biom_table(filtered_otu_table, output_fp)

    # filter mapping file if requested
    if output_mapping_fp:
        mapping_data, mapping_headers, _ = parse_mapping_file(
            open(mapping_fp, 'U'))
        mapping_headers, mapping_data = \
            filter_mapping_file(
                mapping_data,
                mapping_headers,
                filtered_otu_table.sample_ids)
        open(
            output_mapping_fp,
            'w').write(
            format_mapping_file(
                mapping_headers,
                mapping_data))


if __name__ == "__main__":
    main()
