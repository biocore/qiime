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
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from itertools import izip
from numpy import inf, isinf
from biom import load_table

from qiime.util import (parse_command_line_parameters, make_option,
                        write_biom_table, EmptyBIOMTableError)
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
                       ("ID-based filtering", "Keep samples where the id is listed in ids.txt", "%prog -i otu_table.biom -o filtered_otu_table.biom --sample_id_fp ids.txt"),
                       ("ID-based filtering (negation)", "Discard samples where the id is listed in ids.txt", "%prog -i otu_table.biom -o filtered_otu_table.biom --sample_id_fp ids.txt --negate_sample_id_fp")]

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
                help='Path to file listing sample ids to keep. Valid formats for the file are: 1) any white space, newline, or tab delimited list of samples, 2) a mapping file with samples in the first column [default: %default]'),
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
                help="the maximum total observation count in a sample for that sample to be retained [default: infinity]"),
    make_option('--negate_sample_id_fp',
                action='store_true', default=False,
                help='discard samples specified in --sample_id_fp instead of keeping them [default: %default]'),

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

    if (mapping_fp is None and valid_states is not None):
        option_parser.error("--mapping_fp must be provided if --valid_states "
                            "is passed.")

    if not ((mapping_fp and valid_states) or
            min_count != 0 or
            not isinf(max_count) or
            sample_id_fp is not None):
        option_parser.error("No filtering requested. Must provide either "
                            "mapping_fp and valid states, min counts, "
                            "max counts, or sample_id_fp (or some combination "
                            "of those).")
    if (mapping_fp and valid_states) and sample_id_fp:
        option_parser.error("Providing both --sample_id_fp and "
                            "--mapping_fp/--valid_states is not supported.")
    if output_mapping_fp and not mapping_fp:
        option_parser.error("Must provide input mapping file to generate"
                            " output mapping file.")

    otu_table = load_table(opts.input_fp)

    negate_sample_id_fp = opts.negate_sample_id_fp
    if mapping_fp and valid_states:
        sample_ids_to_keep = sample_ids_from_metadata_description(
            open(mapping_fp, 'U'), valid_states)
        negate_sample_id_fp = False
    else:
        sample_ids_to_keep = otu_table.ids()

        if sample_id_fp is not None:
            o = open(sample_id_fp, 'U')
            sample_id_f_ids = set([l.strip().split()[0] for l in o if not
                                   l.startswith('#')])
            o.close()
            sample_ids_to_keep = set(sample_ids_to_keep) & sample_id_f_ids

    filtered_otu_table = filter_samples_from_otu_table(
        otu_table, sample_ids_to_keep, min_count, max_count,
        negate_ids_to_keep=negate_sample_id_fp)

    try:
        write_biom_table(filtered_otu_table, output_fp)
    except EmptyBIOMTableError:
        option_parser.error(
            "Filtering resulted in an empty BIOM table. "
            "This indicates that no samples remained after filtering.")

    # filter mapping file if requested
    if output_mapping_fp:
        mapping_data, mapping_headers, _ = parse_mapping_file(
            open(mapping_fp, 'U'))
        mapping_headers, mapping_data = \
            filter_mapping_file(
                mapping_data,
                mapping_headers,
                filtered_otu_table.ids())
        open(
            output_mapping_fp,
            'w').write(
            format_mapping_file(
                mapping_headers,
                mapping_data))


if __name__ == "__main__":
    main()
