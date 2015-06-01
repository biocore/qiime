#!/usr/bin/env python
# File created on 20 Dec 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


from biom import load_table

from qiime.util import (parse_command_line_parameters,
                        get_options_lookup, make_option)
from qiime.filter import (filter_samples_from_distance_matrix,
                          sample_ids_from_metadata_description,
                          get_seqs_to_keep_lookup_from_seq_id_file)
from qiime.parse import parse_distmat

options_lookup = get_options_lookup()

script_info = {}
script_info[
    'brief_description'] = "Filter a distance matrix to contain only a specified set of samples."
script_info[
    'script_description'] = "Remove samples from a distance matrix based on a mapping file or an otu table or a list of sample ids."
script_info['script_usage'] = [("",
                                "Filter samples ids listed in sample_id_list.txt from dm.txt",
                                "%prog -i dm.txt -o dm_out_sample_list.txt --sample_id_fp sample_id_list.txt"),
                               ("",
                                "Filter samples ids in otu_table.biom from dm.txt",
                                "%prog -i dm.txt -o dm_out_otu_table.txt -t otu_table.biom"),
                               ("",
                                "Filter samples ids where DOB is 20061218 in Fasting_Map.txt. (Run \"filter_samples_from_otu_table.py -h\" for additional information on how metadata filtering can be specified.)",
                                "%prog -i dm.txt -o dm_out_mapping_file.txt -m Fasting_Map.txt -s \"DOB:20061218\""), ]
script_info['output_description'] = ""
script_info['required_options'] = [
    make_option('-i', '--input_distance_matrix',
                help='the input distance matrix', type="existing_filepath"),
    make_option('-o', '--output_distance_matrix',
                help='path to store the output distance matrix', type="new_filepath")]

script_info['optional_options'] = [
    make_option('--sample_id_fp', type="existing_filepath",
                help='A list of sample identifiers (or tab-delimited lines with'
                ' a sample identifier in the first field) which should be retained'),
    make_option('-t', '--otu_table_fp',
                type="existing_filepath", help='the otu table filepath'),
    make_option(
        '-m',
        '--mapping_fp',
        help='path to the mapping file',
        type="existing_filepath"),
    make_option(
        '-s',
        '--valid_states',
        type='string',
        help="string containing valid states, e.g. 'STUDY_NAME:DOB'"),
    make_option('--negate', default=False,
                action='store_true',
                help="discard specified samples (instead of keeping them) [default: %default]")]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    output_f = open(opts.output_distance_matrix, 'w')
    if opts.otu_table_fp:
        otu_table = load_table(opts.otu_table_fp)
        samples_to_keep = otu_table.ids()
        # samples_to_keep = \
        # sample_ids_from_otu_table(open(opts.otu_table_fp,'U'))
    elif opts.sample_id_fp:
        samples_to_keep = \
            get_seqs_to_keep_lookup_from_seq_id_file(
                open(opts.sample_id_fp, 'U'))
    elif opts.mapping_fp and opts.valid_states:
        try:
            samples_to_keep = sample_ids_from_metadata_description(
                open(opts.mapping_fp, 'U'), opts.valid_states)
        except ValueError as e:
            option_parser.error(e.message)
    else:
        option_parser.error('must pass either --sample_id_fp, -t, or -m and '
                            '-s')
    # note that negate gets a little weird here. The function we're calling
    # removes the specified samples from the distance matrix, but the other
    # QIIME filter scripts keep these samples specified.  So, the interface of
    # this script is designed to keep the specified samples, and therefore
    # negate=True is passed to filter_samples_from_distance_matrix by default.
    d = filter_samples_from_distance_matrix(
        parse_distmat(
            open(opts.input_distance_matrix, 'U')),
        samples_to_keep,
        negate=not opts.negate)
    output_f.write(d)
    output_f.close()


if __name__ == "__main__":
    main()
