#!/usr/bin/env python
# File created on 24 Feb 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from skbio.parse.sequences import parse_fasta
from skbio.util.misc import create_dir
from qiime.parse import parse_mapping_file
from qiime.filter import filter_mapping_file, sample_ids_from_metadata_description
from qiime.format import format_mapping_file
from biom.table import TableException


def get_mapping_values(mapping_f, mapping_field):
    mapping_data, mapping_headers, comments = parse_mapping_file(mapping_f)

    try:
        field_index = mapping_headers.index(mapping_field)
    except ValueError:
        raise KeyError("Field is not in mapping file (search is case " +
                       "and white-space sensitive). \n\tProvided field: " +
                       "%s. \n\tValid fields: %s" % (mapping_field, ' '.join(mapping_headers)))
    return set([e[field_index] for e in mapping_data])


def split_mapping_file_on_field(mapping_f,
                                mapping_field,
                                column_rename_ids=None,
                                include_repeat_cols=True):
    """ split mapping file based on value in field """

    mapping_f = list(mapping_f)
    mapping_values = get_mapping_values(mapping_f, mapping_field)

    mapping_data, mapping_headers, _ = parse_mapping_file(mapping_f)

    if column_rename_ids:
        try:
            column_rename_ids = mapping_headers.index(column_rename_ids)
        except ValueError:
            raise KeyError("Field is not in mapping file (search is case " +
                           "and white-space sensitive). \n\tProvided field: " +
                           "%s. \n\tValid fields: %s" % (mapping_field, ' '.join(mapping_headers)))

    for v in mapping_values:
        v_fp_str = v.replace(' ', '_')
        sample_ids_to_keep = sample_ids_from_metadata_description(
            mapping_f, valid_states_str="%s:%s" % (mapping_field, v))

        # parse mapping file each time though the loop as filtering operates on
        # values
        mapping_data, mapping_headers, _ = parse_mapping_file(mapping_f)
        mapping_headers, mapping_data = filter_mapping_file(
            mapping_data,
            mapping_headers,
            sample_ids_to_keep,
            include_repeat_cols=include_repeat_cols,

            column_rename_ids=column_rename_ids)
        yield v_fp_str, format_mapping_file(mapping_headers, mapping_data)


def split_otu_table_on_sample_metadata(otu_table, mapping_f, mapping_field):
    """ split otu table into sub otu tables where each represent samples corresponding to only a certain value in mapping_field
    """
    mapping_f = list(mapping_f)
    mapping_values = get_mapping_values(mapping_f, mapping_field)

    for v in mapping_values:
        v_fp_str = v.replace(' ', '_')
        sample_ids_to_keep = sample_ids_from_metadata_description(
            mapping_f, valid_states_str="%s:%s" % (mapping_field, v))

        try:
            # filtering cannot be inplace otherwise we lose data
            filtered_otu_table = otu_table.filter(
                lambda values, id_, metadata: id_ in sample_ids_to_keep,
                axis='observation', inplace=False)
        except TableException:
            # all samples are filtered out, so no otu table to write
            continue
        yield v_fp_str, filtered_otu_table


def split_fasta(infile, seqs_per_file, outfile_prefix, working_dir=''):
    """ Split infile into files with seqs_per_file sequences in each

        infile: list of fasta lines or open file object
        seqs_per_file: the number of sequences to include in each file
        out_fileprefix: string used to create output filepath - output filepaths
         are <out_prefix>.<i>.fasta where i runs from 0 to number of output files
        working_dir: directory to prepend to temp filepaths (defaults to
         empty string -- files written to cwd)

        List of output filepaths is returned.

    """
    if seqs_per_file <= 0:
        raise ValueError("seqs_per_file must be > 0!")

    seq_counter = 0
    out_files = []
    if working_dir and not working_dir.endswith('/'):
        working_dir += '/'
        create_dir(working_dir)

    for seq_id, seq in parse_fasta(infile):
        if seq_counter == 0:
            current_out_fp = '%s%s.%d.fasta' \
                % (working_dir, outfile_prefix, len(out_files))
            current_out_file = open(current_out_fp, 'w')
            out_files.append(current_out_fp)
        current_out_file.write('>%s\n%s\n' % (seq_id, seq))
        seq_counter += 1

        if seq_counter == seqs_per_file:
            current_out_file.close()
            seq_counter = 0

    if not current_out_file.closed:
        current_out_file.close()

    return out_files
