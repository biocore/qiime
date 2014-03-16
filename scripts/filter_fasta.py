#!/usr/bin/env python
# File created on 18 May 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Jens Reeder"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


from qiime.util import make_option
from skbio.parse.sequences import fasta_parse
from cogent.parse.fastq import fastq_parse
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.parse import fields_to_dict
from qiime.filter import (filter_fasta, filter_fastq,
                          get_seqs_to_keep_lookup_from_seq_id_file,
                          get_seqs_to_keep_lookup_from_fasta_file,
                          sample_ids_from_metadata_description,
                          get_seqs_to_keep_lookup_from_biom)

options_lookup = get_options_lookup()

script_info = {}
script_info[
    'brief_description'] = "This script can be applied to remove sequences from a fasta or fastq file based on input criteria."
script_info['script_description'] = ""
script_info['script_usage'] = []

script_info[
    'script_usage'].append(("OTU map-based filtering", "Keep all sequences that show up in an OTU map.",
                            "%prog -f inseqs.fasta -o otu_map_filtered_seqs.fasta -m otu_map.txt"))

script_info[
    'script_usage'].append(("Chimeric sequence filtering", "Discard all sequences that show up in chimera checking output. NOTE: It is very important to pass -n here as this tells the script to negate the request, or discard all sequences that are listed via -s. This is necessary to remove the identified chimeras from inseqs.fasta.",
                            "%prog -f inseqs.fasta -o non_chimeric_seqs.fasta -s chimeric_seqs.txt -n"))

script_info[
    'script_usage'].append(("Sequence list filtering", "Keep all sequences from as fasta file that are listed in a text file.",
                            "%prog -f inseqs.fasta -o list_filtered_seqs.fasta -s seqs_to_keep.txt"))

script_info[
    'script_usage'].append(("biom-based filtering", "Keep all sequences that are listed as observations in a biom file.",
                            "%prog -f inseqs.fastq -o biom_filtered_seqs.fastq -b otu_table.biom"))

script_info[
    'script_usage'].append(("fastq filtering", "Keep all sequences from a fastq file that are listed in a text file (note: file name must end with .fastq to support fastq filtering).",
                            "%prog -f inseqs.fastq -o list_filtered_seqs.fastq -s seqs_to_keep.txt"))

script_info[
    'script_usage'].append(("sample id list filtering", "Keep all sequences from a fasta file where the sample id portion of the sequence identifier is listed in a text file (sequence identifiers in fasta file must be in post-split libraries format: sampleID_seqID).",
                            "%prog -f sl_inseqs.fasta -o sample_id_list_filtered_seqs.fasta --sample_id_fp map.txt"))

script_info['output_description'] = ""
script_info['required_options'] = [
    options_lookup['input_fasta'],
    make_option(
        '-o',
        '--output_fasta_fp',
        type='new_filepath',
        help='the output fasta filepath')
]
script_info['optional_options'] = [
    make_option('-m', '--otu_map', type='existing_filepath',
                help='an OTU map where sequences ids are those which should be retained'),
    make_option('-s', '--seq_id_fp', type='existing_filepath',
                help='A list of sequence identifiers (or tab-delimited lines with'
                ' a seq identifier in the first field) which should be retained'),
    make_option('-b', '--biom_fp', type='existing_filepath',
                help='A biom file where otu identifiers should be retained'),
    make_option('-a', '--subject_fasta_fp', type='existing_filepath',
                help='A fasta file where the seq ids should be retained.'),
    make_option('-p', '--seq_id_prefix', type='string',
                help='keep seqs where seq_id starts with this prefix'),
    make_option('--sample_id_fp', type='existing_filepath',
                help='keep seqs where seq_id starts with a sample id listed in this file'),
    make_option('-n', '--negate', help='discard passed seq ids rather than'
                ' keep passed seq ids [default: %default]', default=False,
                action='store_true'),
    make_option('--mapping_fp', type='existing_filepath',
                help='mapping file path (for use with --valid_states) [default: %default]'),
    make_option('--valid_states', type='string',
                help='description of sample ids to retain (for use with --mapping_fp) [default: %default]')
]
script_info['version'] = __version__


def filter_fasta_fp(input_seqs_fp, output_seqs_fp, seqs_to_keep, negate=False):
    """Filter a fasta file to include only sequences listed in seqs_to_keep """
    input_seqs = fasta_parse(open(input_seqs_fp, 'U'))
    output_f = open(output_seqs_fp, 'w')
    return filter_fasta(input_seqs, output_f, seqs_to_keep, negate)


def filter_fastq_fp(input_seqs_fp, output_seqs_fp, seqs_to_keep, negate=False):
    """Filter a fastq file to include only sequences listed in seqs_to_keep """
    input_seqs = fastq_parse(open(input_seqs_fp, 'U'), strict=False)
    output_f = open(output_seqs_fp, 'w')
    return filter_fastq(input_seqs, output_f, seqs_to_keep, negate)


def get_seqs_to_keep_lookup_from_otu_map(seqs_to_keep_f):
    """Generate a lookup dictionary from an OTU map"""
    otu_map = fields_to_dict(seqs_to_keep_f)
    seqs_to_keep = []
    for seq_ids in otu_map.values():
        seqs_to_keep += seq_ids
    return {}.fromkeys(seqs_to_keep)


def get_seqs_to_keep_lookup_from_prefix(fasta_f, prefix):
    seqs_to_keep = [seq_id
                    for seq_id, seq in fasta_parse(fasta_f)
                    if seq_id.startswith(prefix)]
    return {}.fromkeys(seqs_to_keep)


def get_seqs_to_keep_lookup_from_sample_ids(fasta_f, sample_ids):
    sample_ids = set(sample_ids)
    seqs_to_keep = set()
    for seq_id, seq in fasta_parse(fasta_f):
        if seq_id.split('_')[0] in sample_ids:
            seqs_to_keep.add(seq_id)
    return {}.fromkeys(seqs_to_keep)


def get_seqs_to_keep_lookup_from_mapping_file(
        fasta_f, mapping_f, valid_states):
    sample_ids = {}.fromkeys(
        sample_ids_from_metadata_description(mapping_f, valid_states))
    seqs_to_keep = []
    for seq_id, seq in fasta_parse(fasta_f):
        if seq_id.split('_')[0] in sample_ids:
            seqs_to_keep.append(seq_id)
        else:
            continue
    return {}.fromkeys(seqs_to_keep)


def main():
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)

    negate = opts.negate
    error_msg = "Must pass exactly one of -a, -b, -s, -p, -m, or --valid_states and --mapping_fp."
    if 1 != sum(map(bool, [opts.otu_map,
                           opts.seq_id_fp,
                           opts.subject_fasta_fp,
                           opts.seq_id_prefix,
                           opts.biom_fp,
                           opts.sample_id_fp,
                           opts.mapping_fp and opts.valid_states])):
        option_parser.error(error_msg)

    if opts.otu_map:
        seqs_to_keep_lookup =\
            get_seqs_to_keep_lookup_from_otu_map(
                open(opts.otu_map, 'U'))
    elif opts.seq_id_fp:
        seqs_to_keep_lookup =\
            get_seqs_to_keep_lookup_from_seq_id_file(
                open(opts.seq_id_fp, 'U'))
    elif opts.subject_fasta_fp:
        seqs_to_keep_lookup =\
            get_seqs_to_keep_lookup_from_fasta_file(
                open(opts.subject_fasta_fp, 'U'))
    elif opts.seq_id_prefix:
        seqs_to_keep_lookup =\
            get_seqs_to_keep_lookup_from_prefix(
                open(opts.input_fasta_fp), opts.seq_id_prefix)
    elif opts.mapping_fp and opts.valid_states:
        seqs_to_keep_lookup =\
            get_seqs_to_keep_lookup_from_mapping_file(
                open(opts.input_fasta_fp, 'U'),
                open(opts.mapping_fp, 'U'),
                opts.valid_states)
    elif opts.biom_fp:
        seqs_to_keep_lookup = \
            get_seqs_to_keep_lookup_from_biom(open(opts.biom_fp, 'U'))
    elif opts.sample_id_fp:
        sample_ids = set([e.strip().split()[0]
                         for e in open(opts.sample_id_fp, 'U')])
        seqs_to_keep_lookup = get_seqs_to_keep_lookup_from_sample_ids(
            open(opts.input_fasta_fp), sample_ids)
    else:
        option_parser.error(error_msg)

    if opts.input_fasta_fp.endswith('.fastq'):
        filter_fp_f = filter_fastq_fp
    else:
        filter_fp_f = filter_fasta_fp

    filter_fp_f(opts.input_fasta_fp,
                opts.output_fasta_fp,
                seqs_to_keep_lookup,
                negate)

if __name__ == "__main__":
    main()
