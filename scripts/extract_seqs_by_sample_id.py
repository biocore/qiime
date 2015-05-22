#!/usr/bin/env python
# File created on 08 Jan 2010.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


from skbio.parse.sequences import parse_fasta
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from qiime.util import extract_seqs_by_sample_id
from qiime.parse import parse_mapping_file
from qiime.filter import (parse_metadata_state_descriptions,
                          get_sample_ids)

options_lookup = get_options_lookup()

script_info = {}
script_info[
    'brief_description'] = """Extract sequences based on the SampleID"""
script_info[
    'script_description'] = """This script creates a fasta file which will contain only sequences that ARE associated with a set of sample IDs, OR all sequences that are NOT associated with a set of sample IDs (-n)"""
script_info['script_usage'] = []
script_info['script_usage'].append(("""Examples:""", """Create the file outseqs.fasta (-o), which will be a subset of inseqs.fasta (-i) containing only the sequences THAT ARE associated with sample ids S2, S3, \
S4 (-s). As always, sample IDs are case-sensitive:""", """%prog -i inseqs.fasta -o outseqs_by_sample.fasta -s S2,S3,S4"""))
script_info['script_usage'].append(
    ("""""",
     """Create the file outseqs.fasta (-o), which will be a subset of inseqs.fasta (-i) containing only the sequences THAT ARE NOT (-n) associated with sample ids S2, S3, S4 (-s). As always, sample IDs are case-sensitive:""",
     """%prog -i inseqs.fasta -o outseqs_by_sample_negated.fasta -s S2,S3,S4 -n"""))

script_info['script_usage'].append(
    ("""""",
     """Create the file outseqs.fasta (-o), which will be a subset of inseqs.fasta (-i) containing only the sequences THAT ARE associated with sample ids whose "Treatment" value is "Fast" in the mapping file:""",
     """%prog -i inseqs.fasta -o outseqs_by_mapping_field.fasta -m map.txt -s "Treatment:Fast" """))

script_info[
    'output_description'] = """The script produces a fasta file containing containing only the specified SampleIDs."""

script_info['required_options'] = [
    options_lookup['fasta_as_primary_input'],
    make_option('-o', '--output_fasta_fp', type="new_filepath",
                help='the output fasta file')
]

script_info['optional_options'] = [
    make_option('-n', '--negate', action='store_true', default=False,
                help='negate the sample ID list (i.e., output sample ' +
                'ids not passed via -s) [default: %default]'),
    make_option('-s', '--sample_ids', type='string',
                help="comma-separated sample_ids to include in output fasta file" +
                " (or exclude if --negate), or string describing mapping file states" +
                " defining sample ids (mapping_fp must be provided for the latter)"),
    options_lookup['mapping_fp']]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    negate = opts.negate
    sample_ids = opts.sample_ids
    mapping_fp = opts.mapping_fp
    input_fasta_fp = opts.input_fasta_fp
    output_fasta_fp = opts.output_fasta_fp

    if not mapping_fp:
        sample_ids = sample_ids.split(',')
    else:
        map_data, map_header, map_comments = parse_mapping_file(mapping_fp)
        sample_ids = get_sample_ids(
            map_data,
            map_header,
            parse_metadata_state_descriptions(sample_ids))
        if len(sample_ids) == 0:
            raise ValueError(
                "No samples match the search criteria: %s" %
                valid_states)

    if opts.verbose:
        # This is useful when using the --valid_states feature so you can
        # find out if a search query didn't work as you expected before a
        # lot of time is spent
        print "Extracting samples: %s" % ', '.join(sample_ids)

    try:
        seqs = parse_fasta(open(input_fasta_fp))
    except IOError:
        option_parser.error(
            'Cannot open %s. Does it exist? Do you have read access?' %
            input_fasta_fp)
        exit(1)

    try:
        output_fasta_f = open(output_fasta_fp, 'w')
    except IOError:
        option_parser.error(
            "Cannot open %s. Does path exist? Do you have write access?" %
            output_fasta_fp)
        exit(1)

    for r in extract_seqs_by_sample_id(seqs, sample_ids, negate):
        output_fasta_f.write('>%s\n%s\n' % r)
    output_fasta_f.close()

if __name__ == "__main__":
    main()
