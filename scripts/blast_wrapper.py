#!/usr/bin/env python
# File created on 20 Dec 2009.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from skbio.parse.sequences import parse_fasta
from qiime.util import qiime_blast_seqs

options_lookup = get_options_lookup()

# blast_wrapper.py
script_info = {}
script_info['brief_description'] = """Blast Interface"""
script_info[
    'script_description'] = """This script is a functionally-limited interface to the qiime.util.qiime_blast_seqs function, primarily useful for testing purposes. Once that function has been integrated into qiime as the primary blast interface it will move to PyCogent. An expanded version of this command line interface may replace the script functionality of cogent.app.blast at that point."""
script_info['script_usage'] = []
script_info['script_usage'].append(("""Example:""", """Blast all sequences in inseqs.fasta (-i) against a BLAST db constructed \
from refseqs.fasta (-r).""", """%prog -i $PWD/inseqs.fasta -r $PWD/refseqs.fasta"""))
script_info[
    'output_description'] = """This is a utility program, which returns BLAST results."""
script_info['required_options'] = [
    options_lookup['fasta_as_primary_input'],
    make_option('-r', '--refseqs_fp', type='string',
                help='path to blast database as a fasta file')
]
script_info['optional_options'] = [
    make_option('-n', '--num_seqs_per_blast_run', type='int', default='1000',
                help='number of sequences passed to each blast call ' +
                "- useful for very large sequence collections [default: %default]")
]

script_info['version'] = __version__


def main():
    option_parser, options, args = parse_command_line_parameters(**script_info)

    blast_results = qiime_blast_seqs(
        seqs=parse_fasta(open(options.input_fasta_fp)),
        refseqs_fp=options.refseqs_fp,
        seqs_per_blast_run=options.num_seqs_per_blast_run)

    for query_id, blast_result in blast_results.items():
        first_blast_result = blast_result[0][0]
        print '%s: %s %s %s' % (
            query_id,
            first_blast_result['SUBJECT ID'],
            first_blast_result['E-VALUE'],
            first_blast_result['% IDENTITY'])

if __name__ == "__main__":
    main()
