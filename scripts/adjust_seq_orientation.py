#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"

from os.path import split, splitext

from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from qiime.adjust_seq_orientation import (rc_fasta_file,
    append_rc, null_seq_desc_mapper)

options_lookup = get_options_lookup()

script_info = {}
script_info[
    'brief_description'] = """Get the reverse complement of all sequences"""
script_info[
    'script_description'] = """Write the reverse complement of all seqs in seqs.fasta (-i) to seqs_rc.fasta (default, change output_fp with -o). Each sequence description line will have ' RC' appended to the end of it (default,
leave sequence description lines untouched by passing -r):"""
script_info['script_usage'] = []
script_info['script_usage'].append(("""Example:""",
                                    """Reverse complement all sequences in seqs.fna and write result to seqs_rc.fna""",
                                    """%prog -i seqs.fna"""))
script_info['output_description'] = """"""
script_info['required_options'] = [
    options_lookup['fasta_as_primary_input']
]
script_info['optional_options'] = [
    options_lookup['output_fp'],
    make_option('-r', '--retain_seq_id', action='store_true',
                help='leave seq description lines untouched' +
                ' [default: append " RC" to seq description lines]')
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    verbose = opts.verbose

    input_fasta_fp = opts.input_fasta_fp
    output_fp = opts.output_fp
    retain_seq_id = opts.retain_seq_id

    if retain_seq_id:
        seq_desc_mapper = null_seq_desc_mapper
    else:
        seq_desc_mapper = append_rc

    if not output_fp:
        input_file_basename, input_file_ext = \
            splitext(split(input_fasta_fp)[1])
        output_fp = '%s_rc%s' % (input_file_basename, input_file_ext)

    rc_fasta_file(input_fasta_fp, output_fp, seq_desc_mapper)


if __name__ == "__main__":
    main()
