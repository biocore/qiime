#!/usr/bin/env python
from __future__ import division

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@gmail.com"

from os.path import split, splitext

from qiime.util import parse_command_line_parameters, get_options_lookup,\
    make_option, subsample_fasta

options_lookup = get_options_lookup()

script_info = {}
script_info[
    'brief_description'] = """Randomly subsample sequences from a given fasta file"""
script_info[
    'script_description'] = """Subsample the seqs.fna file, randomly select 5% of the sequences:"""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Example:""",
     """Subsample seqs.fasta to approximately 5%""",
     """%prog -i $PWD/seqs.fna -p 0.05 -o $PWD/subsampled_seqs.fna"""))
script_info['output_description'] = """"""
script_info['required_options'] = [
    options_lookup['fasta_as_primary_input'],
    make_option('-p', '--percent_subsample', action='store', type='float',
                help='Specify the percentage (as a fraction between 0 and 1) '
                'of sequences to subsample')
]
script_info['optional_options'] = [
    options_lookup['output_fp']
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    verbose = opts.verbose

    input_fasta_fp = opts.input_fasta_fp
    output_fp = opts.output_fp
    percent_subsample = opts.percent_subsample

    if percent_subsample > 1 or percent_subsample <= 0:
        raise ValueError(('percent_subsample must be in range of 0-1'))

    if not output_fp:
        input_file_basename, input_file_ext = \
            splitext(split(input_fasta_fp)[1])
        output_fp = '%s_subsample_%3.2f%s' % (input_file_basename,
                                              percent_subsample, input_file_ext)

    subsample_fasta(input_fasta_fp, output_fp, percent_subsample)


if __name__ == "__main__":
    main()
