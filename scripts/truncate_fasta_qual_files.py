#!/usr/bin/env python
# File created Sept 30, 2010
from __future__ import division

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Walters"
__email__ = "William.A.Walters@colorado.edu"

from qiime.util import make_option

from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.truncate_fasta_qual_files import truncate_fasta_qual
from qiime.util import create_dir


options_lookup = get_options_lookup()

script_info = {}
script_info[
    'brief_description'] = """Generates filtered fasta and quality score files by truncating at the specified base position."""
script_info['script_description'] = """This module is designed to remove regions of poor quality in
454 sequence data.  Drops in quality can be visualized with the
quality_scores_plot.py module.  The base position specified will
be used as an index to truncate the sequence and quality scores, and
all data at that base position and to the end of the sequence will be
removed in the output filtered files."""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Example:""",
     """Truncate the input fasta and quality files at base position 100, output to the filtered_seqs directory:""",
     """%prog -f seqs.fna -q seqs.qual -b 100 -o filtered_seqs/"""))
script_info[
    'output_description'] = """Filtered versions of the input fasta and qual file (based on input name with '_filtered' appended) will be generated in the output directory"""
script_info['required_options'] = [
    make_option('-f', '--fasta_fp',
                type='existing_filepath',
                help='Input fasta filepath to be truncated.'),

    make_option('-q', '--qual_fp',
                type='existing_filepath',
                help='Input quality scores filepath to be truncated.'),

    make_option('-b', '--base_pos', type='int',
                help='Nucleotide position to truncate the fasta and quality score ' +
                'files at.')
]

script_info['optional_options'] = [
    make_option('-o', '--output_dir',
                type='new_path',
                help='Output directory.  Will be created if does not exist.  ' +
                '[default: %default]', default="."),
]

script_info['version'] = __version__


def main():
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)

    fasta_fp = opts.fasta_fp
    qual_fp = opts.qual_fp
    output_dir = opts.output_dir
    base_pos = int(opts.base_pos)

    create_dir(output_dir)

    truncate_fasta_qual(fasta_fp, qual_fp, output_dir, base_pos)


if __name__ == "__main__":
    main()
