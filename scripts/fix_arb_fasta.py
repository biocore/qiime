#!/usr/bin/env python
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Daniel McDonald", "Yoshiki Vazquez Baeza"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

from skbio.parse.sequences import parse_fasta
from qiime.util import parse_command_line_parameters, get_options_lookup
from sys import stdout
from qiime.util import make_option

options_lookup = get_options_lookup()

# fix_arb_fasta.py
script_info = {}
script_info['brief_description'] = """Reformat ARB FASTA files"""
script_info[
    'script_description'] = """This script fixes ARB FASTA formatting by repairing incorrect line break chararcters, stripping spaces and replacing "." with "-" characters."""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Example:""",
     """Fix the input ARB FASTA format file arb.fasta and print the result to stdout:""",
     """%prog -f arb.fasta"""))
script_info['script_usage'].append(
    ("""Example saving to an output file:""",
     """Fix the input ARB FASTA format file arb.fasta and print the result to fixed.fasta:""",
     """%prog -f arb.fasta -o fixed.fasta"""))
script_info[
    'output_description'] = """The reformatted sequences are written to stdout or to the file path provided with -o."""
script_info['required_options'] = [options_lookup['input_fasta']]
script_info['optional_options'] = [
    make_option('-o', '--output_fp', type="new_filepath",
                help="path where output will be written [default: print to screen]",
                default=None)
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    output_fp = opts.output_fp
    # if no output file path is provided then print to stdout
    if output_fp:
        fd = open(output_fp, 'w')
    else:
        fd = stdout

    for label, seq in parse_fasta(open(opts.input_fasta_fp, 'U')):
        print >> fd, '>%s\n%s' % (
            label, seq.replace(' ', '').replace('.', '-'))

    fd.close()

if __name__ == '__main__':
    main()
