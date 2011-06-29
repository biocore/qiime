#!/usr/bin/env python
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Release"

from cogent.parse.fasta import MinimalFastaParser
from qiime.util import parse_command_line_parameters, get_options_lookup
from sys import argv
from qiime.util import make_option

options_lookup = get_options_lookup()

#fix_arb_fasta.py
script_info={}
script_info['brief_description']="""Reformat ARB FASTA files"""
script_info['script_description']="""This script fixes ARB FASTA formatting by repairing incorrect line break chararcters, stripping spaces and replacing "." with "-" characters."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""Fix the input ARB FASTA format file arb.fasta and print the result to stdout:""","""fix_arb_fasta.py -i arb.fasta"""))
script_info['output_description']="""The reformatted sequences are written to stdout."""
script_info['required_options']=[options_lookup['input_fasta']]
script_info['optional_options']=[]
script_info['version'] = __version__
def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    for label, seq in MinimalFastaParser(open(opts.input_fasta_fp, 'U')):
        print '>%s\n%s' % (label, seq.replace(' ','').replace('.','-'))

if __name__ == '__main__':
    main()
