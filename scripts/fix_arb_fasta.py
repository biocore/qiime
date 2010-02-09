#!/usr/bin/env python
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Pre-release"

from cogent.parse.fasta import MinimalFastaParser
from qiime.util import parse_command_line_parameters
from sys import argv
from optparse import make_option

script_description = """Fix arb fasta format by repairing incorrect line break chars, and stripping spaces and replacing . with - chars."""

script_usage = """Fix input arb fasta format file arb.fasta and print result to stdout:

fix_arb_fasta.py -i arb.fasta
"""

required_options = [\
    make_option('-i', '--input_file', help='input file in Arb fasta format'),
]

optional_options = []

def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)

    for label, seq in MinimalFastaParser(open(opts.input_file, 'U')):
        print '>%s\n%s' % (label, seq.replace(' ','').replace('.','-'))

if __name__ == '__main__':
    main()
