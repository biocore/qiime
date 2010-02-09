#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Rob Knight"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Pre-release"
 

from qiime.util import parse_command_line_parameters
from optparse import make_option
from qiime.split_libraries import qual_scores
from qiime.make_fastq import make_fastq_single, make_fastq_multi

script_description = """Generates fastq file for ERA submission from paired fasta and qual files.

The fastq format for each record is as follows:
@seq_id [and optional description]
seq as bases
+[optionally with repeat of seq_id and repeat line]
qual scores as string of chr(33+qual)

The ERA currently requires a separate fastq file for each library, split by
library id. This code takes the output from split_libraries.py and the
corresponding qual files, pulls the qual info by id, and writes everything
either to one file or to per-library files."""

script_usage = """usage: Take input fasta file input_fasta_filepath and qual file input_qual_filepath: make separate file for each library (with the -s option: assumes that the fasta file is the output of split_libraries.py or similar script).

make_fasta.py -f input_fasta_filepath -q input_qual_filepath -s"""

required_options = [\
    make_option('-f', '--fasta', dest='fasta_fp',
        help='name of fasta output of split_libraries.py'),
    make_option('-q', '--qual', dest='qual_fps',
        help='names of qual files, comma-delimited'),
]

optional_options = [\
    make_option('-o','--result_fp',action='store',default=None,\
          type='string',dest='result_fp',help='Path to store '+\
          'results [default: <input_sequences_filename>.fastq]'),
    make_option('-s', '--split', dest='split', action='store_true',
        help='make separate file for each library [default:%default]',
        default=False)
]


def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)

    in_fasta = open(opts.fasta_fp, 'U')
    quals = qual_scores([open(f, 'U') for f in opts.qual_fps.split(',')])
    if not opts.result_fp:
        opts.result_fp = opts.fasta_fp + '.fastq'

    if opts.split:
        make_fastq_multi(in_fasta, quals, opts.result_fp)
    else:
        make_fastq_single(in_fasta, quals, opts.result_fp)


if __name__ == "__main__":
    main()
