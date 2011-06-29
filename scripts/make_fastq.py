#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Jens Reeder"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"
__status__ = "Release"
 

from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from qiime.parse import parse_qual_scores
from qiime.make_fastq import make_fastq_single, make_fastq_multi

options_lookup = get_options_lookup()

#make_fastq.py - Could use an expanded output description.
script_info={}
script_info['brief_description']="""Make fastq file for ERA submission from paired fasta and qual files"""
script_info['script_description']="""The ERA currently requires a separate fastq file for each library, split by library id. This code takes the output from split_libraries.py and the corresponding qual files, pulls the qual info by id, and writes everything either to one file or to per-library files.

The fastq format for each record is as follows:

- @seq_id [and optional description]
- seq as bases + [optionally with repeat of seq_id and repeat line]
- qual scores as string of chr(33+qual)
"""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""Take input fasta file input_fasta_filepath and qual file input_qual_filepath: make separate file for each library (with the -s option: assumes that the fasta file is the output of split_libraries.py or similar script):""","""make_fasta.py -f input_fasta_filepath -q input_qual_filepath -s"""))
script_info['output_description']="""This script creates separate fastq files for each library."""
script_info['required_options']=[options_lookup['input_fasta'],\
    make_option('-q', '--qual', dest='qual_fps',
        help='names of qual files, comma-delimited'),
]
script_info['optional_options']=[\
    make_option('-o','--result_fp',action='store',default=None,\
          type='string',dest='result_fp',help='Path to store '+\
          'results [default: <input_sequences_filename>.fastq]'),
    make_option('-s', '--split', dest='split', action='store_true',
        help='make separate file for each library [default:%default]',
        default=False)
]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    in_fasta = open(opts.input_fasta_fp, 'U')
    quals = parse_qual_scores([open(f, 'U') for f in opts.qual_fps.split(',')])
    if not opts.result_fp:
        opts.result_fp = opts.input_fasta_fp + '.fastq'

    if opts.split:
        make_fastq_multi(in_fasta, quals, opts.result_fp)
    else:
        make_fastq_single(in_fasta, quals, opts.result_fp)


if __name__ == "__main__":
    main()
