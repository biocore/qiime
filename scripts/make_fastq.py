#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Jens Reeder"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Kyle Bittinger"
__email__ = "kylebittinger@gmail.com"


from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from qiime.parse import parse_qual_scores
from qiime.make_fastq import make_fastq_single, make_fastq_multi

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Make FASTQ file for ERA submission from \
paired FASTA and QUAL files"
script_info['script_description'] = "The ERA currently requires a separate \
FASTQ file for each library, split by library id. This code takes the output \
from split_libraries.py and the corresponding QUAL files and produces \
ERA-compatible FASTQ files."
script_info['script_usage'] = []
script_info['script_usage'].append(("Example:", "Take input FASTA file \
input_fasta_filepath and QUAL file input_qual_filepath: make separate file \
for each library (with the -s option: assumes that the FASTA file is the \
output of split_libraries.py or similar script):", "%prog -f $PWD/seqs.fna -q \
$PWD/Fasting_Example.qual -s"))
script_info['output_description'] = """Matches QUAL info to FASTA entries by id,\
 and writes FASTQ output to one file or to per-library files.

The FASTQ format for each record is as follows:

@seq_id [and optional description]
seq as bases
+ [and optionally with repeat of seq_id and repeat line]
qual scores as string of chr(33+qual)
"""

script_info['required_options'] = [options_lookup['input_fasta'],
                                   make_option('-q', '--qual', dest='qual_fps',
                                               type='existing_filepaths',
                                               help='names of QUAL files, comma-delimited'),
                                   ]
script_info['optional_options'] = [
    make_option('-o', '--result_fp', default=None,
                type='new_filepath', help='Path to store results '
                '[default: <input_sequences_filename>.fastq]'),
    make_option('-s', '--split', dest='split', action='store_true',
                help='make separate file for each library [default:%default]',
                default=False)
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    in_fasta = open(opts.input_fasta_fp, 'U')
    quals = parse_qual_scores([open(f, 'U') for f in opts.qual_fps])
    if not opts.result_fp:
        opts.result_fp = opts.input_fasta_fp + '.fastq'

    if opts.split:
        make_fastq_multi(in_fasta, quals, opts.result_fp)
    else:
        make_fastq_single(in_fasta, quals, opts.result_fp)


if __name__ == "__main__":
    main()
