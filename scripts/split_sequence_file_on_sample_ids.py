#!/usr/bin/env python
# File created on 20 Oct 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


from skbio.parse.sequences import parse_fasta
from qiime.util import (parse_command_line_parameters,
                        make_option,
                        split_sequence_file_on_sample_ids_to_files)

script_info = {}
script_info[
    'brief_description'] = "Split a single post-split_libraries.py fasta (or post-split_libraries_fastq.py fastq) file into per-sample files."
script_info[
    'script_description'] = "Split a single post-split_libraries.py fasta (or post-split_libraries_fastq.py fastq) file into per-sample fasta files. This script requires that the sequences identitifers are in post-split_libraries.py format (i.e., SampleID_SeqID). A file will be created for each unique SampleID."
script_info['script_usage'] = [(
    "",
    "Split seqs.fna into one fasta file per sample and store the resulting fasta files in 'out'",
    "%prog -i seqs.fna -o out/"),
    ("",
    "Split seqs.fastq into one fastq file per sample and store the resulting fastq files in 'out_fastq'",
    "%prog -i seqs.fastq --file_type fastq -o out_fastq/")]
script_info['script_usage_output_to_remove'] = ['$PWD/out/']
script_info[
    'output_description'] = "This script will produce an output directory with as many files as samples."
script_info['required_options'] = [
    make_option(
        '-i',
        '--input_seqs_fp',
        type="existing_filepath",
        help='the input fasta file to split'),
    make_option(
        '-o',
        '--output_dir',
        type="new_dirpath",
        help='the output directory [default: %default]'),
]
script_info['optional_options'] = [
    make_option('--buffer_size', type="int", default=500,
                help="the number of sequences to read into memory before writing to file (you usually won't need to change this) [default: %default]"),
    make_option('--file_type', type=str, default='fasta',
                help="Type of file. Either fasta or fastq")
]
script_info['version'] = __version__


def main():
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)

    with open(opts.input_seqs_fp, 'U') as input_seqs_f:
        split_sequence_file_on_sample_ids_to_files(
            input_seqs_f,
            opts.file_type,
            opts.output_dir,
            opts.buffer_size)


if __name__ == "__main__":
    main()
