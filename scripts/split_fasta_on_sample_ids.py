#!/usr/bin/env python
# File created on 20 Oct 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"


from skbio.parse.sequences import fasta_parse
from qiime.util import (parse_command_line_parameters,
                        make_option,
                        split_fasta_on_sample_ids_to_files)

script_info = {}
script_info[
    'brief_description'] = "Split a single post-split_libraries.py fasta file into per-sample fasta files."
script_info[
    'script_description'] = "Split a single post-split_libraries.py fasta file into per-sample fasta files. This script requires that the sequences identitifers are in post-split_libraries.py format (i.e., SampleID_SeqID). A fasta file will be created for each unique SampleID."
script_info['script_usage'] = [(
    "",
    "Split seqs.fna into one fasta file per sample and store the resulting fasta files in 'out'",
    "%prog -i seqs.fna -o out/")]
script_info['script_usage_output_to_remove'] = ['$PWD/out/']
script_info[
    'output_description'] = "This script will produce an output directory with as many files as samples."
script_info['required_options'] = [
    make_option(
        '-i',
        '--input_fasta_fp',
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
]
script_info['version'] = __version__


def main():
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)

    split_fasta_on_sample_ids_to_files(
        MinimalFastaParser(open(opts.input_fasta_fp, 'U')),
        opts.output_dir,
        opts.buffer_size)


if __name__ == "__main__":
    main()
