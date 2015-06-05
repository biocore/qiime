#!/usr/bin/env python
from __future__ import division

__author__ = "William Walters"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["William Walters"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@gmail.com"

from skbio.util import create_dir

from qiime.util import (parse_command_line_parameters, get_options_lookup,
                        make_option, qiime_open)
from qiime.extract_barcodes import extract_barcodes

options_lookup = get_options_lookup()

script_info = {}
script_info[
    'brief_description'] = "This script is designed to format fastq sequence and barcode data so they are compatible with split_libraries_fastq.py (see http://qiime.org/tutorials/processing_illumina_data.html)."
script_info[
    'script_description'] = "A variety of data formats are possible, depending upon how one utilized sequencing primers, designed primer constructs (e.g., partial barcodes on each end of the read), or processed the data (e.g., barcodes were put into the sequence labels rather than the reads). See various input examples below."
script_info['script_usage'] = []

script_info[
    'script_usage'].append(("Parse barcodes of 12 base pairs from the beginning of a single read. Will create an output fastq file of the barcodes and an output file of the reads supplied with the barcodes removed.", "",
                            "%prog -f inseqs.fastq -c barcode_single_end --bc1_len 12 -o processed_seqs"))

script_info[
    'script_usage'].append(("Parse barcodes of 12 base pairs from the beginning of a single read, reverse complement the barcodes before writing. Will create an output fastq file of the barcodes and an output file of the reads supplied with the barcodes removed", "",
                            "%prog -f inseqs.fastq -c barcode_single_end --bc1_len 12 -o processed_seqs --rev_comp_bc1"))

script_info[
    'script_usage'].append(("Parse barcodes of 6 base pairs from the beginning of paired reads. Will create an output fastq file of the barcodes and an output file of each of the reads supplied with the barcodes removed. The order of the barcodes written is determined by the order of the files passed (-f is written first, followed by -r)", "",
                            "%prog -f inseqs_R1.fastq -r inseqs_R2.fastq -c barcode_paired_end --bc1_len 6 --bc2_len 6 -o processed_seqs"))

script_info[
    'script_usage'].append(("Parse barcodes of 6 base pairs from the beginning of paired reads, attempt to orient reads based upon detection of forward and reverse primers in the mapping file. Will create an output fastq file of the barcodes and an output file of each of the reads supplied with the barcodes removed. The order of the barcodes written is determined by the order of the files passed (-f is written first, followed by -r)", "",
                            "%prog -f inseqs_R1.fastq -r inseqs_R2.fastq -c barcode_paired_end --map_fp mapping_data.txt --attempt_read_reorientation --bc1_len 6 --bc2_len 6 -o processed_seqs"))

script_info[
    'script_usage'].append(("Parse barcodes of 6 base pairs from the beginning, 8 base pairs at the end of a stitched read. Will create an output fastq file of the barcodes and an output fastq file of the stitched read supplied with the barcodes removed. The barcode at the beginning of the stitched read is written first, followed by the barcode at the end, unless reversed by the --switch_bc_order option is used", "",
                            "%prog -f inseqs_R1.fastq -c barcode_paired_stitched --bc1_len 6 --bc2_len 8 -o processed_seqs"))

script_info[
    'script_usage'].append(("Parse barcodes of 12 base pairs from labels of the input fastq file. Example label (note that the desired character preceding the barcode is '#'): @MCIC-SOLEXA_0051_FC:1:1:14637:1026#CGATGTGATTTC/1 This will create an output fastq file of the barcodes (no other sequence are written). A second file with barcodes in the label can be passed with -r, and if this is done, the combined barcodes from -f and -r will be written together", "",
                            "%prog -f inseqs_R1.fastq -c barcode_in_label --char_delineator '#' --bc1_len 12 -o processed_seqs"))


script_info[
    'output_description'] = "In the output directory, there will be fastq files (barcode file, and one or two reads files)"

script_info['required_options'] = [
    make_option('-f', '--fastq1', type='existing_filepath',
                help='input fastq filepath. This file is considered read 1.')]

script_info['optional_options'] = [

    make_option('-r', '--fastq2', type='existing_filepath', default=None,
                help='input fastq filepath. This file is considered read 2. '
                '[default: %default]'),

    make_option('-o', '--output_dir', default='.', type='new_dirpath',
                help='directory prefix for output files [default: %default]'),

    make_option('-c', '--input_type', default="barcode_single_end",
                type='choice', choices=["barcode_single_end", "barcode_paired_end",
                                        "barcode_paired_stitched", "barcode_in_label"], help='Specify '
                'the input type. barcode_single_end: Input is a single fastq file, '
                'that starts with the barcode sequence. barcode_paired_end: Input is '
                'a pair of fastq files (--fastq1 and --fastq2) that each begin with '
                'a barcode sequence. The barcode for fastq1 will be written first, '
                'followed by the barcode from fastq2. barcode_paired_stitched: '
                'Input is a single fastq file that has barcodes at the beginning and '
                'end. The barcode from the beginning of the read will be written first '
                'followed by the barcode from the end of the read, unless the order '
                'is switched with --switch_bc_order. barcode_in_label: Input is a '
                'one (--fastq1) or two (--fastq2) fastq files with the barcode written '
                'in the labels. [default: %default]'),

    make_option('-l', '--bc1_len', default=6, type='int', action='store',
                help='Specify the length, in base pairs, of barcode 1. This applies '
                'to the --fastq1 file and all options specified by --input_type '
                '[default: %default]'),

    make_option('-L', '--bc2_len', default=6, type='int', action='store',
                help='Specify the length, in base pairs, of barcode 2. This applies '
                'to the --fastq2 file and options "barcode_paired_end", '
                '"barcode_paired_stitched", and "barcode_in_label" for the --input_type'
                ' [default: %default]'),

    make_option('--rev_comp_bc1', default=False, action='store_true',
                help='Reverse complement barcode 1 before writing'
                ' [default: %default]'),

    make_option('--rev_comp_bc2', default=False, action='store_true',
                help='Reverse complement barcode 2 before writing'
                ' [default: %default]'),

    make_option(
        '-s', '--char_delineator', default=":", action='store', type='str',
        help='Character in fastq label that should immediately precede the '
        'barcode sequence. The length of the barcode is specified by the '
        '--bc1_len (and optionally --bc2_len if paired end files are used) '
        'parameter. [default: %default]'),

    make_option('--switch_bc_order', default=False, action='store_true',
                help='Reverse barcode order written when using the '
                '-c barcode_paired_stitched option. [default: %default]'),

    make_option('-m', '--mapping_fp', type='existing_filepath', default=None,
                help='Filepath of mapping file. NOTE: Must contain a header'
                ' line indicating SampleID in the first column and'
                ' BarcodeSequence in the second, LinkerPrimerSequence in the third '
                'and a ReversePrimer column before the final Description column. '
                'Needed for --attempt_read_orientation option. [default: %default]'),

    make_option('-a', '--attempt_read_reorientation', default=False,
                action='store_true', help='Will attempt to search for the forward and'
                ' reverse primer in the read and adjust the sequence orientation to '
                'match the orientation of the forward primer. An exact match for the  '
                'forward and reverse complemented versions of the primers are tested '
                'for, and sequences are reverse complemented, if necessary, before '
                'writing. Sequences without an exact match are written to a separate '
                'output fastq file, labeled as _no_primer_match.fastq. '
                '[default: %default]'),
    make_option('-d', '--disable_header_match', default=False,
                action='store_true', help='Enable this option to suppress header '
                'matching between input fastq files.'
                '[default: %default]')

]
script_info['version'] = __version__


def main():
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)

    if opts.attempt_read_reorientation:
        if not opts.mapping_fp:
            option_parser.error("To use --attempt_read_reorientation, one must "
                                "supply a mapping file that contains both LinkerPrimerSequence "
                                "and ReversePrimer columns.")
    if opts.input_type == "barcode_paired_end":
        if not opts.fastq2:
            option_parser.error("To use input_type of barcode_paired_end, "
                                "a second fastq file must be specified with --fastq2")

    if not opts.fastq2:
        disable_header_match = True
    else:
        disable_header_match = opts.disable_header_match

    fastq1 = qiime_open(opts.fastq1)
    if opts.fastq2:
        fastq2 = qiime_open(opts.fastq2)
    else:
        fastq2 = None
    create_dir(opts.output_dir)
    if opts.mapping_fp:
        map_fp = qiime_open(opts.mapping_fp)
    else:
        map_fp = None

    extract_barcodes(fastq1, fastq2, opts.output_dir, opts.input_type,
                     opts.bc1_len, opts.bc2_len, opts.rev_comp_bc1, opts.rev_comp_bc2,
                     opts.char_delineator, opts.switch_bc_order, map_fp,
                     opts.attempt_read_reorientation, disable_header_match)


if __name__ == "__main__":
    main()
