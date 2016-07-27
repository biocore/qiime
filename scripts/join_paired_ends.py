#!/usr/bin/env python

__author__ = "Mike Robeson"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Mike Robeson"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Mike Robeson"
__email__ = "robesonms@ornl.gov"

from qiime.join_paired_ends import (join_method_names,
                                    join_method_constructors,
                                    write_synced_barcodes_fastq)
from qiime.util import (parse_command_line_parameters, get_options_lookup,
                        make_option, load_qiime_config, create_dir)
import os
import gzip


options_lookup = get_options_lookup()
qiime_config = load_qiime_config()

script_info = {}
script_info['brief_description'] = """Joins paired-end Illumina reads."""
script_info['script_description'] = """This script takes forward and reverse Illumina reads and joins them using the method chosen. Will optionally create an updated index reads file containing index reads for the surviving joined paired end reads. If the option to write an updated index file is chosen, be sure that the order and header format of the index reads is the same as the order and header format of reads in the files that will be joined (this is the default for reads generated on the Illumina instruments).

Currently, there are two methods that can be selected by the user to join paired-end data:

1. fastq-join - Erik Aronesty, 2011. ea-utils : "Command-line tools for processing biological sequencing data" (http://code.google.com/p/ea-utils)

2. SeqPrep - (https://github.com/jstjohn/SeqPrep)
"""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Join paired-ends with \'fastq-join\':""",
     """This is the default method to join paired-end Illumina data:""",
     """ %prog -f $PWD/forward_reads.fastq -r $PWD/reverse_reads.fastq -o $PWD/fastq-join_joined"""))
script_info['script_usage'].append(
    ("""Join paired-ends with \'SeqPrep\':""",
     """Produces similar output to the \'fastq-join\' but returns data in gzipped format.""",
     """ %prog -m SeqPrep -f $PWD/forward_reads.fastq -r $PWD/reverse_reads.fastq -o $PWD/SeqPrep_joined"""))
script_info['script_usage'].append(
    ("""Update the index / barcode reads file to match the surviving joined pairs.""",
     """This is required if you will be using split_libraries_fastq.py.""",
     """ %prog -f $PWD/forward_reads.fastq -r $PWD/reverse_reads.fastq -b $PWD/barcodes.fastq -o $PWD/fastq-join_joined"""))
script_info['output_description'] = """All paired-end joining software will return a joined / merged / assembled paired-end fastq file. Depending on the method chosen, additional files may be written to the user-specified output directory.


1. fastq-join will output fastq-formatted files as:

   - \"\*.join\": assembled / joined reads output
   - \"\*.un1\": unassembled / unjoined reads1 output
   - \"\*.un2\": unassembled / unjoined reads2 output

2. SeqPrep will output fastq-formatted gzipped files as:

   - \"\*_assembled.gz\": unassembled / unjoined reads1 output
   - \"\*_unassembled_R1.gz\": unassembled / unjoined reads1 output
   - \"\*_unassembled_R2.gz\": unassembled / unjoined reads2 output

3. If a barcode / index file is provided via the \'-b\' option, an updated
   barcodes file will be output as:

   - \"..._barcodes.fastq\": This barcode / index file must be used in
     conjunction with the joined
     paired-ends file as input to split_libraries_fastq.py. Except for
     missing reads that may result from failed merging of paired-ends, the
     index-reads and joined-reads must be in the same order.

"""
script_info['required_options'] = [
    make_option('-f', '--forward_reads_fp', type="existing_filepath",
                help='Path to input forward reads in FASTQ format.'),
    make_option('-r', '--reverse_reads_fp', type="existing_filepath",
                help='Path to input reverse reads in FASTQ format.'),
    make_option('-o', '--output_dir', type='new_dirpath',
                help='Directory to store result files')]
script_info['optional_options'] = [
    make_option('-m', '--pe_join_method', type='choice',
                choices=list(join_method_names.keys()),
                help='Method to use for joining paired-ends. Valid choices' +
                      ' are: ' + ', '.join(join_method_names.keys()) +
                      ' [default: %default]', default='fastq-join'),
    make_option('-b', '--index_reads_fp', type='existing_filepath',
                help='Path to the barcode / index reads in FASTQ format.'
                ' Will be filtered based on surviving joined pairs.'),
    make_option('-j', '--min_overlap', type='int',
                help='Applies to both fastq-join and SeqPrep methods.' +
                      ' Minimum allowed overlap in base-pairs required to join pairs.' +
                      ' If not set, progam defaults will be used. For example, for fastq-join (6 bp) will be used.'
                      ' Must be an integer. [default: %default]', default=None),
    make_option('-p', '--perc_max_diff', type='int',
                help='Only applies to fastq-join method, otherwise ignored. ' +
                     'Maximum allowed % differences within region of overlap.' +
                      ' If not set, progam defaults will be used. For example, for fastq-join (8%) will be used.' +
                      ' Must be an integer between 1-100 [default: %default]',
                default=None),
    make_option('-y', '--max_ascii_score',
                help='Only applies to SeqPrep method, otherwise ignored.' +
                      ' Maximum quality score / ascii code allowed to appear within' +
                      ' joined pairs output. For more information, please see:' +
                      ' http://en.wikipedia.org/wiki/FASTQ_format.'
                      ' [default: %default]', default='J'),
    make_option('-n', '--min_frac_match', type='float',
                help='Only applies to SeqPrep method, otherwise ignored.' +
                      ' Minimum allowed fraction of matching bases required' +
                      ' to join reads. Must be a float between 0-1.' +
                      ' If not set, progam defaults will be used.' +
                      ' [default: %default]', default=None),
    make_option('-g', '--max_good_mismatch', type='float',
                help='Only applies to SeqPrep method, otherwise ignored.' +
                      ' Maximum mis-matched high quality bases allowed' +
                      ' to join reads. Must be a float between 0-1.' +
                      ' If not set, progam defaults will be used.' +
                      ' [default: %default]', default=None),
    make_option('-6', '--phred_64',
                help='Only applies to SeqPrep method, otherwise ignored.' +
                      ' Set if input reads are in phred+64 format. Output will '
                      'always be phred+33. [default: %default]',
                default=False)]

script_info['version'] = __version__


def main():
    # parse command line parameters
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # Create local copy of options
    forward_reads_fp = opts.forward_reads_fp
    reverse_reads_fp = opts.reverse_reads_fp
    pe_join_method = opts.pe_join_method
    output_dir = opts.output_dir
    # fastq-join only options:
    perc_max_diff = opts.perc_max_diff
    # SeqPrep only options:
    max_ascii_score = opts.max_ascii_score
    min_frac_match = opts.min_frac_match
    max_good_mismatch = opts.max_good_mismatch
    phred_64 = opts.phred_64
    # both fastq-join & SeqPrep options
    min_overlap = opts.min_overlap

    create_dir(output_dir, fail_on_exist=False)

    # send parameters to appropriate join method
    # currently only two join methods exist:
    # 'fastq-join' and 'SeqPrep'
    if pe_join_method == "fastq-join":
        join_func = join_method_names["fastq-join"]
        paths = join_func(forward_reads_fp,
                          reverse_reads_fp,
                          perc_max_diff=perc_max_diff,
                          min_overlap=min_overlap,
                          working_dir=output_dir)

    if pe_join_method == "SeqPrep":
        join_func = join_method_names["SeqPrep"]
        paths = join_func(forward_reads_fp,
                          reverse_reads_fp,
                          max_overlap_ascii_q_score=max_ascii_score,
                          min_overlap=min_overlap,
                          max_mismatch_good_frac=max_good_mismatch,
                          min_frac_matching=min_frac_match,
                          phred_64=phred_64,
                          working_dir=output_dir)

    # If index / barcode file is supplied, filter unused barcode reads
    # and write them to a new file. Name based on joined-pairs / assembled
    # outfile
    if opts.index_reads_fp:
        index_reads = opts.index_reads_fp
        assembly_fp = paths['Assembled']  # grab joined-pairs output path
        write_synced_barcodes_fastq(assembly_fp, index_reads)


if __name__ == "__main__":
    main()
