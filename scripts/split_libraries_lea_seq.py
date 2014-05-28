#!/usr/bin/env python
from __future__ import division

__author__ = "Charudatta Navare"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Charudatta Navare", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Charudatta Navare"
__email__ = "charudatta.navare@gmail.com"

import tempfile
import os
from qiime.util import create_dir
from qcli import parse_command_line_parameters, make_option
from qiime.golay import decode_golay_12
from qiime.split_libraries_lea_seq import read_input_file


script_info = {}
script_info[
    'brief_description'] = "Implements Low-Error Amplicon Sequencing (LEA-Seq)"
script_info[
    'script_description'] = """Implements Low-Error Amplicon Sequencing (LEA-Seq) method, \
described in: Faith, Jeremiah J., et al. \
The long-term stability of the human gut microbiota.\
Science 341.6141 (2013).\
This method is based on redundant sequencing of a set of linear PCR\
template extensionsof 16S rRNA genes. The oligonucleotide primer\
that is used for PCR template extensions is labeled with a random barcode\
5' to the universal 16S rRNA primer sequence. This PCR pool is then\
amplified with exponential PCR, using primers that specifically\
amplify only the linear PCR molecules. An index primer is added to\
the amplicons along with a primer specific for each sample.
This exponential PCR pool is then sequenced redundantly (20x coverage).\
The resulting sequences are separated by sample, using the index sequence.\
The amplicon sequences within each sample are separated by the random\
barcodes. The large number of reads for each barcode helps to\
create an error-corrected consensus sequence for the\
initial template molecule.
"""
script_info['script_usage'] = []
script_info['script_usage'].append(("""Example:""",
                                    """Specify forward read and reverse read fasta files,
                                    use the metadata mapping file map.txt,
                                    and output the data to output_dir""",
                                    """%prog -i fwd_read.fq,rev_read.fq -m map.txt -o output_dir --barcode_type=7"""))

script_info['output_description'] = """The %prog generates:\
A fasta file called seqs.fna which contains\
error corrected consensus sequence for the template DNA\
"""
script_info['required_options'] = [
    make_option('-i', '--sequence_read_fps', type='existing_filepaths',
                help='the forward and reverse sequence read fastq files '
                '(comma-separated)'),
    make_option('-o', '--output_dir', type='new_dirpath',
                help='directory to store output files'),
    make_option('-m', '--mapping_fp', type='existing_filepath',
                help='metadata mapping file')
]
script_info['optional_options'] = [
    make_option('--barcode_type', type='string',
                help='the type of barcode used. This can be an integer, e.g. '
                '6 for length 6 barcodes, or golay_12 for golay error-'
                'correcting barcodes. Error correction will only be '
                'applied for golay_12 barcodes [default: %default]',
                default='golay_12'),
    make_option('--max_barcode_errors', type='float',
                help='the maximum allowable number of errors in the barcode '
                'if passing --barcode_type golay_12 [default: %default]',
                default=1.5),
    make_option('--min_consensus', type='float',
                help='threshold for consensus score'
                'the minimum score allowable at any position in sequence'
                'where the score is calulated as:'
                'occurence of base in consensus sequence/ total sequences'
                '[default: %default]',
                default=0.66),
    make_option('--max_cluster_ratio', type='float',
                help='threshold for cluster ratio'
                'the maximum allowable cluster ratio'
                'above which you need to find the consensus sequence'
                'for the given sequences'
                '[default: %default]',
                default=2.5),
    make_option('--min_difference_in_bcs', type='float',
                help='threshold for selecting unique barcodes.'
                'Barcodes that are more similar to each other'
                'than this value will be discarded'
                '[default: %default]',
                default=0.86)


]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    barcode_type = opts.barcode_type
    max_barcode_errors = opts.max_barcode_errors
    mapping_fp = opts.mapping_fp
    sequence_read_fps = opts.sequence_read_fps
    min_consensus = opts.min_consensus
    max_cluster_ratio = opts.max_cluster_ratio
    output_dir = opts.output_dir
    min_difference_in_bcs = opts.min_difference_in_bcs
    create_dir(output_dir)
    consensus_outfile = open(os.path.join(output_dir, "seqs.fna"), "w")
   
    if barcode_type == 'golay_12':
        barcode_correction_fn = decode_golay_12
        barcode_len = 12
    else:
        barcode_correction_fn = None

        try:
            barcode_len = int(barcode_type)
        except ValueError:
            option_parser.error("Invalid barcode type '%s'. The barcode type "
                                "must be either golay_12 or a positive "
                                "integer indicating the barcode length." %
                                barcode_type)

    if max_barcode_errors < 0:
        option_parser.error("--max_barcode_errors must be greater than or "
                            "equal to zero. You provided %.4f." %
                            max_barcode_errors)

    if barcode_len < 1:
        option_parser.error("Invalid barcode length: %d. Must be greater "
                            "than zero." % barcode_len)



    if len(sequence_read_fps) != 2:
        option_parser.error("You must provide exactly two sequence read "
                            "filepaths, the first for forward reads and "
                            "second for reverse reads. You specified %d "
                            "filepaths." % len(sequence_read_fps))


    consensus_seq_lookup = read_input_file(sequence_read_fps, mapping_fp,
                                           output_dir, barcode_type, barcode_correction_fn,
                                           max_barcode_errors, min_consensus,
                                           max_cluster_ratio, min_difference_in_bcs)

    for sample_id in consensus_seq_lookup:
        for random_bc in consensus_seq_lookup[sample_id]:
            consensus_seq = consensus_seq_lookup[sample_id][random_bc]
            consensus_outfile.write(">" + sample_id + random_bc
                                    + "\n" + consensus_seq + "\n")

if __name__ == "__main__":
    main()
