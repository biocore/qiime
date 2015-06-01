#!/usr/bin/env python
from __future__ import division

__author__ = "Charudatta Navare"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Charudatta Navare", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Charudatta Navare"
__email__ = "charudatta.navare@gmail.com"


from os import path

from qiime.util import create_dir
from qiime.golay import decode_golay_12
from qiime.split_libraries_lea_seq import get_LEA_seq_consensus_seqs
from qcli import parse_command_line_parameters, make_option


script_info = {}
script_info['brief_description'] = "Demultiplexes Low-Error Amplicon Sequencing (LEA-Seq) data"
script_info['script_description'] = """
Implements Low-Error Amplicon Sequencing (LEA-Seq) method, described in:

Faith, Jeremiah J., et al.
The long-term stability of the human gut microbiota.Science 341.6141 (2013).

This method is based on redundant sequencing of a set of linear PCR template
extensions of 16S rRNA genes. The oligonucleotide primer that is used for
PCR template extensions is labeled with a random barcode
5' to the universal 16S rRNA primer sequence. This PCR pool is then
amplified with exponential PCR, using primers that specifically
amplify only the linear PCR molecules. An index primer is added to
the amplicons along with a primer specific for each sample.
This exponential PCR pool is then sequenced redundantly (20x coverage).
The resulting sequences are separated by sample, using the index sequence.
The amplicon sequences within each sample are separated by the random
barcodes. The large number of reads for each barcode helps to create an
error-corrected consensus sequence for the initial template molecule.
"""
script_info['script_usage'] = []
script_info['script_usage'].append((
    """General Example: Specify forward read and reverse read \
fasta files, use the metadata mapping file map.txt,\
and output the data to output_dir""",
    """output_dir""",
    """%prog -i fwd_read.fq,rev_read.fq -m map.txt -o output --b 7"""
))
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
    make_option('-b', '--barcode_type', type='string',
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
                help='threshold for consensus score: '
                'the minimum score allowable at any position in sequence. '
                'where the score is calulated as: '
                'occurence of base in consensus sequence/ total sequences'
                '[default: %default]',
                default=6.6),
    make_option('--max_cluster_ratio', type='float',
                help='threshold for cluster ratio: '
                'the maximum allowable cluster ratio '
                'above which you need to find the consensus sequence '
                'for the given sequences.'
                '[default: %default]',
                default=2.5),
    make_option('--min_difference_in_bcs', type='float',
                help='threshold for selecting unique barcodes: '
                'Barcodes that are more similar to each other '
                'than this value will be discarded.'
                '[default: %default]',
                default=0.86),
    make_option('--fwd_length', type='int',
                help='removes phasing from forward read'
                'by truncating it to standard length for the region'
                '[default: %default]',
                default=64),
    make_option('--rev_length', type='int',
                help='removes phasing from reverse read'
                'by truncating it to standard length for the region'
                '[default: %default]',
                default=77),
    make_option('--min_difference_in_clusters', type='float',
                help='the percent identity threshold while using '
                'uclust to cluster sequence reads, which is helpful'
                'in measuring quality of sequencing.'
                '[default: %default]',
                default=0.98),
    make_option('--min_reads_per_random_bc', type='int',
                help='minimum number of reads per random'
                'barcode, attempts to remove random barcodes'
                ' that are sequencing errors of true barcodes'
                'might be useful in saving memory and time'
                '[default: %default]',
                default=1),
    make_option('--header_barcode_column', type='string',
                help='header of barcode column'
                '[default: %default]',
                default='BarcodeSequence'),
    make_option('--reverse_primer_column', type='string',
                help='header of reverse primer column'
                '[default: %default]',
                default='ReversePrimer'),
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
    fwd_length = opts.fwd_length
    rev_length = opts.rev_length
    min_reads_per_random_bc = opts.min_reads_per_random_bc
    min_diff_in_clusters = opts.min_difference_in_clusters
    barcode_column = opts.header_barcode_column
    reverse_primer_column = opts.reverse_primer_column
    create_dir(output_dir)
    fwd_consensus_outfile = open(path.join(output_dir, "fwd.fna"), "w")
    rev_consensus_outfile = open(path.join(output_dir, "rev.fna"), "w")
    log_file = open(path.join(output_dir, "log.txt"), "w")

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

    if min_diff_in_clusters < 0 or min_diff_in_clusters > 1:
        option_parser.error("--min_difference_in_clusters must be "
                            "between 0 to 1. You provided %.4f." %
                            min_diff_in_clusters)

    if min_difference_in_bcs < 0 or min_difference_in_bcs > 1:
        option_parser.error("--min_difference_in_bcs must be between 0 to 1."
                            " You provided %.4f." %
                            min_difference_in_bcs)

    if barcode_len < 1:
        option_parser.error("Invalid barcode length: %d. Must be greater "
                            "than zero." % barcode_len)

    if len(sequence_read_fps) != 2:
        option_parser.error("You must provide exactly two sequence read "
                            "filepaths, the first for forward reads and "
                            "second for reverse reads. You specified %d "
                            "filepaths." % len(sequence_read_fps))

    fwd_read_f = open(sequence_read_fps[0], 'U')
    rev_read_f = open(sequence_read_fps[1], 'U')

    map_f = open(mapping_fp, 'U')

    (consensus_seq_lookup,
     log_out) = get_LEA_seq_consensus_seqs(fwd_read_f,
                                           rev_read_f,
                                           map_f,
                                           output_dir,
                                           barcode_type,
                                           barcode_len,
                                           barcode_correction_fn,
                                           max_barcode_errors,
                                           min_consensus,
                                           max_cluster_ratio,
                                           min_difference_in_bcs,
                                           fwd_length,
                                           rev_length,
                                           min_reads_per_random_bc,
                                           min_diff_in_clusters,
                                           barcode_column,
                                           reverse_primer_column)

    for sample_id in consensus_seq_lookup:
        for bc_index, rand_bc in enumerate(consensus_seq_lookup[sample_id]):
            consensus_seq = consensus_seq_lookup[sample_id][rand_bc]
            fwd_consensus, rev_consensus = consensus_seq.split('^')
            fwd_consensus_outfile.write(">{}_{}\n{}\n".format(
                sample_id, bc_index, fwd_consensus))
            rev_consensus_outfile.write(">{}_{}\n{}\n".format(
                sample_id, bc_index, rev_consensus))

    log_file.write(log_out)
    log_file.close()

    fwd_read_f.close()
    rev_read_f.close()
    fwd_consensus_outfile.close()
    rev_consensus_outfile.close()
    map_f.close()


if __name__ == "__main__":
    main()
