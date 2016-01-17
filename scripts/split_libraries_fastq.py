#!/usr/bin/env python
# File created on 07 Jun 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Emily TerAvest", "Yoshiki Vazquez Baeza",
               "Rob Knight"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from os import rename

from skbio.util import safe_md5, create_dir
from skbio.sequence import DNA
from skbio.format.sequences import format_fastq_record

from qiime.util import parse_command_line_parameters, make_option, gzip_open
from qiime.parse import parse_mapping_file, parse_items
from qiime.split_libraries_fastq import (process_fastq_single_end_read_file,
                                         BARCODE_DECODER_LOOKUP, process_fastq_single_end_read_file_no_barcode)
from qiime.split_libraries import check_map
from qiime.split_libraries_fastq import get_illumina_qual_chars
from qiime.golay import get_invalid_golay_barcodes

script_info = {}

script_info['brief_description'] = ("This script performs demultiplexing of "
                                    "Fastq sequence data where barcodes and sequences are contained in two "
                                    "separate fastq files (common on Illumina runs).")

script_info['script_description'] = ""

script_info['script_usage'] = []

script_info['script_usage'].append(("Demultiplex and quality filter "
                                    "(at Phred >= Q20) one lane of Illumina fastq data and write results to "
                                    "./slout_q20.", "", "%prog -i lane1_read1.fastq.gz "
                                    "-b lane1_barcode.fastq.gz --rev_comp_mapping_barcodes -o slout_q20/ "
                                    "-m map.txt -q 19"))

script_info['script_usage'].append(("Demultiplex and quality filter "
                                    "(at Phred >= Q20) one lane of Illumina fastq data and write results to "
                                    "./slout_q20. Store trimmed quality scores in addition to sequence data.",
                                    "", "%prog -i lane1_read1.fastq.gz -b lane1_barcode.fastq.gz "
                                    "--rev_comp_mapping_barcodes -o slout_q20/ -m map.txt "
                                    "--store_qual_scores -q 19"))

script_info['script_usage'].append(("Demultiplex and quality filter "
                                    "(at Phred >= Q20) two lanes of Illumina fastq data and write results to "
                                    "./slout_q20.", "", "%prog -i lane1_read1.fastq.gz,lane2_read1.fastq.gz "
                                    "-b lane1_barcode.fastq.gz,lane2_barcode.fastq.gz "
                                    "--rev_comp_mapping_barcodes -o slout_q20/ -m map.txt,map.txt "
                                    "--store_qual_scores -q 19"))

script_info['script_usage'].append(("Quality filter (at Phred >= Q20) one "
                                    "non-multiplexed lane of Illumina fastq data and write results "
                                    "to ./slout_single_sample_q20.", "", "%prog -i lane1_read1.fastq.gz "
                                    "--sample_ids my.sample.1 -o slout_single_sample_q20/ "
                                    "-q 19 --barcode_type 'not-barcoded'"))

script_info['script_usage'].append(("Quality filter (at Phred >= Q20) two "
                                    "non-multiplexed lanes of Illumina fastq data with different samples in "
                                    "each and write results to ./slout_not_multiplexed_q20.", "",
                                    "%prog -i lane1_read1.fastq.gz,lane2_read1.fastq.gz "
                                    "--sample_ids my.sample.1,my.sample.2 -o slout_not_multiplexed_q20/ "
                                    "-q 19 --barcode_type 'not-barcoded'"))

script_info['output_description'] = ""
script_info['required_options'] = [
    make_option('-i', '--sequence_read_fps', type="existing_filepaths",
                help='the sequence read fastq files (comma-separated if more than '
                'one)'),
    make_option('-o', '--output_dir', type="new_dirpath", help='directory to '
                'store output files'),
]

script_info['optional_options'] = [
    make_option('-m', '--mapping_fps', type="existing_filepaths",
                help='metadata mapping files (comma-separated if more than'
                ' one) [default: %default]', default=None),
    make_option('-b', '--barcode_read_fps', type="existing_filepaths",
                default=None, help='the barcode read fastq files (comma-separated '
                'if more than one) [default: %default]'),
    make_option("--store_qual_scores", default=False, action='store_true',
                help='store qual strings in .qual files [default: %default]'),
    make_option("--sample_ids", default=None, help='comma-separated list of '
                'samples ids to be applied to all sequences, must be one per input '
                'file path (used when data is not multiplexed) [default: %default]'),
    make_option("--store_demultiplexed_fastq", default=False,
                action='store_true', help='write demultiplexed fastq files '
                '[default: %default]'),
    make_option("--retain_unassigned_reads", default=False,
                action='store_true', help="retain sequences which don't map to a "
                'barcode in the mapping file (sample ID will be "Unassigned") '
                '[default: %default]'),
    make_option('-r', '--max_bad_run_length', type='int', help='max number '
                'of consecutive low quality base calls allowed before truncating a '
                'read [default: %default]', default=3),
    make_option('-p', '--min_per_read_length_fraction', type='float',
                default=0.75, help='min number of consecutive high quality base calls '
                'to include a read (per single end read) as a fraction of the input '
                'read length [default: %default]'),
    make_option('-n', '--sequence_max_n', type='int', help='maximum number '
                'of N characters allowed in a sequence to retain it -- this is '
                'applied after quality trimming, and is total over combined paired '
                'end reads if applicable [default: %default]', default=0),
    make_option('-s', '--start_seq_id', type='int', help='start seq_ids as '
                'ascending integers beginning with start_seq_id [default: %default]',
                default=0),
    make_option('--rev_comp_barcode', action='store_true', help='reverse '
                'complement barcode reads before lookup [default: %default]',
                default=False),
    make_option('--rev_comp_mapping_barcodes', action='store_true',
                help='reverse complement barcode in mapping before lookup (useful if '
                'barcodes in mapping file are reverse complements of golay codes) '
                '[default: %default]', default=False),
    make_option('--rev_comp', action='store_true', help='reverse complement '
                'sequence before writing to output file (useful for reverse-'
                'orientation reads) [default: %default]', default=False),
    make_option('-q', '--phred_quality_threshold', type='int', help='the '
                'maximum unacceptable Phred quality score (e.g., for Q20 and better, '
                'specify -q 19) [default: %default]', default=3),
    make_option('--last_bad_quality_char', help='DEPRECATED: use -q instead. '
                'This method of setting is not robust to different versions of '
                'CASAVA.'),
    make_option('--barcode_type', type='string', help='The type of barcode '
                'used. This can be an integer, e.g. for length 6 barcodes, or '
                '"golay_12" for golay error-correcting barcodes. Error correction will '
                'only be applied for "golay_12" barcodes. If data is not barcoded, pass '
                '"not-barcoded". [default: %default]',
                default='golay_12'),
    make_option('--max_barcode_errors', default=1.5, type='float',
                help='maximum number of errors in barcode [default: %default]'),
    make_option('--phred_offset', default=None, type="choice",
                choices=['33', '64'], help="the ascii offset to use when "
                "decoding phred scores (either 33 or 64). Warning: in most "
                "cases you don't need to pass this value "
                "[default: determined automatically]"),
    make_option('--read_arguments_from_file', default=False,
                action='store_true', help='If this flag is enabled, then the '
                'inputs to "-i" or "--sequence_read_fps", "-b" or '
                '"--barcode_read_fps", "-m" or "--mapping_fps" and '
                '"--sample_ids" will each be interpreted as a single text file'
                ', where the contents are one file-path or sample identifier '
                'per line (depending on the flag). NOTE: In most cases regular'
                ' users don\'t need to use this flag, as it is intended for '
                'use in multiple_split_libraries_fastq.py [default: %default]')
    # NEED TO FIX THIS FUNCTIONALITY - CURRENTLY READING THE WRONG FIELD
    # make_option('--filter_bad_illumina_qual_digit',
    #    action='store_true',
    #    help='filter sequences which are tagged as not passing the Illumina'+
    #    ' quality filter [default: %default]',
    #    default=False),
]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    read_arguments_from_file = opts.read_arguments_from_file

    # these arguments can optionally be read from a file, reasoning is to
    # allow arguments that would span over hundreds of samples and would be
    # prohibitive to execute as a command line call
    if read_arguments_from_file:
        # sample_ids is the only one of these arguments that's returned as a
        # string, the rest of them are lists
        if opts.sample_ids:
            opts.sample_ids = ','.join(parse_items(opts.sample_ids))
        if opts.sequence_read_fps:
            opts.sequence_read_fps = parse_items(opts.sequence_read_fps[0])
        if opts.barcode_read_fps:
            opts.barcode_read_fps = parse_items(opts.barcode_read_fps[0])
        if opts.mapping_fps:
            opts.mapping_fps = parse_items(opts.mapping_fps[0])

    sequence_read_fps = opts.sequence_read_fps
    barcode_read_fps = opts.barcode_read_fps
    sample_ids = None
    if opts.sample_ids is not None:
        sample_ids = opts.sample_ids.split(',')
    mapping_fps = opts.mapping_fps
    phred_quality_threshold = opts.phred_quality_threshold
    retain_unassigned_reads = opts.retain_unassigned_reads
    min_per_read_length_fraction = opts.min_per_read_length_fraction
    max_bad_run_length = opts.max_bad_run_length
    rev_comp = opts.rev_comp
    rev_comp_barcode = opts.rev_comp_barcode
    rev_comp_mapping_barcodes = opts.rev_comp_mapping_barcodes
    seq_max_N = opts.sequence_max_n
    start_seq_id = opts.start_seq_id
    # NEED TO FIX THIS FUNCTIONALITY - CURRENTLY READING THE WRONG FIELD
    # opts.filter_bad_illumina_qual_digit
    filter_bad_illumina_qual_digit = False
    store_qual_scores = opts.store_qual_scores
    store_demultiplexed_fastq = opts.store_demultiplexed_fastq
    barcode_type = opts.barcode_type
    max_barcode_errors = opts.max_barcode_errors

    # if this is not a demultiplexed run,
    if barcode_type == 'not-barcoded':
        if sample_ids is None:
            option_parser.error("If not providing barcode reads (because "
                                "your data is not multiplexed), must provide --sample_ids.")
        if len(sample_ids) != len(sequence_read_fps):
            option_parser.error("If providing --sample_ids (because "
                                "your data is not multiplexed), must provide the same number "
                                "of sample ids as sequence read filepaths.")
        barcode_read_fps = [None] * len(sequence_read_fps)
        mapping_fps = [None] * len(sequence_read_fps)
    elif barcode_read_fps is None:
        option_parser.error("Must provide --barcode_read_fps if "
                            "--barcode_type is not 'not-barcoded'")
    elif mapping_fps is None:
        option_parser.error("Must provide --mapping_fps if "
                            "--barcode_type is not 'not-barcoded'")

    phred_offset = opts.phred_offset
    if phred_offset is not None:
        try:
            phred_offset = int(phred_offset)
        except ValueError:
            # shouldn't be able to get here...
            option_parser.error(
                "If --phred_offset is provided, it must be a valid integer.")

    if opts.last_bad_quality_char is not None:
        option_parser.error('--last_bad_quality_char is no longer supported. '
                            'Use -q instead (see option help text by passing -h)')

    if not (0 < min_per_read_length_fraction <= 1):
        option_parser.error('--min_per_read_length_fraction must be greater '
                            'than 0 and less than or equal to 1. You passed '
                            '%1.5f.' % min_per_read_length_fraction)

    barcode_correction_fn = BARCODE_DECODER_LOOKUP.get(barcode_type, None)

    if len(mapping_fps) == 1 and len(sequence_read_fps) > 1:
        mapping_fps = mapping_fps * len(sequence_read_fps)

    if len(set([len(sequence_read_fps), len(barcode_read_fps),
                len(mapping_fps)])) > 1:
        option_parser.error("Same number of sequence, barcode, and mapping "
                            "files must be provided.")

    output_dir = opts.output_dir
    create_dir(output_dir)

    output_fp_temp = '%s/seqs.fna.incomplete' % output_dir
    output_fp = '%s/seqs.fna' % output_dir
    output_f = open(output_fp_temp, 'w')
    qual_fp_temp = '%s/qual.fna.incomplete' % output_dir
    qual_fp = '%s/seqs.qual' % output_dir
    output_fastq_fp_temp = '%s/seqs.fastq.incomplete' % output_dir
    output_fastq_fp = '%s/seqs.fastq' % output_dir

    if store_qual_scores:
        qual_f = open(qual_fp_temp, 'w')
        # define a qual writer whether we're storing
        # qual strings or not so we don't have to check
        # every time through the for loop below

        def qual_writer(h, q):
            qual_f.write('>%s\n%s\n' % (h, q))
    else:
        def qual_writer(h, q):
            pass

    if store_demultiplexed_fastq:
        output_fastq_f = open(output_fastq_fp_temp, 'w')
        # define a fastq writer whether we're storing
        # qual strings or not so we don't have to check
        # every time through the for loop below

        def fastq_writer(h, s, q):
            output_fastq_f.write(format_fastq_record(h, s, q))
    else:
        def fastq_writer(h, s, q):
            pass

    log_fp = '%s/split_library_log.txt' % output_dir
    log_f = open(log_fp, 'w')
    histogram_fp = '%s/histograms.txt' % output_dir
    histogram_f = open(histogram_fp, 'w')

    for i in range(len(sequence_read_fps)):
        sequence_read_fp = sequence_read_fps[i]
        barcode_read_fp = barcode_read_fps[i]
        mapping_fp = mapping_fps[i]
        if mapping_fp is not None:
            mapping_f = open(mapping_fp, 'U')
            _, _, barcode_to_sample_id, _, _, _, _ = check_map(mapping_f,
                disable_primer_check=True,
                has_barcodes=barcode_read_fp is not None)
        else:
            mapping_f = None
            barcode_to_sample_id = {}

        if rev_comp_mapping_barcodes:
            barcode_to_sample_id = {str(DNA(k).rc()): v for k, v in
                                    barcode_to_sample_id.iteritems()}

        if barcode_type == 'golay_12':
            invalid_golay_barcodes = get_invalid_golay_barcodes(
                barcode_to_sample_id.keys())
            if len(invalid_golay_barcodes) > 0:
                option_parser.error("Some or all barcodes are not valid golay "
                                    "codes. Do they need to be reverse complemented? If these "
                                    "are not golay barcodes pass --barcode_type 12 to disable "
                                    "barcode error correction, or pass --barcode_type # if "
                                    "the barcodes are not 12 base pairs, where # is the size "
                                    "of the barcodes. Invalid codes:\n\t%s" %
                                    ' '.join(invalid_golay_barcodes))

        log_f.write("Input file paths\n")
        if mapping_fp is not None:
            log_f.write('Mapping filepath: %s (md5: %s)\n' %
                        (mapping_fp, safe_md5(open(mapping_fp)).hexdigest()))
        log_f.write('Sequence read filepath: %s (md5: %s)\n' %
                    (sequence_read_fp,
                     str(safe_md5(open(sequence_read_fp)).hexdigest())))

        if sequence_read_fp.endswith('.gz'):
            sequence_read_f = gzip_open(sequence_read_fp)
        else:
            sequence_read_f = open(sequence_read_fp, 'U')

        seq_id = start_seq_id

        if barcode_read_fp is not None:
            log_f.write('Barcode read filepath: %s (md5: %s)\n\n' %
                        (barcode_read_fp,
                         safe_md5(open(barcode_read_fp)).hexdigest()))

            if barcode_read_fp.endswith('.gz'):
                barcode_read_f = gzip_open(barcode_read_fp)
            else:
                barcode_read_f = open(barcode_read_fp, 'U')

            seq_generator = process_fastq_single_end_read_file(
                sequence_read_f, barcode_read_f, barcode_to_sample_id,
                store_unassigned=retain_unassigned_reads,
                max_bad_run_length=max_bad_run_length,
                phred_quality_threshold=phred_quality_threshold,
                min_per_read_length_fraction=min_per_read_length_fraction,
                rev_comp=rev_comp, rev_comp_barcode=rev_comp_barcode,
                seq_max_N=seq_max_N, start_seq_id=start_seq_id,
                filter_bad_illumina_qual_digit=filter_bad_illumina_qual_digit,
                log_f=log_f, histogram_f=histogram_f,
                barcode_correction_fn=barcode_correction_fn,
                max_barcode_errors=max_barcode_errors,
                phred_offset=phred_offset)
        else:
            seq_generator = process_fastq_single_end_read_file_no_barcode(
                sequence_read_f, sample_ids[i],
                store_unassigned=retain_unassigned_reads,
                max_bad_run_length=max_bad_run_length,
                phred_quality_threshold=phred_quality_threshold,
                min_per_read_length_fraction=min_per_read_length_fraction,
                rev_comp=rev_comp, seq_max_N=seq_max_N,
                start_seq_id=start_seq_id,
                filter_bad_illumina_qual_digit=filter_bad_illumina_qual_digit,
                log_f=log_f, histogram_f=histogram_f,
                phred_offset=phred_offset)

        for fasta_header, sequence, quality, seq_id in seq_generator:
            output_f.write('>%s\n%s\n' % (fasta_header, sequence))
            qual_writer(fasta_header, quality)
            fastq_writer(fasta_header, sequence, quality)

        start_seq_id = seq_id + 1
        log_f.write('\n---\n\n')

    output_f.close()
    rename(output_fp_temp, output_fp)

    # process the optional output files, as necessary
    if store_qual_scores:
        qual_f.close()
        rename(qual_fp_temp, qual_fp)

    if store_demultiplexed_fastq:
        output_fastq_f.close()
        rename(output_fastq_fp_temp, output_fastq_fp)

if __name__ == "__main__":
    main()
