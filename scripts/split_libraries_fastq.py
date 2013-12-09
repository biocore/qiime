#!/usr/bin/env python
# File created on 07 Jun 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Emily TerAvest"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from os import rename

from cogent import DNA
from cogent.util.misc import safe_md5, create_dir
from qiime.util import parse_command_line_parameters, make_option, gzip_open
from qiime.parse import parse_mapping_file
from qiime.split_libraries_fastq import (process_fastq_single_end_read_file,
    BARCODE_DECODER_LOOKUP, process_fastq_single_end_read_file_no_barcode)
from qiime.split_libraries import check_map
from qiime.split_libraries_fastq import get_illumina_qual_chars
from qiime.golay import get_invalid_golay_barcodes
from qiime.quality import phred_to_ascii33, phred_to_ascii64

phred_to_ascii_fs = {'33':phred_to_ascii33, '64':phred_to_ascii64}

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
    "--sample_id my.sample -o slout_single_sample_q20/ "
    "-m map_not_multiplexed.txt  -q 19 --barcode_type 'not-barcoded'"))

script_info['script_usage'].append(("Quality filter (at Phred >= Q20) one "
    "non-multiplexed lane of Illumina fastq data and write results "
    "to ./slout_single_sample_q20.", "", "%prog -i lane1_read1.fastq.gz "
    "--sample_id my.sample.1 -o slout_single_sample_q20/ "
    "-m map_not_multiplexed.txt -q 19 --barcode_type 'not-barcoded'"))

script_info['script_usage'].append(("Quality filter (at Phred >= Q20) two "
    "non-multiplexed lanes of Illumina fastq data with different samples in "
    "each and write results to ./slout_not_multiplexed_q20.", "", 
    "%prog -i lane1_read1.fastq.gz,lane2_read1.fastq.gz "
    "--sample_id my.sample.1,my.sample.2 -o slout_not_multiplexed_q20/ "
    "-m map_not_multiplexed.txt -q 19 --barcode_type 'not-barcoded'"))

script_info['output_description']= ""
script_info['required_options'] = [
    make_option('-i', '--sequence_read_fps', type="existing_filepaths",
        help='the sequence read fastq files (comma-separated if more than '
        'one)'),
    make_option('-o', '--output_dir', type="new_dirpath", help='directory to '
        'store output files'),
    make_option('-m', '--mapping_fps', type="existing_filepaths",
        help='metadata mapping files (comma-separated if more than one)')
]

script_info['optional_options'] = [
     make_option('-b', '--barcode_read_fps', type="existing_filepaths",
        default=None, help='the barcode read fastq files (comma-separated '
        'if more than one) [default: %default]'),
     make_option("--store_qual_scores", default=False, action='store_true',
        help='store qual strings in .qual files [default: %default]'),
     make_option("--sample_ids", default=None, help='comma-separated list of '
        'samples id to be applied to all sequences, must be one per input '
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
        'golay_12 for golay error-correcting barcodes. Error correction will '
         'only be applied for golay_12 barcodes. [default: %default]',
        default='golay_12'),
     make_option('--max_barcode_errors', default=1.5, type='float',
        help='maximum number of errors in barcode [default: %default]'),
     make_option('--phred_offset', default=None, type="choice",
        choices=phred_to_ascii_fs.keys(), help="the ascii offset to use when "
        "decoding phred scores - warning: in most cases you don't need to "
        "pass this value [default: determined automatically]")
     # NEED TO FIX THIS FUNCTIONALITY - CURRENTLY READING THE WRONG FIELD
     # make_option('--filter_bad_illumina_qual_digit',
     #    action='store_true',\
     #    help='filter sequences which are tagged as not passing the Illumina'+\
     #    ' quality filter [default: %default]',
     #    default=False),\
]

script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

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
    filter_bad_illumina_qual_digit = False #opts.filter_bad_illumina_qual_digit
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
    elif barcode_read_fps is None:
        option_parser.error("Must provide --barcode_read_fps if "
            "--barcode_type is not 'not-barcoded'")

    phred_offset = opts.phred_offset
    if phred_offset is not None:
        try:
            phred_to_ascii_f = phred_to_ascii_fs[phred_offset]
        except KeyError:
            # shouldn't be able to get here, but we'll stay on the
            # safe side
            option_parser.error("Only valid phred offsets are: %s" %
                ' '.join(phred_to_ascii_fs.keys()))
    else:
        # let split_libraries_fastq.process_fastq_single_end_read_file 
        # figure it out...
        phred_to_ascii_f = None

    if opts.last_bad_quality_char is not None:
        option_parser.error('--last_bad_quality_char is no longer supported. '
         'Use -q instead (see option help text by passing -h)')

    if not (0 <= min_per_read_length_fraction <= 1):
        option_parser.error('--min_per_read_length_fraction must be between '
            '0 and 1 (inclusive). You passed %1.5f' %
            min_per_read_length_fraction)

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
    output_f = open(output_fp_temp,'w')
    qual_fp_temp = '%s/qual.fna.incomplete' % output_dir
    qual_fp = '%s/seqs.qual' % output_dir
    output_fastq_fp_temp = '%s/seqs.fastq.incomplete' % output_dir
    output_fastq_fp = '%s/seqs.fastq' % output_dir

    if store_qual_scores:
        qual_f = open(qual_fp_temp,'w')
        # define a qual writer whether we're storing
        # qual strings or not so we don't have to check
        # every time through the for loop below
        def qual_writer(h,q):
            qual_f.write('>%s\n%s\n' % (h,q))
    else:
        def qual_writer(h,q):
            pass

    if store_demultiplexed_fastq:
        output_fastq_f = open(output_fastq_fp_temp,'w')
        # define a fastq writer whether we're storing
        # qual strings or not so we don't have to check
        # every time through the for loop below
        def fastq_writer(h,s,q):
            output_fastq_f.write('@%s\n%s\n+\n%s\n' % (h,s,q))
    else:
        def fastq_writer(h,s,q):
            pass

    log_fp = '%s/split_library_log.txt' % output_dir
    log_f = open(log_fp,'w')
    histogram_fp = '%s/histograms.txt' % output_dir
    histogram_f = open(histogram_fp,'w')

    for i in range(len(sequence_read_fps)):
        sequence_read_fp = sequence_read_fps[i]
        barcode_read_fp = barcode_read_fps[i]
        mapping_fp = mapping_fps[i]
        mapping_f = open(mapping_fp, 'U')
        _, _, barcode_to_sample_id, _, _, _, _ = check_map(mapping_f,
            disable_primer_check=True, has_barcodes=barcode_read_fp is not
            None)

        if rev_comp_mapping_barcodes:
            barcode_to_sample_id = {DNA.rc(k):v for k, v in
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
        log_f.write('Mapping filepath: %s (md5: %s)\n' %
            (mapping_fp,safe_md5(open(mapping_fp)).hexdigest()))
        log_f.write('Sequence read filepath: %s (md5: %s)\n' %
            (sequence_read_fp,
            str(safe_md5(open(sequence_read_fp)).hexdigest())))

        if sequence_read_fp.endswith('.gz'):
            sequence_read_f = gzip_open(sequence_read_fp)
        else:
            sequence_read_f = open(sequence_read_fp,'U')

        seq_id = start_seq_id
        
        if barcode_read_fp is not None:
            log_f.write('Barcode read filepath: %s (md5: %s)\n\n' %
                (barcode_read_fp,safe_md5(open(barcode_read_fp)).hexdigest()))

            if barcode_read_fp.endswith('.gz'):
                barcode_read_f = gzip_open(barcode_read_fp)
            else:
                barcode_read_f = open(barcode_read_fp,'U')

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
               phred_to_ascii_f=phred_to_ascii_f)
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
               phred_to_ascii_f=phred_to_ascii_f)

        for fasta_header, sequence, quality, seq_id in seq_generator:
            output_f.write('>%s\n%s\n' % (fasta_header,sequence))
            qual_writer(fasta_header,quality)
            fastq_writer(fasta_header,sequence,quality)

        start_seq_id = seq_id + 1
        log_f.write('\n---\n\n')

    output_f.close()
    rename(output_fp_temp, output_fp)

    # process the optional output files, as necessary
    if store_qual_scores:
        qual_f.close()
        rename(qual_fp_temp,qual_fp)

    if store_demultiplexed_fastq:
        output_fastq_f.close()
        rename(output_fastq_fp_temp, output_fastq_fp)

if __name__ == "__main__":
    main()
