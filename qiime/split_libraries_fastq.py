#!/usr/bin/env python
# File created on 22 Mar 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Jose Antonio Navas Molina", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from itertools import izip, cycle
from os.path import split, splitext, join
from os import makedirs

import numpy as np

from skbio.parse.sequences import parse_fastq
from skbio.sequence import DNA
from skbio.format.sequences import format_fastq_record

from qiime.format import (format_histogram_one_count,
                          format_split_libraries_fastq_log)
from qiime.parse import is_casava_v180_or_later
from qiime.hamming import decode_hamming_8
from qiime.golay import decode_golay_12
from qiime.util import qiime_open


class FastqParseError(Exception):
    pass


def get_illumina_qual_chars():
    # pulled from stack overflow (url wrapped over two lines):
    # http://stackoverflow.com/questions/5891453/
    # is-there-a-python-library-that-contains-a-list-of-all-the-ascii-characters
    return ['\t', '\n', '\r', ' ', '!', '"', '#', '$', '%', '&', "'", '(', ')',
            '*', '+', ',', '-', '.', '/', '0', '1', '2', '3', '4', '5', '6', '7',
            '8', '9', ':', ';', '<', '=', '>', '?', '@', 'A', 'B', 'C', 'D', 'E',
            'F', 'G', 'H', 'I', 'J', 'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S',
            'T', 'U', 'V', 'W', 'X', 'Y', 'Z', '[', '\\', ']', '^', '_', '`', 'a',
            'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j', 'k', 'l', 'm', 'n', 'o',
            'p', 'q', 'r', 's', 't', 'u', 'v', 'w', 'x', 'y', 'z', '{', '|', '}',
            '~']


def bad_chars_from_threshold(first_bad_char):
    if first_bad_char == '':
        return {}
    else:
        all_chars = get_illumina_qual_chars()
        first_bad_char_index = all_chars.index(first_bad_char)
        bad_chars = list(get_illumina_qual_chars()[:first_bad_char_index + 1])
        return {}.fromkeys(bad_chars)

def _contiguous_regions(condition):
    """Finds contiguous True regions of the boolean array "condition". Returns
    a 2D array where the first column is the start index of the region and the
    second column is the end index.

    This function was derived from SO:
    http://stackoverflow.com/a/4495197/579416
    """

    # Find the indicies of changes in "condition"
    d = np.diff(condition)
    idx, = d.nonzero()

    # We need to start things after the change in "condition". Therefore,
    # we'll shift the index by 1 to the right.
    idx += 1

    if condition[0]:
        # If the start of condition is True prepend a 0
        idx = np.r_[0, idx]

    if condition[-1]:
        # If the end of condition is True, append the length of the array
        idx = np.r_[idx, condition.size]

    # Reshape the result into two columns
    idx.shape = (-1,2)
    return idx

def read_qual_score_filter(seq, qual, max_run_length, threshold):
    """slices illumina sequence and quality line based on quality filter
    """
    last_good_slice_end_pos = 0
    bad_run_length = 0
    mask = qual <= threshold
    for starts, ends in _contiguous_regions(mask):
        if ends - starts > max_run_length:
            return seq[:starts], qual[:starts]

    return seq, qual


def quality_filter_sequence(header,
                            sequence,
                            quality,
                            max_bad_run_length,
                            phred_quality_threshold,
                            min_per_read_length,
                            seq_max_N,
                            filter_bad_illumina_qual_digit):
    if filter_bad_illumina_qual_digit:
        h = header.split()[0]
        try:
            # this block is a little strange because each of these
            # can throw a ValueError. The same thing needs to be done
            # in either case, so it doesn't really make sense to split
            # into two separate try/excepts, particulary because that would
            # complicate the logic
            quality_char = header[h.index('#') + 1]
            illumina_quality_digit = int(quality_char)
        except ValueError:
            pass
        else:
            if illumina_quality_digit == 0:
                return 3, sequence, quality

    sequence, quality = read_qual_score_filter(sequence,
                                               quality,
                                               max_bad_run_length,
                                               phred_quality_threshold)

    if (len(sequence) < min_per_read_length):
        return 1, sequence, quality
    elif (sequence.count('N') > seq_max_N):
        return 2, sequence, quality
    else:
        return 0, sequence, quality


def check_header_match_pre180(header1, header2):

    # split on '#' and '/' to handle cases with and without the
    # Illumina quality digit
    header1 = header1.split('#')[0].split('/')[0]
    header2 = header2.split('#')[0].split('/')[0]

    return header1 == header2


def check_header_match_180_or_later(header1, header2):
    """ Confirm headers are compatible in CASAVA 1.8.0 or later

        These contain information on the read number, so can differ
    """
    header1 = header1.split(':')
    header2 = header2.split(':')
    for e1, e2 in zip(header1, header2):
        if e1.split(' ')[0] != e2.split(' ')[0]:
            return False

    return True

BARCODE_DECODER_LOOKUP = {
    'golay_12': decode_golay_12,
    #'hamming_8':decode_hamming_8,
}


def correct_barcode(barcode, barcode_to_sample_id, correction_fn):
    """Correct barcode given barcode, dict of valid barcodes, and correction fn

       return value: (number of errors,
                      corrected barcode,
                      correction was attempted (bool),
                      sample id [None if can't be determined])

    """
    # Map the barcode if possible
    try:
        sample_id = barcode_to_sample_id[barcode]
    except KeyError:
        sample_id = None

    if sample_id is not None or correction_fn is None or 'N' in barcode:
        # barcode isn't corrected, either because is maps directly to
        # a sample ID, because we're not correcting barcodes, or because
        # it contains an 'N' character
        return 0, barcode, False, sample_id
    else:
        # correct the barcode
        corrected_barcode, num_errors = correction_fn(barcode)
        try:
            sample_id = barcode_to_sample_id[corrected_barcode]
        except KeyError:
            sample_id = None

        return num_errors, corrected_barcode, True, sample_id


def process_fastq_single_end_read_file_no_barcode(
        fastq_read_f,
        sample_id,
        store_unassigned=False,
        max_bad_run_length=0,
        phred_quality_threshold=2,
        min_per_read_length_fraction=0.75,
        rev_comp=False,
        seq_max_N=0,
        start_seq_id=0,
        filter_bad_illumina_qual_digit=False,
        log_f=None,
        histogram_f=None,
        phred_offset=None):
    """ Quality filtering when a single sample has been run in a lane

        This code simulates a barcode file to allow us to re-use the quality
         filtering code when not demultiplexing. A post-split-libraries file
         is generated by assigning all sequences to sample_id
    """
    # simulate a barcode fastq file
    fake_barcodes = cycle(["@", "AAAAAAAAAAAA", "+", "CCCCCCCCCCCC"])
    # make a fake barcode mapping
    barcode_to_sample_id = {"AAAAAAAAAAAA": sample_id}
    for e in process_fastq_single_end_read_file(
            fastq_read_f,
            fake_barcodes,
            barcode_to_sample_id,
            store_unassigned=store_unassigned,
            max_bad_run_length=max_bad_run_length,
            phred_quality_threshold=phred_quality_threshold,
            min_per_read_length_fraction=min_per_read_length_fraction,
            rev_comp=rev_comp,
            rev_comp_barcode=False,
            seq_max_N=seq_max_N,
            start_seq_id=start_seq_id,
            filter_bad_illumina_qual_digit=filter_bad_illumina_qual_digit,
            log_f=log_f,
            histogram_f=histogram_f,
            barcode_correction_fn=None,
            max_barcode_errors=0,
            strict_header_match=False,
            phred_offset=phred_offset):
        yield e


def process_fastq_single_end_read_file(fastq_read_f,
                                       fastq_barcode_f,
                                       barcode_to_sample_id,
                                       store_unassigned=False,
                                       max_bad_run_length=0,
                                       phred_quality_threshold=2,
                                       min_per_read_length_fraction=0.75,
                                       rev_comp=False,
                                       rev_comp_barcode=False,
                                       seq_max_N=0,
                                       start_seq_id=0,
                                       filter_bad_illumina_qual_digit=False,
                                       log_f=None,
                                       histogram_f=None,
                                       barcode_correction_fn=None,
                                       max_barcode_errors=1.5,
                                       strict_header_match=True,
                                       phred_offset=None):
    """parses fastq single-end read file
    """
    header_index = 0
    sequence_index = 1
    quality_index = 2

    seq_id = start_seq_id
    # grab the first lines and then seek back to the beginning of the file
    try:
        fastq_read_f_line1 = fastq_read_f.readline()
        fastq_read_f_line2 = fastq_read_f.readline()
        fastq_read_f.seek(0)
    except AttributeError:
        fastq_read_f_line1 = fastq_read_f[0]
        fastq_read_f_line2 = fastq_read_f[1]

    if phred_offset is None:
        post_casava_v180 = is_casava_v180_or_later(fastq_read_f_line1)
        if post_casava_v180:
            phred_offset = 33
        else:
            phred_offset = 64

    if phred_offset == 33:
        check_header_match_f = check_header_match_180_or_later
    elif phred_offset == 64:
        check_header_match_f = check_header_match_pre180
    else:
        raise ValueError("Invalid PHRED offset: %d" % phred_offset)

    # compute the barcode length, if they are all the same.
    # this is useful for selecting a subset of the barcode read
    # if it's too long (e.g., for technical reasons on the sequencer)
    barcode_lengths = set([len(bc)
                          for bc, sid in barcode_to_sample_id.items()])
    if len(barcode_lengths) == 1:
        barcode_length = barcode_lengths.pop()
    else:
        barcode_length = None

    # compute the minimum read length as a fraction of the length of the input
    # read
    min_per_read_length = min_per_read_length_fraction * \
        len(fastq_read_f_line2)

    # prep data for logging
    input_sequence_count = 0
    count_barcode_not_in_map = 0
    count_too_short = 0
    count_too_many_N = 0
    count_bad_illumina_qual_digit = 0
    count_barcode_errors_exceed_max = 0
    sequence_lengths = []
    seqs_per_sample_counts = {}
    for bc_data, read_data in izip(
            parse_fastq(fastq_barcode_f, strict=False, phred_offset=phred_offset),
            parse_fastq(fastq_read_f, strict=False, phred_offset=phred_offset)):
        input_sequence_count += 1
        # Confirm match between barcode and read headers
        if strict_header_match and \
           (not check_header_match_f(bc_data[header_index], read_data[header_index])):
            raise FastqParseError("Headers of barcode and read do not match. Can't continue. "
                                  "Confirm that the barcode fastq and read fastq that you are "
                                  "passing match one another.")
        else:
            header = read_data[header_index]

        # Grab the barcode sequence
        if barcode_length:
            # because thirteen cycles are sometimes used for
            # techical reasons, this step looks only at the
            # first tweleve bases. note that the barcode is
            # rev-comp'ed after this step if requested since
            # the thirteen base is a technical artefact, not
            # barcode sequence.
            barcode = bc_data[sequence_index][:barcode_length]
        else:
            barcode = bc_data[sequence_index]
        if rev_comp_barcode:
            barcode = str(DNA(barcode).rc())
        # Grab the read sequence
        sequence = read_data[1]
        # Grab the read quality
        quality = read_data[2]

        # correct the barcode (if applicable) and map to sample id
        num_barcode_errors, corrected_barcode, correction_attempted, sample_id = \
            correct_barcode(
                barcode,
                barcode_to_sample_id,
                barcode_correction_fn)
        # skip samples with too many errors
        if (num_barcode_errors > max_barcode_errors):
            count_barcode_errors_exceed_max += 1
            continue

        # skip unassignable samples unless otherwise requested
        if sample_id is None:
            if not store_unassigned:
                count_barcode_not_in_map += 1
                continue
            else:
                sample_id = 'Unassigned'

        quality_filter_result, sequence, quality =\
            quality_filter_sequence(header,
                                    sequence,
                                    quality,
                                    max_bad_run_length,
                                    phred_quality_threshold,
                                    min_per_read_length,
                                    seq_max_N,
                                    filter_bad_illumina_qual_digit)

        # process quality result
        if quality_filter_result != 0:
            # if the quality filter didn't pass record why and
            # move on to the next record
            if quality_filter_result == 1:
                count_too_short += 1
            elif quality_filter_result == 2:
                count_too_many_N += 1
            elif quality_filter_result == 3:
                count_bad_illumina_qual_digit += 1
            else:
                raise ValueError(
                    "Unknown quality filter result: %d" %
                    quality_filter_result)
            continue

        sequence_lengths.append(len(sequence))

        try:
            seqs_per_sample_counts[sample_id] += 1
        except KeyError:
            seqs_per_sample_counts[sample_id] = 1

        if rev_comp:
            sequence = str(DNA(sequence).rc())
            quality = quality[::-1]

        fasta_header = '%s_%s %s orig_bc=%s new_bc=%s bc_diffs=%d' %\
            (sample_id, seq_id, header, barcode,
             corrected_barcode, num_barcode_errors)
        yield fasta_header, sequence, quality, seq_id
        seq_id += 1

    # Add sample IDs with zero counts to dictionary for logging
    for curr_sample_id in barcode_to_sample_id.values():
        if curr_sample_id not in seqs_per_sample_counts.keys():
            seqs_per_sample_counts[curr_sample_id] = 0

    if log_f is not None:
        log_str = format_split_libraries_fastq_log(count_barcode_not_in_map,
                                                   count_too_short,
                                                   count_too_many_N,
                                                   count_bad_illumina_qual_digit,
                                                   count_barcode_errors_exceed_max,
                                                   input_sequence_count,
                                                   sequence_lengths,
                                                   seqs_per_sample_counts)
        log_f.write(log_str)

    if len(sequence_lengths) and histogram_f is not None:
        counts, bin_edges = make_histograms(sequence_lengths)
        histogram_str = format_histogram_one_count(counts, bin_edges)
        histogram_f.write(histogram_str)
        histogram_f.write('\n--\n\n')


def make_histograms(lengths, binwidth=10):
    """Makes histogram data for pre and post lengths"""
    min_len = min(lengths)
    max_len = max(lengths)
    floor = (min_len / binwidth) * binwidth
    ceil = ((max_len / binwidth) + 2) * binwidth
    bins = np.arange(floor, ceil, binwidth)
    hist, bin_edges = np.histogram(lengths, bins)
    return hist, bin_edges


def extract_reads_from_interleaved(
        input_fp, forward_id, reverse_id, output_dir):
    """Parses a single fastq file and creates two new files: forward and reverse, based on
    the two values (comma separated) in read_direction_identifiers

    input_fp: file path to input
    read_direction_identifiers: comma separated values to identify forward and reverse reads
    output_folder: file path to the output folder
    """
    forward_fp = join(output_dir, "forward_reads.fastq")
    reverse_fp = join(output_dir, "reverse_reads.fastq")
    ffp = open(forward_fp, 'w')
    rfp = open(reverse_fp, 'w')

    for label, seq, qual in parse_fastq(qiime_open(input_fp), strict=False,
                                        enforce_qual_range=False):
        fastq_string = format_fastq_record(label, seq, qual)
        if forward_id in label:
            ffp.write(fastq_string)
        elif reverse_id in label and forward_id not in label:
            rfp.write(fastq_string)
        else:
            ffp.close()
            rfp.close()
            raise ValueError("One of the input sequences doesn't have either identifier "
                             "or it has both.\nLabel: %s\nForward: %s\n Reverse: %s" %
                             (label, forward_id, reverse_id))
    ffp.close()
    rfp.close()
