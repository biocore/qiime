#!/usr/bin/env python
# File created on 22 Mar 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

from itertools import izip
from os.path import split, splitext
from os import makedirs
from numpy import log10, arange, histogram
from cogent import DNA
from cogent.parse.fastq import MinimalFastqParser
from qiime.format import (format_histogram_one_count,
                          format_split_libraries_fastq_log)
from qiime.hamming import decode_hamming_8
from qiime.golay import decode_golay_12

class FastqParseError(Exception):
    pass

def get_illumina_qual_chars():
    return '@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'

def bad_chars_from_threshold(first_bad_char):
    if first_bad_char == '':
        return {}
    else:
        all_chars = get_illumina_qual_chars()
        first_bad_char_index = all_chars.index(first_bad_char)
        bad_chars = list(get_illumina_qual_chars()[:first_bad_char_index+1])
        return {}.fromkeys(bad_chars)

def read_qual_score_filter(seq,qual,max_run_length,threshold):
    """slices illumina sequence and quality line based on quality filter
    """
    bad_chars = bad_chars_from_threshold(threshold)
    last_good_slice_end_pos = 0
    bad_run_length = 0
    for i in range(len(seq)):
        if qual[i] in bad_chars:
            bad_run_length += 1
        else:
            bad_run_length = 0
            last_good_slice_end_pos = i + 1
            
        if bad_run_length > max_run_length:
            return seq[:last_good_slice_end_pos],\
                   qual[:last_good_slice_end_pos]
    
    # There were no runs that were too bad for too long 
    return seq, qual

def quality_filter_sequence(header,
                            sequence,
                            quality,
                            max_bad_run_length,
                            last_bad_quality_char,
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
            quality_char = header[h.index('#')+1]
            illumina_quality_digit = int(quality_char)
        except ValueError:
            pass
        else:
            if illumina_quality_digit == 0:
                return 3, sequence, quality
        
    sequence, quality = read_qual_score_filter(sequence,
                                       quality,
                                       max_bad_run_length, 
                                       last_bad_quality_char)
                                       
    if (len(sequence) < min_per_read_length):
        return 1, sequence, quality
    elif (sequence.count('N') > seq_max_N):
        return 2, sequence, quality
    else:
        return 0, sequence, quality

def check_header_match(header1,header2):
    
    # split on '#' and '/' to handle cases with and without the
    # Illumina quality digit
    header1 = header1.split('#')[0].split('/')[0]
    header2 = header2.split('#')[0].split('/')[0]
    
    return header1 == header2
    
BARCODE_DECODER_LOOKUP = {
 'golay_12':decode_golay_12,
 #'hamming_8':decode_hamming_8,
}
    
def correct_barcode(barcode,barcode_to_sample_id,correction_fn):
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
    
    if sample_id != None or correction_fn == None or 'N' in barcode:
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


def process_fastq_single_end_read_file(fastq_read_f,
                                       fastq_barcode_f,
                                       barcode_to_sample_id,
                                       store_unassigned=False,
                                       max_bad_run_length=0,
                                       last_bad_quality_char='B',
                                       min_per_read_length=75,
                                       rev_comp=False,
                                       rev_comp_barcode=False,
                                       seq_max_N=0,
                                       start_seq_id=0,
                                       filter_bad_illumina_qual_digit=True,
                                       log_f=None,
                                       histogram_f=None,
                                       barcode_correction_fn=None,
                                       max_barcode_errors=1.5):
    """parses fastq single-end read file
    """
    header_index = 0
    sequence_index = 1
    quality_index = 2
    
    seq_id = start_seq_id
    
    # prep data for logging
    input_sequence_count = 0
    count_barcode_not_in_map = 0
    count_too_short = 0
    count_too_many_N = 0
    count_bad_illumina_qual_digit = 0
    count_barcode_errors_exceed_max = 0
    sequence_lengths = []
    seqs_per_sample_counts = {}
    
    for bc_data,read_data in izip(MinimalFastqParser(fastq_barcode_f,strict=False),
                                  MinimalFastqParser(fastq_read_f,strict=False)):
        input_sequence_count += 1
        # Confirm match between barcode and read headers
        if not check_header_match(bc_data[header_index],
                                  read_data[header_index]):
            raise FastqParseError,\
             ("Headers of barcode and read do not match. Can't continue. "
              "Confirm that the barcode fastq and read fastq that you are "
              "passing match one another.")
        else:
            header = read_data[header_index]
        
        # Grab the barcode sequence
        barcode = bc_data[sequence_index]
        if rev_comp_barcode:
            barcode = DNA.rc(barcode)
        # Grab the read sequence
        sequence = read_data[1]
        # Grab the read quality
        quality = read_data[2]
        
        # correct the barcode (if applicable) and map to sample id
        num_barcode_errors, corrected_barcode, correction_attempted, sample_id = \
         correct_barcode(barcode,barcode_to_sample_id,barcode_correction_fn)
        # skip samples with too many errors
        if (num_barcode_errors > max_barcode_errors):
          count_barcode_errors_exceed_max += 1
          continue
        
        # skip unassignable samples unless otherwise requested
        if sample_id == None:
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
                                  last_bad_quality_char,
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
                raise ValueError,\
                 "Unknown quality filter result: %d" % quality_filter_result
            continue
        
        sequence_lengths.append(len(sequence))
        
        try:
            seqs_per_sample_counts[sample_id] += 1
        except KeyError:
            seqs_per_sample_counts[sample_id] = 1
        
        if rev_comp:
            sequence = DNA.rc(sequence)
            quality = quality[::-1]
        
        fasta_header = '%s_%s %s orig_bc=%s new_bc=%s bc_diffs=%d' %\
          (sample_id,seq_id,header,barcode,corrected_barcode,num_barcode_errors)
        yield fasta_header, sequence, quality, seq_id
        seq_id += 1

    if log_f != None:
        log_str = format_split_libraries_fastq_log(count_barcode_not_in_map,
                             count_too_short,
                             count_too_many_N,
                             count_bad_illumina_qual_digit,
                             count_barcode_errors_exceed_max,
                             input_sequence_count,
                             sequence_lengths,
                             seqs_per_sample_counts)
        log_f.write(log_str)
    
    if len(sequence_lengths) and histogram_f != None:
        counts, bin_edges = make_histograms(sequence_lengths)
        histogram_str = format_histogram_one_count(counts,bin_edges)
        histogram_f.write(histogram_str)
        histogram_f.write('\n--\n\n')
    
def make_histograms(lengths, binwidth=10):
    """Makes histogram data for pre and post lengths"""
    min_len = min(lengths)
    max_len = max(lengths)
    floor = (min_len/binwidth)*binwidth
    ceil = ((max_len/binwidth)+2)*binwidth
    bins = arange(floor, ceil, binwidth)
    hist, bin_edges = histogram(lengths,bins)
    return hist, bin_edges
