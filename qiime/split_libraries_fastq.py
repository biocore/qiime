#!/usr/bin/env python
# File created on 22 Mar 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.2.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from itertools import izip
from numpy import log10
from cogent import DNA
from cogent.parse.fastq import MinimalFastqParser
from os.path import split, splitext
from os import makedirs

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
                            quality_threshold,
                            min_per_read_length,
                            seq_max_N,
                            filter_bad_illumina_qual_digit):
    if filter_bad_illumina_qual_digit:
        illumina_quality_digit = int(header.split()[0].split('/')[1])
        if illumina_quality_digit == 0:
            return 3, sequence, quality
        
    sequence, quality = read_qual_score_filter(sequence,
                                       quality,
                                       max_bad_run_length, 
                                       quality_threshold)
                                       
    if (len(sequence) < min_per_read_length):
        return 1, sequence, quality
    elif (sequence.count('N') > seq_max_N):
        return 2, sequence, quality
    else:
        return 0, sequence, quality
    
def extract_barcode(header,
                    sequence,
                    quality,
                    barcode_length,
                    rev_comp_barcode,
                    barcode_in_seq):
    if barcode_in_seq:
        barcode = sequence[:barcode_length]
        sequence = sequence[barcode_length:]
        quality = quality[barcode_length:]
    else:
        try:
            barcode = header.split()[1]
        except IndexError:
            raise IndexError, \
             ("Barcode is not present in header - is it included in "
              "the sequence?\n Bad header: %s" % header)
        
    if rev_comp_barcode:
        barcode = DNA.rc(barcode)
    
    return barcode, sequence, quality


def process_fastq_single_end_read_file(fastq_f,
                                       barcode_to_sample_id,
                                       barcode_length,
                                       store_unassigned=False,
                                       max_bad_run_length=0,
                                       quality_threshold='B',
                                       min_per_read_length=75,
                                       rev_comp=False,
                                       rev_comp_barcode=False,
                                       barcode_in_seq=False,
                                       seq_max_N=0,
                                       start_seq_id=0,
                                       filter_bad_illumina_qual_digit=True):
    """parses fastq single-end read file
    """
    
    seq_id = start_seq_id
    
    count_not_barcode_in_map = 0
    count_too_short = 0
    count_too_many_N = 0
    
    for header,sequence,quality in MinimalFastqParser(fastq_f,strict=False):
        barcode, sequence, quality = extract_barcode(header,
                                                      sequence,
                                                      quality,
                                                      barcode_length,
                                                      rev_comp_barcode,
                                                      barcode_in_seq)
        try:
          sample_id = barcode_to_sample_id[barcode]
        except KeyError:
          if not store_unassigned:
              count_not_barcode_in_map += 1
              continue
          else:
              sample_id = 'Unassigned'
        
        quality_filter_result, sequence, quality =\
          quality_filter_sequence(header,
                                  sequence,
                                  quality,
                                  max_bad_run_length,
                                  quality_threshold,
                                  min_per_read_length,
                                  seq_max_N,
                                  filter_bad_illumina_qual_digit)
        if quality_filter_result != 0:
            # if the quality filter didn't pass record why and 
            # move on to the next record
            if quality_filter_result == 1:
                count_too_short += 1
            elif quality_filter_result == 2:
                count_too_many_N += 1
            else:
                raise ValueError,\
                 "Unknown quality filter result: %d" % quality_filter_result
            continue
        
        if rev_comp:
            sequence = DNA.rc(sequence)
            quality = quality[::-1]
        
        fasta_header = '%s_%s %s' % (sample_id,seq_id,header)
        yield fasta_header, sequence, quality, seq_id
        seq_id += 1

