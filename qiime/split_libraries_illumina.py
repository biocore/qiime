#!/usr/bin/env python
# File created on 22 Mar 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from itertools import izip
from numpy import log10
from cogent.util.misc import revComp
from os.path import split, splitext
from os import makedirs
from qiime.parse import (parse_illumina_line, IlluminaParseError)

# UNTESTED FUNCTIONS FOR WORKING WITH ILLUMINA QUAL LINES 
# def qual_char_to_prob(c):
#     return 1 / (10**(ord(c)-64))
#     
# def qual_line_to_probs(l):
#     return [qual_char_to_prob(c) for c in l]

def get_illumina_qual_chars():
    return '@ABCDEFGHIJKLMNOPQRSTUVWXYZ[\]^_`abcdefghijklmnopqrstuvwxyz{|}~'

def bad_chars_from_threshold(threshold):
    i = -1 * int(log10(threshold))
    return {}.fromkeys(list(get_illumina_qual_chars()[:i]))

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


def illumina_read_description_from_read_data(parsed_read):
    """Create a read description from a parsed illumina read object
    """
    return ':'.join(map(str,[\
     parsed_read['Machine Name'],\
     parsed_read['Channel Number'],\
     parsed_read['Tile Number'],\
     parsed_read['X Position'],\
     parsed_read['Full Y Position Field'][:-2]]))

def parse_illumina_paired_end_read_files(read1_file,read2_file,barcode_length,\
    max_bad_run_length,quality_threshold,min_per_read_length,rev_comp_barcode,\
    barcode_max_N=0,seq_max_N=0):
    """Parses Illumina paired-end read file pair
    """
    
    for read1_line, read2_line in izip(read1_file,read2_file):
        read1 = parse_illumina_line(read1_line,barcode_length,rev_comp_barcode)
        read2 = parse_illumina_line(read2_line,barcode_length,rev_comp_barcode)
        
        read1_desc = illumina_read_description_from_read_data(read1)
        read2_desc = illumina_read_description_from_read_data(read2)
        
        read1_barcode = read1['Barcode']
        read2_barcode = read2['Barcode']
        
        if (read1_barcode.count('N') > barcode_max_N) or \
           (read2_barcode.count('N') > barcode_max_N):
           continue
        
        if read1_desc != read2_desc:
            raise IlluminaParseError, \
              "Error in sequence files, descriptions of"+\
              " corresponding lines are not compatible: %s != %s" %\
              (read1_desc, read2_desc)
        assert read1_barcode == read2_barcode
        
        seq1, qual1 = read_qual_score_filter(\
         read1['Sequence'], read1['Quality Score'],\
         max_bad_run_length, quality_threshold)
         
        if (len(seq1) < min_per_read_length) or (seq1.count('N') > seq_max_N):
            continue
            
        seq2, qual2 = read_qual_score_filter(\
         read2['Sequence'], read2['Quality Score'],\
         max_bad_run_length, quality_threshold)
         
        if (len(seq2) < min_per_read_length) or (seq2.count('N') > seq_max_N):
            continue
            
        seq = seq1 + revComp(seq2)
        qual = qual1 + qual2[::-1]
        
        yield read1_desc, read1_barcode, seq, qual


def process_illumina_paired_end_read_files(
    read1_fp,read2_fp,output_seqs_fp,output_qual_fp,
    barcode_to_sample_id,barcode_length,
    store_unassigned,max_bad_run_length,
    quality_threshold,min_per_read_length,
    rev_comp_barcode,start_seq_id=0):
    """parses Ilimuna paired-end read file
    """
    read1_file = open(read1_fp)
    read2_file = open(read2_fp)
    output_seqs_file = open(output_seqs_fp,'w')
    output_qual_file = open(output_qual_fp,'w')
    
    seq_id = start_seq_id
    
    for seq_desc,barcode,seq,qual in\
      parse_illumina_paired_end_read_files(read1_file,read2_file,barcode_length,\
      max_bad_run_length,quality_threshold,min_per_read_length,rev_comp_barcode):
      
      try:
          sample_id = barcode_to_sample_id[barcode]
      except KeyError:
          if not store_unassigned:
              continue
          else:
              sample_id = 'Unassigned'
              
      fasta_header = '%s_%s %s' % (sample_id,seq_id,seq_desc)
      output_seqs_file.write('>%s\n%s\n' % (fasta_header,seq))
      output_qual_file.write('>%s\n%s\n' % (fasta_header,qual))
      seq_id += 1
      
    output_seqs_file.close()
    output_qual_file.close()
    
    return seq_id
    
def parse_illumina_single_end_read_file(read_file,barcode_length,\
    max_bad_run_length,quality_threshold,min_per_read_length,
    rev_comp,rev_comp_barcode,barcode_max_N=0,seq_max_N=0):
    """Parses Illumina single-end read file
    """
    
    for read_line in read_file:
        read = parse_illumina_line(read_line,barcode_length,rev_comp_barcode)
        
        read_desc = illumina_read_description_from_read_data(read)
        
        read_barcode = read['Barcode']
        
        if read_barcode.count('N') > barcode_max_N:
           continue
        
        seq, qual = read_qual_score_filter(\
         read['Sequence'], read['Quality Score'],\
         max_bad_run_length, quality_threshold)
         
        if (len(seq) < min_per_read_length) or (seq.count('N') > seq_max_N):
            continue
            
        if rev_comp:
            seq = revComp(seq)
            qual = qual[::-1]
        
        yield read_desc, read_barcode, seq, qual
    
def process_illumina_single_end_read_file(read_fp,output_seqs_fp,output_qual_fp,
    barcode_to_sample_id,barcode_length,
    store_unassigned,max_bad_run_length,
    quality_threshold,min_per_read_length, rev_comp, rev_comp_barcode,
    start_seq_id=0):
    """parses Ilimuna single-end read file
    """
    read_file = open(read_fp)
    output_seqs_file = open(output_seqs_fp,'w')
    output_qual_file = open(output_qual_fp,'w')
    
    seq_id = start_seq_id
    
    for seq_desc,barcode,seq,qual in\
      parse_illumina_single_end_read_file(read_file,barcode_length,\
      max_bad_run_length,quality_threshold,min_per_read_length,
      rev_comp,rev_comp_barcode):
        try:
          sample_id = barcode_to_sample_id[barcode]
        except KeyError:
          if not store_unassigned:
              continue
          else:
              sample_id = 'Unassigned'
        fasta_header = '%s_%s %s' % (sample_id,seq_id,seq_desc)
        output_seqs_file.write('>%s\n%s\n' % (fasta_header,seq))
        output_qual_file.write('>%s\n%s\n' % (fasta_header,qual))
        seq_id += 1

    output_seqs_file.close()
    output_qual_file.close()
    
    return seq_id

def mapping_data_to_barcode_map(mapping_data):
    return dict([(d[1].upper(),d[0]) for d in mapping_data])
