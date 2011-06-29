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

from os.path import split, splitext
from os import makedirs
from qiime.util import make_option
from qiime.split_libraries_illumina import (
    mapping_data_to_barcode_map,
    read_qual_score_filter, bad_chars_from_threshold, IlluminaParseError,
    process_illumina_paired_end_read_files, process_illumina_single_end_read_file, 
    mapping_data_to_barcode_map)
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.parse import parse_mapping_file

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Script for processing raw Illumina Genome Analyzer II data."
script_info['script_description'] = "Script for parsing, library splitting, and quality filtering of raw Illumina Genome Analyzer II data."
script_info['script_usage'] = [\
 ("Parse paired-end read data (-5 and -3 provided), write output to s_1_seqs.fasta","","%prog -5 s_1_1_sequences.fasta -3 s_1_2_sequences.fasta -b barcode_map_6bp.txt"),\
 ("Parse 5' read only (-5 only provided), write output to s_1_5prime_seqs.fasta","","%prog -5 s_1_1_sequences.fasta -b barcode_map_6bp.txt"),\
 ("Parse 3' read only (-3 only provided), write output to  s_1_3prime_seqs.fasta","","%prog -3 s_1_2_sequences.fasta -b barcode_map_6bp.txt"),\
 ("Parse multiple 5' read only files (multiple -5 values provided), write output to s_1_5prime_seqs.fasta, s_2_5primer_seqs.fasta","","%prog -5 s_1_1_sequences.fasta,s_2_1_sequences.fasta -b barcode_map_6bp.txt")
]
script_info['output_description']= ""
script_info['required_options'] = [options_lookup['mapping_fp']]
script_info['optional_options'] = [\
 make_option('-5','--five_prime_read_fp',\
  help='the 5\' read filepath [default: %default]'),\
 make_option('-3','--three_prime_read_fp',\
  help='the 3\' read filepath [default: %default]'),\
 make_option('-o','--output_dir',\
    help='output directory [default: %default]', default='./'),\
 make_option('-u','--store_unassigned',action='store_true',\
    help='store seqs which can\'t be assigned to samples'+\
    ' because of unknown barcodes [default: %default]', default=False),\
 make_option('-q','--quality_threshold',type='float',\
    help='max base call error probability to consider high-quality '+\
    '(probability of base call being error, so values closer to 1 mean that'+\
    ' the base call is more likely to be erroneous) [default: %default]',\
    default=1e-5),\
 make_option('-r','--max_bad_run_length',type='int',\
    help='max number of consecutive low quality base calls allowed '+\
    'before truncating a read [default: 1; the read is trucated at the'+\
    'second low quality call]',default=1),\
 make_option('-p','--min_per_read_length',type='int',\
    help='min number of consecutive high quality base calls to include'+\
    'a read (per single end read) [default: %default]',default=75),\
 make_option('-n','--sequence_max_n',type='int',\
    help='maximum number of N characters allowed in a sequence to retain it -- '
    'this is applied after quality trimming, and is total over combined paired '
    'end reads if applicable [default: %default]',default=0),
 make_option('-s','--start_seq_id',type='int',\
    help='start seq_ids as ascending integers beginning with start_seq_id'+\
    '[default: %default]',default=0),\
 make_option('--rev_comp_barcode',action='store_true',\
    help='reverse compliment barcodes before lookup'+\
    '[default: %default]',default=False),\
 make_option('--barcode_in_header',action='store_true',\
    help='barcode is in header line (rather than beginning of sequence)'+\
    '[default: %default]',default=False)
]
script_info['version'] = __version__



def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
       
    if not (opts.five_prime_read_fp or opts.three_prime_read_fp):
        parser.error(\
            'Must provide five_prime_read_fp and/or three_prime_read_fp.')
            
    try:
        five_prime_read_fps = opts.five_prime_read_fp.split(',')
    except AttributeError:
        five_prime_read_fps = []
    try:
        three_prime_read_fps = opts.three_prime_read_fp.split(',')
    except AttributeError:
        three_prime_read_fps = []
        
    if five_prime_read_fps and three_prime_read_fps:
        assert len(five_prime_read_fps) == len(three_prime_read_fps),\
         "Must pass same number of five_prime_read and three_prime_read files."
    output_dir = opts.output_dir
    mapping_fp = opts.mapping_fp
    max_bad_run_length = opts.max_bad_run_length
    quality_threshold = opts.quality_threshold
    min_per_read_length = opts.min_per_read_length
    store_unassigned= opts.store_unassigned
    seq_max_N = opts.sequence_max_n
    start_seq_id = opts.start_seq_id
    barcode_in_seq = not opts.barcode_in_header
    
    # if barcode_in_seq and three_prime_read_fps:
    #     option_parser.error(\
    #         '--barcode_in_header is only supported option '
    #         'for 5\' single end read data -- contact qiime.help@colorado.edu '
    #         'with questions')
    
    try:
        makedirs(output_dir)
    except OSError:
        pass
    
    mapping_data, mapping_headers, mapping_comments =\
     parse_mapping_file(open(mapping_fp,'U'))
    barcode_to_sample_id = mapping_data_to_barcode_map(mapping_data)
    barcode_length = len(barcode_to_sample_id.keys()[0])
    
    if five_prime_read_fps and three_prime_read_fps:
        next_seq_id = start_seq_id
        for five_prime_read_fp, three_prime_read_fp in\
         zip(five_prime_read_fps,three_prime_read_fps):
            five_prime_read_basename = splitext(split(five_prime_read_fp)[1])[0]
            five_prime_read_basename_fields = five_prime_read_basename.split('_')
            output_basename = '%s_%s' %\
             (five_prime_read_basename_fields[0],five_prime_read_basename_fields[1])
            output_seqs_fp = '%s/%s_seqs.fasta' % (output_dir,output_basename)
            log_fp = '%s/%s_seqs.log' % (output_dir,output_basename)
            output_qual_fp = '%s/%s_qual.txt' % (output_dir,output_basename)
        
            next_seq_id = process_illumina_paired_end_read_files(
             five_prime_read_fp,
             three_prime_read_fp,
             output_seqs_fp,
             output_qual_fp,
             barcode_to_sample_id,
             barcode_length=barcode_length,
             store_unassigned=store_unassigned,
             max_bad_run_length=max_bad_run_length,
             quality_threshold=quality_threshold,
             min_per_read_length=min_per_read_length,
             rev_comp_barcode=opts.rev_comp_barcode,
             seq_max_N=seq_max_N,
             start_seq_id=next_seq_id)
    else:
        if five_prime_read_fps:
            read_fps = five_prime_read_fps
            output_fp_str = '5prime'
            rev_comp = False
        else:
            read_fps = three_prime_read_fps
            output_fp_str = '3prime'
            rev_comp = True
            
        next_seq_id = start_seq_id
        for read_fp in read_fps:
            read_basename = splitext(split(read_fp)[1])[0]
            read_basename_fields = read_basename.split('_')
            output_basename = '%s_%s_%s' %\
             (read_basename_fields[0],read_basename_fields[1],output_fp_str)
            output_seqs_fp = '%s/%s_seqs.fasta' % (output_dir,output_basename)
            log_fp = '%s/%s_seqs.log' % (output_dir,output_basename)
            output_qual_fp = '%s/%s_qual.txt' % (output_dir,output_basename)
        
            next_seq_id = process_illumina_single_end_read_file(\
             read_fp,\
             output_seqs_fp,\
             output_qual_fp,\
             barcode_to_sample_id,\
             barcode_length=barcode_length,\
             store_unassigned=store_unassigned,\
             max_bad_run_length=max_bad_run_length,\
             quality_threshold=quality_threshold,\
             min_per_read_length=min_per_read_length,\
             rev_comp=rev_comp,
             rev_comp_barcode=opts.rev_comp_barcode,
             barcode_in_seq=barcode_in_seq,
             seq_max_N=seq_max_N,
             start_seq_id=next_seq_id)


if __name__ == "__main__":
    main()