#!/usr/bin/env python
# File created on 07 Jun 2011
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"
 
from os import rename
from cogent.util.misc import safe_md5, create_dir
from qiime.util import parse_command_line_parameters, make_option
from qiime.parse import parse_mapping_file
from qiime.split_libraries_fastq import (process_fastq_single_end_read_file,
                                         BARCODE_DECODER_LOOKUP)
from qiime.split_libraries import check_map
from qiime.split_libraries_illumina import get_illumina_qual_chars
from qiime.golay import get_invalid_golay_barcodes

script_info = {}
script_info['brief_description'] = "This script performs demultiplexing of Fastq sequence data where barcodes and sequences are contained in two separate fastq files (common on Illumina runs)."
script_info['script_description'] = ""
script_info['script_usage'] = [("Demultiplex and quality filter two lanes of Illumina run results and write results to ./sl_out/.","","%prog -i s_7_2_sequence.txt,s_8_2_sequence.txt -b s_7_1_sequence.txt,s_8_1_sequence.txt -o ./sl_out/ -m lane7_map.txt,lane8_map.txt")]
script_info['output_description']= ""
script_info['required_options'] = [\
 # Example required option
 make_option('-i',
             '--sequence_read_fps',
             type="existing_filepaths",
             help='the sequence read fastq files (comma-separated if more than one)'),
 make_option('-b',
             '--barcode_read_fps',
             type="existing_filepaths",
             help='the barcode read fastq files (comma-separated if more than one)'),
 make_option('-o',
             '--output_dir',
             type="new_dirpath",
             help='directory to store output files'),
 make_option('-m',
             '--mapping_fps',
             type="existing_filepaths",
             help='metadata mapping files (comma-separated if more than one)'),
]
script_info['optional_options'] = [
     make_option("--retain_unassigned_reads",
        default=False,
        action='store_true',
        help='retain sequences which don\'t map to a barcode in the '+\
        ' mapping file (sample ID will be "Unassigned") [default: %default]'),
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
     make_option('--rev_comp',action='store_true',\
        help='reverse compliment sequence before writing to output file'+\
        ' (useful for reverse-orientation reads) [default: %default]',
        default=False),\
     make_option('--last_bad_quality_char',type='choice',
        choices = list(get_illumina_qual_chars()),
        help='the last character to be considered low quality '+\
        '(i.e., these character and those before it will be considered '+\
        'low quality base calls) [default: %default]',
        default='B'),\
     make_option('--barcode_type',type='string',
        help='The type of barcode used. This can be an integer, '+\
             'e.g. for length 6 barcodes, or golay_12 for golay '+\
             'error-correcting barcodes. Error correction will '+\
             'only be applied for golay_12 barcodes. [default: %default]',
        default='golay_12'),
     make_option('--max_barcode_errors',
        default=1.5, type=float,
        help='maximum number of errors in barcode [default: %default]'),
     # NEED TO FIX THIS FUNCTIONALITY - CURRENTLY READING THE WRONG FIELD
     # make_option('--filter_bad_illumina_qual_digit',
     #    action='store_true',\
     #    help='filter sequences which are tagged as not passing the Illumina'+\
     #    ' quality filter [default: %default]',
     #    default=False),\
]
script_info['version'] = __version__



def main():
    option_parser, opts, args =\
       parse_command_line_parameters(**script_info)
    
    sequence_read_fps = opts.sequence_read_fps
    barcode_read_fps = opts.barcode_read_fps
    mapping_fps = opts.mapping_fps
    
    retain_unassigned_reads = opts.retain_unassigned_reads
    max_bad_run_length = opts.max_bad_run_length
    last_bad_quality_char = opts.last_bad_quality_char
    min_per_read_length = opts.min_per_read_length
    rev_comp = opts.rev_comp
    rev_comp_barcode = opts.rev_comp_barcode
    seq_max_N = opts.sequence_max_n
    start_seq_id = opts.start_seq_id
    # NEED TO FIX THIS FUNCTIONALITY - CURRENTLY READING THE WRONG FIELD
    filter_bad_illumina_qual_digit = False #opts.filter_bad_illumina_qual_digit
    
    barcode_type = opts.barcode_type
    max_barcode_errors = opts.max_barcode_errors
    
    try:
        barcode_correction_fn = BARCODE_DECODER_LOOKUP[barcode_type]
    except KeyError:
        barcode_correction_fn = None
    
    if len(sequence_read_fps) != len(barcode_read_fps):
        parser.error("Same number of sequence and barcode files must be provided.")
    
    output_dir = opts.output_dir
    create_dir(output_dir)
    
    output_fp_temp = '%s/seqs.fna.incomplete' % output_dir
    output_fp = '%s/seqs.fna' % output_dir
    output_f = open(output_fp_temp,'w')
    log_fp = '%s/split_library_log.txt' % output_dir
    log_f = open(log_fp,'w')
    histogram_fp = '%s/histograms.txt' % output_dir
    histogram_f = open(histogram_fp,'w')
    
    for sequence_read_fp, barcode_read_fp, mapping_fp in\
      zip(sequence_read_fps, barcode_read_fps, mapping_fps):
        mapping_f = open(mapping_fp, 'U')
        h, i, barcode_to_sample_id, warnings, errors, p, a =\
           check_map(mapping_f, disable_primer_check=True)
           
            
        if barcode_type == 'golay_12':
            invalid_golay_barcodes = \
             get_invalid_golay_barcodes(barcode_to_sample_id.keys())
            if len(invalid_golay_barcodes) > 0:
                option_parser.error("Some or all barcodes are not valid golay codes. "+\
                "Do they need to be reverse complimented? If these are not "+\
                "golay barcodes pass --barcode_type 12 to disable barcode "+\
                "error correction. Invalid codes:\n\t%s" % \
                ' '.join(invalid_golay_barcodes))
        
        log_f.write("Input file paths\n")
        log_f.write('Mapping filepath: %s (md5: %s)\n' %\
          (mapping_fp,safe_md5(open(mapping_fp)).hexdigest()))
        log_f.write('Sequence read filepath: %s (md5: %s)\n' %\
          (sequence_read_fp,str(safe_md5(open(sequence_read_fp)).hexdigest())))
        log_f.write('Barcode read filepath: %s (md5: %s)\n\n' %\
          (barcode_read_fp,safe_md5(open(barcode_read_fp)).hexdigest()))
        
        sequence_read_f = open(sequence_read_fp,'U')
        barcode_read_f = open(barcode_read_fp,'U')
        seq_id = start_seq_id
        for fasta_header, sequence, quality, seq_id in \
            process_fastq_single_end_read_file(
               sequence_read_f,
               barcode_read_f,
               barcode_to_sample_id,
               store_unassigned=retain_unassigned_reads,
               max_bad_run_length=max_bad_run_length,
               last_bad_quality_char=last_bad_quality_char,
               min_per_read_length=min_per_read_length,
               rev_comp=rev_comp,
               rev_comp_barcode=rev_comp_barcode,
               seq_max_N=seq_max_N,
               start_seq_id=start_seq_id,
               filter_bad_illumina_qual_digit=\
                filter_bad_illumina_qual_digit,
               log_f=log_f,
               histogram_f=histogram_f,
               barcode_correction_fn=barcode_correction_fn,
               max_barcode_errors=max_barcode_errors):
            output_f.write('>%s\n%s\n' % (fasta_header,sequence))
        start_seq_id = seq_id + 1                                       
        log_f.write('\n---\n\n')
        
    output_f.close()
    rename(output_fp_temp,output_fp)

if __name__ == "__main__":
    main()