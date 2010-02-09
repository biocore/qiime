#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight and Micah Hamady"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ =  ["Rob Knight", "Micah Hamady", "Greg Caporaso", "Kyle Bittinger","Jesse Stombaugh","William Walters", "Jens Reeder"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "William Walters"
__email__ = "william.a.walters@colorado.edu"
__status__ = "Pre-release"
 

from qiime.util import parse_command_line_parameters
from optparse import make_option

from qiime.split_libraries import preprocess

script_description = """Split libraries according to barcodes specified in mapping file."""


script_usage = """
Process sequences from two files (a.fna, b.fna) using the sample-to-
barcode associations in samples.txt:

   %prog -f a.fna,b.fna -m samples.txt

Process sequences using the quality scores (a.qual, b.qual) enabling a more robust filtering:

   %prog -f a.fna,b.fna -q a.qual,b.qual -m samples.txt
 """

required_options = [\
    make_option('-m', '--map', dest='map_fname', 
                help='name of mapping file. NOTE: Must contain a header'+\
                    ' line indicating SampleID in the first column and'+\
                    ' BarcodeSequence in the second,'+\
                    ' LinkerPrimerSequence in the third.'),
    make_option('-f', '--fasta', dest='fasta_fnames', 
                help='names of fasta files, comma-delimited')
    ]

optional_options = [\
    make_option('-q', '--qual', dest='qual_fnames', 
        help='names of qual files, comma-delimited [default: %default]'),

    make_option('-l', '--min-seq-length', dest='min_seq_len',
        type=int, default=200,
        help='minimum sequence length, in nucleotides [default: %default]'),

    make_option('-L', '--max-seq-length', dest='max_seq_len',
        type=int, default=1000,
        help='maximum sequence length, in nucleotides [default: %default]'),

    make_option('-t', '--trim-seq-length', dest='trim_seq_len',
        action='store_true',
        help='calculate sequence lengths after trimming primers and barcodes'+\
         ' [default: %default]', default=False),

    make_option('-s', '--min-qual-score', type=int, default=25,
        help='min average qual score allowed in read [default: %default]'),

    make_option('-k', '--keep-primer', action='store_true',
        help='do not remove primer from sequences', default=False),

    make_option('-B', '--keep-barcode', action='store_true',
        help='do not remove barcode from sequences', default=False),

    make_option('-a', '--max-ambig', type=int, default=0,
        help='maximum number of ambiguous bases [default: %default]'),

    make_option('-H', '--max-homopolymer', type=int, default=6,
        help='maximum length of homopolymer run [default: %default]'),

    make_option('-M', '--max-primer-mismatch', dest='max_primer_mm',
        type=int, default=0,
        help='maximum number of primer mismatches [default: %default]'),

    make_option('-b', '--barcode-type', default='golay_12', 
        help=\
        'barcode type, e.g. 4 or hamming_8 or golay_12 [default: %default]'),

    make_option('-o', '--dir-prefix', default='.',
        help='directory prefix for output files [default: %default]'),

    make_option('-e', '--max-barcode-errors', dest='max_bc_errors',
        default=1.5, type=float,
        help='maximum number of errors in barcode [default: %default]'),

    make_option('-n', '--start-numbering-at', dest='start_index',
        default=1, type=int,
        help='seq id to use for the first sequence [default: %default]'),

    make_option('-r', '--remove_unassigned', default=False,
        action='store_true', help='remove sequences which are Unassigned from \
            output [default: %default]'),

    make_option('-c', '--disable_bc_correction', default=False,
        action='store_true', help='Disable attempts to find nearest '+\
        'corrected barcode.  Can improve performance. [default: %default]'),

    make_option('-w', '--qual_score_window', dest="qual_score_window",
                type=int, default=0,
        action='store', help='Enable sliding window test of quality '+\
        'scores.  If the average score of a continuous set of w nucleotides '+\
        'falls below the threshold (see -s for default), the sequence is '+\
        'discarded. A good value would be 50. 0 (zero) means no filtering. '+\
        'Must pass a .qual file (see -q parameter) if this '+\
        'functionality is enabled. [default: %default]')
]

def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)
    
    if opts.qual_score_window and not opts.qual_fnames:
        option_parser.error('To enable sliding window quality test (-w), .qual '+\
         'files must be included.')
  
    mapping_file = opts.map_fname
    fasta_files = set(opts.fasta_fnames.split(','))
    if opts.qual_fnames:
        qual_files = set(opts.qual_fnames.split(','))
    else:
        qual_files = set()

    for q in qual_files:
        if not q.endswith('qual'):
            stderr.write(
            "Qual file does not end with .qual: is it really a qual file?\n%s\n" 
            % q)

    for f in fasta_files:
        if not (f.endswith('fasta') or f.endswith('fna')):
            stderr.write(
            "Fasta file does not end with .fna: is it really a seq file?\n%s\n" 
            % f)
    
    preprocess(fasta_files, qual_files, mapping_file,
               barcode_type=opts.barcode_type,
               starting_ix = opts.start_index,
               min_seq_len = opts.min_seq_len,
               max_seq_len = opts.max_seq_len, 
               min_qual_score=opts.min_qual_score,
               keep_barcode=opts.keep_barcode,
               keep_primer=opts.keep_primer,
               max_ambig=opts.max_ambig,
               max_primer_mm=opts.max_primer_mm,
               trim_seq_len=opts.trim_seq_len,
               dir_prefix=opts.dir_prefix,
               max_bc_errors = opts.max_bc_errors,
               max_homopolymer = opts.max_homopolymer,
               remove_unassigned = opts.remove_unassigned,
               attempt_bc_correction = not opts.disable_bc_correction,
               qual_score_window = opts.qual_score_window)
 
if __name__ == "__main__":
    main()
