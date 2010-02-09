#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Greg Caporaso", "Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Pre-release"

from qiime.util import parse_command_line_parameters
from optparse import make_option
from os.path import split, splitext
from cogent.parse.fasta import MinimalFastaParser
from cogent.util.misc import revComp
from qiime.adjust_seq_orientation import rc_fasta_file, append_rc

script_description = """ """

script_usage = """ Write the reverse complement of all seqs in seqs.fasta (-i) to
  seqs_rc.fasta (default, change output_fp with -o). Each sequence
  description line will have ' RC' appended to the end of it (default,
  leave sequence description lines untouched by passing -r):
 
 adjust_seq_orientation.py -i seqs.fasta
 """

required_options = [\
 make_option('-i','--input_fasta_fp',\
        help='the input fasta file [REQUIRED]')
]

optional_options = [\
 make_option('-o','--output_fp',\
        help='the input fasta file [default: generated from input_fasta_fp]'),\
 make_option('-r','--retain_seq_id',action='store_true',\
        help='leave seq description lines untouched'+\
        ' [default: append " RC" to seq description lines]')
]

def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)
      
    verbose = opts.verbose
    
    input_fasta_fp = opts.input_fasta_fp
    output_fp = opts.output_fp
    retain_seq_id = opts.retain_seq_id
    
    if retain_seq_id:
        seq_desc_mapper = null_seq_desc_mapper
    else:
        seq_desc_mapper = append_rc
    
    if not output_fp:
        input_file_basename, input_file_ext = \
         splitext(split(input_fasta_fp)[1])
        output_fp = '%s_rc%s' % (input_file_basename,input_file_ext)
        
    rc_fasta_file(input_fasta_fp,output_fp,seq_desc_mapper)


if __name__ == "__main__":
    main()