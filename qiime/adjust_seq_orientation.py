#!/usr/bin/env python
# File created on 07 Oct 2009.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "1.8.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from os.path import split, splitext
from cogent.parse.fasta import MinimalFastaParser
from cogent import DNA

usage_str = """usage: %prog [options] {-i INPUT_FASTA_FP}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:
 Write the reverse complement of all seqs in seqs.fasta (-i) to
  seqs_rc.fasta (default, change output_fp with -o). Each sequence
  description line will have ' RC' appended to the end of it (default,
  leave sequence description lines untouched by passing -r):
 
 python ~/repo/Qiime/qiime/adjust_seq_orientation.py -i seqs.fasta
"""

def null_seq_desc_mapper(s):
    return s

def append_rc(s):
    return s + ' RC'

def rc_fasta_lines(fasta_lines,seq_desc_mapper=append_rc):
    """ 
    """
    for seq_id, seq in MinimalFastaParser(fasta_lines):
        seq_id = seq_desc_mapper(seq_id)
        seq = DNA.rc(seq.upper())
        yield seq_id, seq
    return
        
def rc_fasta_file(fasta_fp,output_fp,seq_id_mapper=append_rc):
    """
    """
    input_f = open(fasta_fp)
    output_f = open(output_fp,'w')
    
    for s in rc_fasta_lines(input_f,seq_id_mapper):
        output_f.write('>%s\n%s\n' % s)

    input_f.close()
    output_f.close()
