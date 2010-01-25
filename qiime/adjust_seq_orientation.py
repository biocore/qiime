#!/usr/bin/env python
# File created on 07 Oct 2009.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__status__ = "1.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Pre-release"


from optparse import OptionParser
from os.path import split, splitext
from cogent.parse.fasta import MinimalFastaParser
from cogent.util.misc import revComp

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

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = usage_str
    version = 'Version: %prog ' + __version__
    parser = OptionParser(usage=usage, version=version)

    # A binary 'verbose' flag
    parser.add_option('-v','--verbose',action='store_true',\
        dest='verbose',help='Print information during execution -- '+\
        'useful for debugging [default: %default]')

    parser.add_option('-i','--input_fasta_fp',\
        help='the input fasta file [REQUIRED]')

    parser.add_option('-o','--output_fp',\
        help='the input fasta file [default: generated from input_fasta_fp]')

    parser.add_option('-r','--retain_seq_id',action='store_true',\
        help='leave seq description lines untouched'+\
        ' [default: append " RC" to seq description lines]')

    # Set default values here if they should be other than None
    parser.set_defaults(verbose=False,retain_seq_id=False)

    required_options = ['input_fasta_fp']

    opts,args = parser.parse_args()
    
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 

    return opts,args

def null_seq_desc_mapper(s):
    return s

def append_rc(s):
    return s + ' RC'

def rc_fasta_lines(fasta_lines,seq_desc_mapper=append_rc):
    """ 
    """
    for seq_id, seq in MinimalFastaParser(fasta_lines):
        seq_id = seq_desc_mapper(seq_id)
        seq = revComp(seq)
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


if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
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