#!/usr/bin/env python
# File created on 08 Jan 2010.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2009, Qiime"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Prototype"


from optparse import OptionParser
from qiime.util import extract_seqs_by_sample_id
from cogent.parse.fasta import MinimalFastaParser

usage_str = """usage: %prog [options] {-i INPUT_FASTA_FP -o OUTPUT_FASTA_FP -s SAMPLE_IDS}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:

 Create the file outseqs.fasta (-o), which will be a subset of inseqs.fasta 
  (-i) containing only the sequences associated with sample ids S2, S3, 
  S4 (-s). (As always, sample IDs are case-sensitive.)
  python extract_seqs_by_sample_id.py -i inseqs.fasta -o outseqs.fasta -s S2,S3,S4
  
 Create the file outseqs.fasta (-o), which will be a subset of inseqs.fasta 
  (-i) containing only the sequences  THAT ARE NOT (-n) associated with sample 
  ids S2, S3, S4 (-s). (As always, sample IDs are case-sensitive.)
  python extract_seqs_by_sample_id.py -i inseqs.fasta -o outseqs.fasta -s S2,S3,S4 -n

"""

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = usage_str
    version = 'Version: %prog ' + __version__
    parser = OptionParser(usage=usage, version=version)

    parser.add_option('-i','--input_fasta_fp',\
        help='the input fasta file [REQUIRED]')
    parser.add_option('-s','--sample_ids',\
        help="comma-separated sample_ids to include in output fasta file"+\
             " (or exclude if -n=True) [REQUIRED]")
    parser.add_option('-o','--output_fasta_fp',\
        help='the output fasta file [REQUIRED]')
    parser.add_option('-n','--negate',action='store_true',\
        help='negate the sample ID list (i.e., output sample '+\
             'ids not passed via -s) [default: %default]')
    
    # Set default values here if they should be other than None
    parser.set_defaults(verbose=False,negate=False)

    opts,args = parser.parse_args()
    required_options = ['input_fasta_fp','output_fasta_fp','sample_ids']
    
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 

    return opts,args

if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    
    sample_ids = opts.sample_ids.split(',')
    negate = opts.negate
    
    try:
        seqs = MinimalFastaParser(open(opts.input_fasta_fp))
    except IOError:
        print "Cannot open %s. Does it exist? Do you have read access?" %\
         opts.input_fasta_fp
        exit(1)
        
    try:
        output_fasta_f = open(opts.output_fasta_fp,'w')
    except IOError:
        print "Cannot open %s. Does path exist? Do you have write access?" %\
         opts.output_fasta_fp
        exit(1)
    
    for r in extract_seqs_by_sample_id(seqs,sample_ids,negate):
        output_fasta_f.write('>%s\n%s\n' % r)
    output_fasta_f.close()