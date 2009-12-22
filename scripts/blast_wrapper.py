#!/usr/bin/env python
# File created on 20 Dec 2009.
from __future__ import division
from optparse import OptionParser
from cogent.parse.fasta import MinimalFastaParser
from qiime.util import qiime_blast_seqs

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2009, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Prototype"

usage_str = """usage: %prog [options] -i INPUT_FASTA_FP -r REFSEQS_FP

This script is a functionality-limited interface to the 
 qiime.util.qiime_blast_seqs function, primarily useful for testing purposes.
 Once that function has been integrated into qiime as the primary blast interface
 it will move to PyCogent. An expanded version of this command line interface may
 replace the script functionality of cogent.app.blast at that point.

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:
Blast all sequences in inseqs.fasta (-i) against a BLAST db constructed 
 from refseqs.fasta (-r).

 python blast_wrapper.py -i inseqs.fasta -r refseqs.fasta
"""

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = usage_str
    version = 'Version: %prog ' + __version__
    parser = OptionParser(usage=usage, version=version)

    parser.add_option('-i','--input_fasta_fp',\
        help='paths to sequences to blast as a fasta file [REQUIRED]')
    parser.add_option('-r','--refseqs_fp',\
        help='path to blast database as a fasta file [REQUIRED]')
    parser.add_option('-n','--num_seqs_per_blast_run',type='int',\
        help='number of sequences passed to each blast call '+\
        '- useful for very large sequence collections [default: %default]')

    parser.set_defaults(verbose=False,num_seqs_per_blast_run=1000)

    opts,args = parser.parse_args()
    required_options = ['input_fasta_fp','refseqs_fp']
    
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 

    return opts,args

if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    verbose = opts.verbose
    
    blast_results = qiime_blast_seqs(\
     seqs=MinimalFastaParser(open(opts.input_fasta_fp)),\
     refseqs_fp=opts.refseqs_fp,\
     seqs_per_blast_run=opts.num_seqs_per_blast_run)
     
    for query_id, blast_result in blast_results.items():
        first_blast_result = blast_result[0][0]
        print '%s: %s %s %s' % (\
         query_id,
         first_blast_result['SUBJECT ID'],
         first_blast_result['E-VALUE'],
         first_blast_result['% IDENTITY'])
    