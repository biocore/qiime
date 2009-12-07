#!/usr/bin/env python
# File created on 09 Nov 2009.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2009, the Qiime project"
__credits__ = ["CONTRIBUTORS"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "<AUTHOR NAME>"
__email__ = "<AUTHOR EMAIL>"
__status__ = "Prototype"

from numpy import min, max, median, mean
from cogent.parse.fasta import MinimalFastaParser
from optparse import OptionParser
from qiime.parse import parse_otus

usage_str = """usage: %prog [options] {-i INPUT_OTU_TABLE}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Given an otu table compute and print the (min, max, median, mean) 
 number of seqs per library.

Example usage:

 python Qiime/scripts/per_library_stats.py -i otu_table.txt


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

    parser.add_option('-i','--input_otu_table',\
         help='the input otu table [REQUIRED]')

    # Set default values here if they should be other than None
    parser.set_defaults(verbose=False)

    opts,args = parser.parse_args()
    required_options = ['input_otu_table']
    
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 

    return opts,args

# def compute_stats(fasta_f,lib_seq_delimited='_'):
#     """ """
#     d = {}
#     for seq_id, seq in MinimalFastaParser(fasta_f):
#         lib_id = seq_id.split(lib_seq_delimited)[0]
#         try:
#             d[lib_id] +=1
#         except KeyError:
#             d[lib_id] = 1
#             
#     counts = d.values()
#     return min(counts), max(counts), median(counts), mean(counts)

def compute_stats(otu_f):
    counts = []
    otu_table = parse_otus(otu_f)[2]
    for i in range(otu_table.shape[1]):
        counts.append(sum(otu_table[:,i]))
        
    return min(counts), max(counts), median(counts), mean(counts)

if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    verbose = opts.verbose
    
    print compute_stats(open(opts.input_otu_table,'U'))
    
    