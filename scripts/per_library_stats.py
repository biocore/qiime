#!/usr/bin/env python
# File created on 09 Nov 2009.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Pre-release"

from optparse import OptionParser
from qiime.util import compute_seqs_per_library_stats

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

if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    verbose = opts.verbose
    
    min_counts, max_counts, median_counts, mean_counts, counts_per_sample =\
     compute_seqs_per_library_stats(open(opts.input_otu_table,'U'))
    
    print 'Seqs/sample summary:' 
    print ' Min: %s' % str(min_counts)
    print ' Max: %s' % str(max_counts)
    print ' Median: %s' % str(median_counts)
    print ' Mean: %s' % str(mean_counts)
    print ''
    print 'Seqs/sample detail:'
    for k,v in counts_per_sample.items():
        print ' %s: %s' % (k,str(v))
    
    
    
    