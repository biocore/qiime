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

from optparse import make_option
from qiime.util import (compute_seqs_per_library_stats, 
    parse_command_line_parameters)

script_description = """Given an otu table compute and print the (min, max, median, mean) number of seqs per library."""

script_usage = """Example usage: calculate stats on OTU table specified by otu_Table.txt

per_library_stats.py -i otu_table.txt
"""

required_options = [make_option('-i','--input_otu_table',\
         help='the input otu table')
]

optional_options = []

def main():
    option_parser, opts,args = parse_command_line_parameters(
        script_description=script_description,
        script_usage=script_usage,
        version=__version__,
        required_options=required_options,
        optional_options=optional_options)
    
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
    
if __name__ == "__main__":
    main()
