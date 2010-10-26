#!/usr/bin/env python
# File created on 09 Nov 2009.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso", "Daniel McDonald","Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from optparse import make_option
from qiime.util import (compute_seqs_per_library_stats, 
    parse_command_line_parameters, get_options_lookup)

options_lookup = get_options_lookup()

#per_library_stats.py
script_info={}
script_info['brief_description']="""Calculate per library statistics"""
script_info['script_description']="""Given an otu table, compute and print the (min, max, median, mean) number of seqs per library."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""Calculate statistics on an OTU table (otu_table.txt)""","""per_library_stats.py -i otu_table.txt"""))
script_info['output_description']="""The resulting statistics are written to stdout."""
script_info['required_options']=[options_lookup['otu_table_as_primary_input']]
script_info['optional_options']=[]
script_info['version'] = __version__

def main():
    option_parser, opts,args = parse_command_line_parameters(**script_info)
    
    min_counts, max_counts, median_counts, mean_counts, counts_per_sample =\
     compute_seqs_per_library_stats(open(opts.otu_table_fp,'U'))
    
    print 'Num samples: %s\n' % str(len(counts_per_sample))
    
    print 'Seqs/sample summary:' 
    print ' Min: %s' % str(min_counts)
    print ' Max: %s' % str(max_counts)
    print ' Median: %s' % str(median_counts)
    print ' Mean: %s' % str(mean_counts)
    print ''
    print 'Seqs/sample detail:'
    sorted_counts_per_sample = [(v,k) for k,v in counts_per_sample.items()]
    sorted_counts_per_sample.sort()
    for v,k in sorted_counts_per_sample:
        print ' %s: %s' % (k,str(v))
    
if __name__ == "__main__":
    main()
