#!/usr/bin/env python
# File created on 09 Nov 2009.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Daniel McDonald","Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Release"

from qiime.util import make_option
from numpy import std
from qiime.util import (compute_seqs_per_library_stats, 
    parse_command_line_parameters, get_options_lookup, 
    median_absolute_deviation, guess_even_sampling_depth)
from qiime.format import format_mapping_file
from qiime.parse import parse_mapping_file

options_lookup = get_options_lookup()

#per_library_stats.py
script_info={}
script_info['brief_description']="""Calculate per library statistics"""
script_info['script_description']="""Given an otu table, compute and print the (min, max, median, mean) number of seqs per library."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""Calculate statistics on an OTU table (otu_table.txt)""","""per_library_stats.py -i otu_table.txt"""))
script_info['script_usage'].append(("""Example appending results to mapping file:""","""Calculate statistics on an OTU table (otu_table.txt)""","""per_library_stats.py -i otu_table.txt -m old_map.txt -o new_map.txt"""))

script_info['output_description']="""The resulting statistics are written to stdout. If -m is passed, a new mapping file is written to the path specified by -o, in addition to the statistics written to stdout"""
script_info['required_options']=[options_lookup['otu_table_as_primary_input']]
script_info['optional_options']=[
make_option('-m','--mapfile',help='a mapping file. If included, this script will modify the mapping file to include sequences per sample (library) information, and write the modified mapping file to the path specified by -o. The sequences (individuals) per sample is presented in a new column entitled "NumIndividuals", and samples present in the mapping file but not the otu table have the value "na" in this column. Note also that the location of comments is not preserved in the new mapping file'),

make_option('-o','--outputfile',help='the output filepath where the modified mapping file will be written')
]
script_info['version'] = __version__

def main():
    option_parser, opts,args = parse_command_line_parameters(**script_info)

    min_counts, max_counts, median_counts, mean_counts, counts_per_sample =\
     compute_seqs_per_library_stats(open(opts.otu_table_fp,'U'))

    counts_per_sample_values = counts_per_sample.values()
    med_abs_dev = median_absolute_deviation(counts_per_sample_values)[0]
    even_sampling_depth = guess_even_sampling_depth(counts_per_sample_values)

    print 'Num samples: %s\n' % str(len(counts_per_sample))

    print 'Seqs/sample summary:' 
    print ' Min: %s' % str(min_counts)
    print ' Max: %s' % str(max_counts)
    print ' Median: %s' % str(median_counts)
    print ' Mean: %s' % str(mean_counts)
    print ' Std. dev.: %s' % (str(std(counts_per_sample_values)))
    print ' Median Absolute Deviation: %s' % str(med_abs_dev)
    print ' Default even sampling depth in\n  core_qiime_analyses.py (just a suggestion): %s' %\
     str(even_sampling_depth)
    print ''
    print 'Seqs/sample detail:'
    sorted_counts_per_sample = [(v,k) for k,v in counts_per_sample.items()]
    sorted_counts_per_sample.sort()
    total_count = 0
    for v,k in sorted_counts_per_sample:
        total_count += v
        print ' %s: %s' % (k,str(v))
    print '\nTotal observations (sequences): %d' % total_count

    if opts.mapfile:
        if not opts.outputfile:
            raise RuntimeError('input mapping file supplied, but no path to'+\
             ' output file')
        f = open(opts.mapfile,'U')
        mapping_lines, headers, comments = parse_mapping_file(f)
        f.close()
        if len(headers)==1:
            endoffset = 0 # if we only have the sample id, this data -> last col
        else:
            endoffset = 1 # usually make this data the penultimate column.
        headers.insert(len(headers)-endoffset,'NumIndividuals')
        for map_line in mapping_lines:
            sample_id = map_line
            try:
                depth = str(counts_per_sample[map_line[0]])
            except KeyError:
                depth = 'na'
            map_line.insert(len(map_line)-endoffset,depth)

        new_map_str = format_mapping_file(headers, mapping_lines, comments)
        f = open(opts.outputfile, 'w')
        f.write(new_map_str)
        f.close()

if __name__ == "__main__":
    main()
