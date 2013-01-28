#!/usr/bin/env python
# File created on 09 Nov 2009.
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Daniel McDonald","Justin Kuczynski",
               "Jose Carlos Clemente Litran", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Development"

from numpy import std
from qiime.util import (compute_seqs_per_library_stats, 
    parse_command_line_parameters, get_options_lookup, 
    median_absolute_deviation, guess_even_sampling_depth, 
    qiime_open, make_option)
from qiime.format import format_mapping_file
from qiime.parse import parse_mapping_file
from biom.parse import parse_biom_table

options_lookup = get_options_lookup()

#per_library_stats.py
script_info={}
script_info['brief_description']="""Calculate statistics on data in a BIOM-formatted OTU table"""
script_info['script_description']="""Given an otu table in BIOM format, compute and print the (min, max, median, mean) number of seqs per library."""
script_info['script_usage']=[]

script_info['script_usage'].append(("""Write to standard out:""","""Calculate statistics on an OTU table""","""%prog -i otu_table.biom"""))

script_info['script_usage'].append(("""Write to standard out and edit mapping file:""","""Calculate statistics on an OTU table and add sequence/sample count data to mapping file.""","""%prog -i otu_table.biom -m Fasting_Map.txt -o map.txt"""))

script_info['output_description']="""The resulting statistics are written to stdout. If -m is passed, a new mapping file is written to the path specified by -o, in addition to the statistics written to stdout"""
script_info['required_options']=[options_lookup['otu_table_as_primary_input']]
script_info['optional_options']=[
make_option('-m','--mapping_fp',type='existing_filepath',help='a mapping file. If included, this script will modify the mapping file to include sequences per sample (library) information, and write the modified mapping file to the path specified by -o. The sequences per sample is presented in a new column entitled "SequenceCount", and samples present in the mapping file but not the otu table have the value "na" in this column. Note also that the location of comments is not preserved in the new mapping file'),

make_option('-o','--output_mapping_fp',help='the output filepath where the modified mapping file will be written', type='new_filepath'),
make_option('--num_otus',action='store_true',help='Counts are presented as number of observed OTUs per sample, rather than counts of sequences per sample [default: %default]',default=False)
]
script_info['version'] = __version__

def main():
    option_parser, opts,args = parse_command_line_parameters(**script_info)
    otu_table_fp = opts.otu_table_fp
    otu_table = parse_biom_table(qiime_open(otu_table_fp))
    min_counts, max_counts, median_counts, mean_counts, counts_per_sample =\
     compute_seqs_per_library_stats(otu_table, opts.num_otus)
    num_otus = len(otu_table.ObservationIds)
    
    counts_per_sample_values = counts_per_sample.values()
    med_abs_dev = median_absolute_deviation(counts_per_sample_values)[0]
    even_sampling_depth = guess_even_sampling_depth(counts_per_sample_values)
    
    try:
        sample_md_keys = otu_table.SampleMetadata[0].keys()
    except TypeError:
        sample_md_keys = ["None provided"]
    try:
        observation_md_keys = otu_table.ObservationMetadata[0].keys()
    except TypeError:
        observation_md_keys = ["None provided"]
    
    num_samples = len(counts_per_sample)
    print 'Num samples: %s' % str(num_samples)
    print 'Num otus: %s' % str(num_otus)
    if not opts.num_otus:
        num_observations = sum(counts_per_sample_values)
        print 'Num observations (sequences): %s' % str(num_observations)
        print 'Table density (fraction of non-zero values): %1.4f' % \
              otu_table.getTableDensity()
    print

    if opts.num_otus:
        print 'OTUs/sample summary:'
    else:
        print 'Seqs/sample summary:' 
    print ' Min: %s' % str(min_counts)
    print ' Max: %s' % str(max_counts)
    print ' Median: %s' % str(median_counts)
    print ' Mean: %s' % str(mean_counts)
    print ' Std. dev.: %s' % (str(std(counts_per_sample_values)))
    print ' Median Absolute Deviation: %s' % str(med_abs_dev)
    print ' Default even sampling depth in\n  core_qiime_analyses.py (just a suggestion): %s' %\
     str(even_sampling_depth)
    print ' Sample Metadata Categories: %s' % '; '.join(sample_md_keys)
    print ' Observation Metadata Categories: %s' % '; '.join(observation_md_keys)
     
    print ''
    if opts.num_otus:
        print 'OTUs/sample detail:'
    else:
        print 'Seqs/sample detail:'
    sorted_counts_per_sample = [(v,k) for k,v in counts_per_sample.items()]
    sorted_counts_per_sample.sort()
    total_count = 0
    for v,k in sorted_counts_per_sample:
        total_count += v
        print ' %s: %s' % (k,str(v))

    if opts.mapping_fp:
        if not opts.output_mapping_fp:
            raise RuntimeError('input mapping file supplied, but no path to'+\
             ' output file')
        f = open(opts.mapping_fp,'U')
        mapping_lines, headers, comments = parse_mapping_file(f)
        f.close()
        if len(headers)==1:
            endoffset = 0 # if we only have the sample id, this data -> last col
        else:
            endoffset = 1 # usually make this data the penultimate column.
        headers.insert(len(headers)-endoffset,'SequenceCount')
        for map_line in mapping_lines:
            sample_id = map_line
            try:
                depth = str(counts_per_sample[map_line[0]])
            except KeyError:
                depth = 'na'
            map_line.insert(len(map_line)-endoffset,depth)

        new_map_str = format_mapping_file(headers, mapping_lines, comments)
        f = open(opts.output_mapping_fp, 'w')
        f.write(new_map_str)
        f.close()

if __name__ == "__main__":
    main()
