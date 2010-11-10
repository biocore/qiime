#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Rob Knight","Jesse Stombaugh","Dan Knights","Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.2.0-dev"
__maintainer__ = "Tony Walters"
__email__ = "William.A.Walters@colorado.edu"
__status__ = "Development"
 

from qiime.util import parse_command_line_parameters, get_options_lookup, \
                       create_dir
from qiime.filter_otu_table import _filter_table, _filter_table_samples, \
                                   split_tax
from optparse import make_option
from string import strip
from qiime.parse import parse_otu_table
import os

options_lookup = get_options_lookup()

#filter_otu_table.py
script_info={}
script_info['brief_description']="""Filters OTU table by minimum OTU count and number of samples or by taxonomy"""
script_info['script_description']="""After the OTU has been generated, the user may want to filter the table based on the number of samples within each OTU or by the number of sequences per OTU. This step is generally done to reduce the noise within the OTU table and can also reduce the overall size of the table, which is essential when performing analyses on large datasets. Along with filtering based on samples or sequences, the user can include and exclude specific taxon groups."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Examples:""","""To filter the OTU table using the default parameters ("-c 1" and "-s 2"), then write the results to the current working directory, you can use the code as follows:""","""filter_otu_table.py -i otu_table.txt"""))
script_info['script_usage'].append(("""""","""To filter by the number of samples (i.e., 5) within each OTU (keep only OTU's found in at least X number of samples), you can use the code as follows:""","""filter_otu_table.py -i otu_table.txt -s 5"""))
script_info['script_usage'].append(("""""","""To filter by the number of sequences (i.e., 5) within each OTU (keep only OTU's with at least X sequences in the OTU), you can use the code as follows:""","""filter_otu_table.py -i otu_table.txt -c 5"""))
script_info['script_usage'].append(("""""","""To include ("Bacteria") and exclude ("Proteobacteria") certain taxon groups (options -t / -e respectively), you can use the code as follows.  The include and exclude parameters must be used together:""","""filter_otu_table.py -i otu_table.txt -t Bacteria -e Proteobacteria"""))
script_info['script_usage'].append(("""Filter samples by number of sequences""","""A user may want to remove samples that have low sequence coverage. NOTE: this feature is mutually exclusive from the other filtering options, so if you pass this, you will need to perform a subsequent filter to remove by the other options.""","""filter_otu_table.py -i otu_table.txt -p 10"""))
script_info['output_description']="""The result of filter_otu_table.py creates a new OTU table, where the filename uses the input OTU filename and appends "filtered.txt" to the end of the name."""
script_info['required_options'] = [options_lookup['otu_table_as_primary_input']]

script_info['optional_options']=[\
    make_option('-c', '--min_count', default=1, type=int,
        help='retain OTUs with at least this many sequences [default=%default]'),\
    make_option('-s', '--min_samples', default=2, type=int,
        help='retain OTUs found in at least this many samples [default=%default]'),\
    make_option('-t', '--include_taxonomy', default='',
        help='list of taxonomy terms to include [default=%default]'),\
    make_option('-e', '--exclude_taxonomy', default='', 
        help='list of taxonomy terms to exclude [default=%default]'),\
    make_option('-o', '--dir_path', default='./',
        help='directory prefix for all analyses [default=%default]'),
    make_option('-p', '--seqs_per_sample',type=int,
        help='minimum sequences per sample to retain the sample. [default=%default]')
]
script_info['version'] = __version__



def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    # process options (was originally process_options())
    filepath=opts.otu_table_fp
    filename=filepath.strip().split('/')[-1]
    filename=filename.split('.')[0]

    params={}
    params['otu_file'] = opts.otu_table_fp
    params['min_otu_count'] = opts.min_count
    params['min_otu_samples'] = opts.min_samples
    params['seqs_per_sample']=opts.seqs_per_sample
    dir_path=opts.dir_path
    seqs_per_sample=params['seqs_per_sample']
    otu_file=open(params['otu_file'], 'U')

    if not os.path.exists(dir_path):
        create_dir(dir_path,False)
        
    if not dir_path.endswith("/"):
        dir_path=dir_path+"/"
        
        
    filtered_table_path=open(os.path.join(dir_path, \
                                'otu_table_filtered.txt'), 'w')
                                
    if opts.include_taxonomy:
        included_taxa = set(map(strip, split_tax(opts.include_taxonomy)))
    else:
        included_taxa = set()

    if opts.exclude_taxonomy:
        excluded_taxa = set(map(strip, split_tax(opts.exclude_taxonomy)))
    else:
        excluded_taxa=set()

    params['included_taxa']=included_taxa
    params['excluded_taxa']=excluded_taxa

    seq_sample_output=[]
    if seqs_per_sample and params['min_otu_count']==1 and \
            params['min_otu_samples']==2 and opts.include_taxonomy=='' and \
            opts.exclude_taxonomy=='':
        otu_file2=otu_file.readlines()
        filtered_table = _filter_table_samples(otu_file2,seqs_per_sample)
        filtered_table_path.write(filtered_table)
        filtered_table_path.close()
    elif seqs_per_sample and (params['min_otu_count']<>1 or \
            params['min_otu_samples']<>2 or opts.include_taxonomy<>'' or \
            opts.exclude_taxonomy<>''):
        raise ValueError, 'You cannot supply seqs per sample with other filtering options. These features are mutually exclusive.'
    else:
        _filter_table(params,filtered_table_path,otu_file)
        

if __name__ == "__main__":
    main()
