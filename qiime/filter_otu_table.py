#!/usr/bin/env python
#filter_otu_table
"""Filters OTU table according to minimum OTU count and number of samples.

If OTU table has taxonomy assigned, can also use taxonomy to filter.
"""
__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Rob Knight and Jesse Stombaugh"] #remember to add yourself
__license__ = "GPL"
__status__ = "1.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Pre-release"

from qiime.parse import parse_otus
from sys import argv
from string import strip
from numpy import array
from optparse import OptionParser


def strip_quotes(s):
    """splits leading/trailing quotes from string s"""
    if not s or len(s) < 2:
        return s
    if s[0] == s[-1]:
        if s[0] in '"\'':
            s = s[1:-1]
    return s

def split_tax(tax):
    """splits tax string on semicolon and comma"""
    fields = tax.split(';')
    if len(fields) == 1:
        fields = fields[0].split(',')
    return map(strip_quotes, fields)

def process_options(opts):
    
    filepath=opts.otu_fname
    filename=filepath.strip().split('/')[-1]
    filename=filename.split('.')[0]
    
    params={}
    params['otu_file'] = opts.otu_fname
    params['min_otu_count'] = opts.min_count
    params['min_otu_samples'] = opts.min_samples
    
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
    
    filtered_otu_fp = '%s/%s_filtered.txt' % (opts.dir_path,filename)
                                    
    params['filtered_otu_file']=filtered_otu_fp
                                    
    return params

def _filter_table(params):

    filtered_otu_fp=params['filtered_otu_file']
    otu_file=open(params['otu_file'], 'U')
    min_otu_count=params['min_otu_count']
    min_otu_samples=params['min_otu_samples']
    included_taxa=params['included_taxa']
    excluded_taxa=params['excluded_taxa']
    
    filtered_table_path=open(filtered_otu_fp,'w')
    for line in otu_file:
        if line.startswith('#'):
            filtered_table_path.write(line)
        else:
            fields = line.split('\t')
            try:
                vals = array(map(int, fields[1:]), dtype=int)
                taxa = None
            except ValueError:
                vals = array(map(int, fields[1:-1]), dtype=int)
                taxa = set(map(strip, split_tax(fields[-1])))
            if vals.sum() >= min_otu_count and \
                (vals > 0).sum() >= min_otu_samples:
                if not taxa:
                    filtered_table_path.write(line)
                else:
                    if taxa.intersection(included_taxa) and not \
                        taxa.intersection(excluded_taxa):
                        filtered_table_path.write(line)
                    elif not included_taxa and not excluded_taxa:
                        filtered_table_path.write(line)

def make_cmd_parser():
    """Returns the command-line options."""
    parser = OptionParser()
    parser.add_option('-i', '--otu_filename', dest='otu_fname',
        help='otu file name [REQUIRED]')
    parser.add_option('-c', '--min_count', dest='min_count', default=1,
        type=int,
        help='minimum number of sequences to leave OTU in file [default=%default]')
    parser.add_option('-s', '--min_samples', dest='min_samples', default=2,
        type=int,
        help='minimum number of samples to leave OTU in file [default=%default]')
    parser.add_option('-t', '--include_taxonomy', dest='include_taxonomy',
        default='', help='list of taxonomy terms to include [default=%default]')
    parser.add_option('-e', '--exclude_taxonomy', dest='exclude_taxonomy',
        default='', help='list of taxonomy terms to exclude [default=%default]')
    parser.add_option('-o', '--dir_path',\
        help='directory prefix for all analyses [default=%default]',default='./')
    return parser.parse_args()

if __name__ == '__main__':
    from sys import argv, stdout
    opts, args = make_cmd_parser()

    params=process_options(opts)
    
    _filter_table(params)
    