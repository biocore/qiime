#!/usr/bin/env python
#filter_otu_table
"""Filters OTU table according to minimum OTU count and number of samples.

If OTU table has taxonomy assigned, can also use taxonomy to filter.
"""
__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Rob Knight and Jesse Stombaugh"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Pre-release"

from qiime.parse import parse_otus
from string import strip
from numpy import array
import os


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

def _filter_table(params):

    dir_path=params['dir_path']
    otu_file=open(params['otu_file'], 'U')
    min_otu_count=params['min_otu_count']
    min_otu_samples=params['min_otu_samples']
    included_taxa=params['included_taxa']
    excluded_taxa=params['excluded_taxa']
    
    
    if not os.path.exists(dir_path):
        os.mkdir(dir_path)
        
    if not dir_path.endswith("/"):
        dir_path=dir_path+"/"
        
    filtered_table_path=open(dir_path+'otu_table_filtered.txt','w')
    
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
