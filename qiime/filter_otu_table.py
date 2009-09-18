#!/usr/bin/env python
#filter_otu_table
"""Filters OTU table according to minimum OTU count and number of samples.

If OTU table has taxonomy assigned, can also use taxonomy to filter.
"""
__author__ = "Rob Knight"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Rob Knight"] #remember to add yourself
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Prototype"

from qiime.parse import parse_otus
from sys import argv, stdout
from string import strip
from numpy import array
from optparse import OptionParser

def make_cmd_parser():
    """Returns the command-line options."""
    parser = OptionParser()
    parser.add_option('-c', '--min_count', dest='min_count', default=5,
        type=int,
        help='minimum number of sequences to leave OTU in file.')
    parser.add_option('-s', '--min_samples', dest='min_samples', default=2,
        type=int,
        help='minimum number of samples to leave OTU in file.')
    parser.add_option('-i', '--include_taxonomy', dest='include_taxonomy',
        default='', help='list of taxonomy terms to include')
    parser.add_option('-e', '--exclude_taxonomy', dest='exclude_taxonomy',
        default='', help='list of taxonomy terms to exclude')
    parser.add_option('-o', '--otu_filename', dest='otu_fname',
        help='otu file name')
    return parser.parse_args()

def split_tax(tax):
    """splits tax string on semicolon and comma"""
    fields = tax.split(';')
    if len(fields) == 1:
        fields = fields[0].split(',')
    return fields
    

if __name__ == '__main__':
    from sys import argv, stdout
    opts, args = make_cmd_parser()
    min_otu_count = opts.min_count
    min_otu_samples = opts.min_samples
    if opts.include_taxonomy:
        tax = opts.include_taxonomy.split(';')
        if len(tax) == 1:
            included_taxa = set(map(strip, split_tax(opts.include_taxonomy)))
    else:
        included_taxa = set()
    if opts.exclude_taxonomy:
        excluded_taxa = set(map(strip, split_tax(opts.exclude_taxonomy)))
    else:
        excluded_taxa=set()
    
    otu_file = open(opts.otu_fname, 'U')

    for line in otu_file:
        if line.startswith('#'):
            stdout.write(line)
        else:
            fields = line.split('\t')
            try:
                vals = array(map(int, fields[1:]), dtype=int)
                taxa = None
            except ValueError:
                vals = array(map(int, fields[1:-1]), dtype=int)
                taxa = set(map(strip, split_tax(fields[-1])))
            if vals.sum() >= min_otu_count and \
                (vals > 0).sum() > min_otu_samples:
                if not taxa:
                    stdout.write(line)
                else:
                    if taxa.intersection(included_taxa) and not \
                        taxa.intersection(excluded_taxa):
                        stdout.write(line)
