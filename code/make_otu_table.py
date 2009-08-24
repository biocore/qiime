#!/usr/bin/env python
#make_otu_table: makes sample x OTU table
__author__ = "Rob Knight"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Rob Knight"] #remember to add yourself
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Prototype"

"""Makes sample x OTU table from OTU map and taxonomy.

Assumes that in the OTU map, the ids are in the format lib_seq, e.g.
M3FclSwb_1023. Will not work if this assumption is not met. Splits on last
underscore only so should be relatively robust to underscore in sample id.
"""
from collections import defaultdict
from string import strip
from numpy import array
from optparse import OptionParser
from cogent.util.misc import flatten, InverseDict
from numpy import zeros
from pipe454.format import format_otu_table
from pipe454.parse import fields_to_dict

def libs_from_seqids(seq_ids, delim='_'):
    """Returns set of libraries."""
    all_libs = set([i.rsplit(delim, 1)[0] for i in seq_ids])
    return all_libs

def seqids_from_otu_to_seqid(otu_to_seqid):
    """Returns set of all seq ids from libs"""
    return set(flatten(otu_to_seqid.values()))

def make_otu_map(otu_to_seqid, otu_to_taxonomy=None, delim='_'):
    """Makes OTU map from otu_to_seqid and otu_to_taxonomy maps."""
    all_seqs = seqids_from_otu_to_seqid(otu_to_seqid)
    try:
        all_otus = map(str, sorted(map(int, otu_to_seqid.keys())))
    except ValueError:
        all_otus = sorted(otu_to_seqid.keys())
    all_libs = sorted(libs_from_seqids(all_seqs))
    table = zeros((len(all_otus), len(all_libs)), int)
    for o in all_otus:
        row_idx = all_otus.index(o)
        row = table[row_idx]
        seqids = otu_to_seqid[o]
        for s in seqids:
            lib = s.rsplit(delim, 1)[0]
            row[all_libs.index(lib)] += 1

    if otu_to_taxonomy:
        taxonomy = [otu_to_taxonomy.get(o, 'None') for o in all_otus]
    else:
        taxonomy=None

    return format_otu_table(all_libs, all_otus, table, taxonomy)

def make_cmd_parser():
    """Returns command-line options"""
    parser = OptionParser()
    parser.add_option('-o', '--otu', dest='otu_fname',
        help='name of otu file')
    parser.add_option('-t', '--taxonomy', dest='taxonomy_fname',
        help='name of taxonomy file', default=None)
    options, args = parser.parse_args()
    return options, args

if __name__ == "__main__":
    from sys import argv, exit, stderr
    options, args = make_cmd_parser()
    if not options.taxonomy_fname:
        otu_to_taxonomy = None
    else:
        res = {}
        infile = open(options.taxonomy_fname,'U')
        for line in infile:
            fields = line.split('\t')
            if not len(fields) == 3:
                continue
            res[fields[0]] = fields[1]
        otu_to_taxonomy = res

    otu_to_seqid = fields_to_dict(open(options.otu_fname, 'U'))

    print make_otu_map(otu_to_seqid, otu_to_taxonomy)
