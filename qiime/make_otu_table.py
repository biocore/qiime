#!/usr/bin/env python
#make_otu_table: makes sample x OTU table
__author__ = "Rob Knight"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Rob Knight", "Justin Kuczynski"] #remember to add yourself
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
from qiime.format import format_otu_table
from qiime.parse import fields_to_dict

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
    parser.add_option('-i', '--input_otu_fname', dest='otu_fname',
        help='name of otu file [Required]')
    parser.add_option('-t', '--taxonomy', dest='taxonomy_fname',
        help='name of taxonomy file', default=None)
    parser.add_option('-o', '--output_fname', dest='output_fname',
        help='name of output file [Default is stdout]')
    options, args = parser.parse_args()
    return options, args

if __name__ == "__main__":
    from sys import argv, exit, stderr, stdout
    options, args = make_cmd_parser()
    if options.output_fname:
        outfile = open(options.output_fname, 'w')
    else:
        outfile = stdout
    if not options.taxonomy_fname:
        otu_to_taxonomy = None
    else:
        res = {}
        infile = open(options.taxonomy_fname,'U')
        for line in infile:
            fields = line.split('\t')
            # typically this looks like: 3 SAM1_32 \t Root,Bacteria,Fi... \t 0.9
            # implying otu 3; sample 1, seq 32 (the representative of otu 3);
            # followed by the taxonomy and confidence
            if not len(fields) == 3:
                continue
            otu = fields[0].split(' ')[0]
            res[otu] = fields[1]
        otu_to_taxonomy = res

    otu_to_seqid = fields_to_dict(open(options.otu_fname, 'U'))

    outfile.write(make_otu_map(otu_to_seqid, otu_to_taxonomy))
