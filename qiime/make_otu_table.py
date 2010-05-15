#!/usr/bin/env python
#make_otu_table: makes sample x OTU table
__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project" 
__credits__ = ["Rob Knight", "Justin Kuczynski"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Development"

"""Makes sample x OTU table from OTU map and taxonomy.

Assumes that in the OTU map, the ids are in the format lib_seq, e.g.
M3FclSwb_1023. Will not work if this assumption is not met. Splits on last
underscore only so should be relatively robust to underscore in sample id.
"""

from collections import defaultdict
from string import strip
from numpy import array
from cogent.util.misc import flatten, InverseDict
from numpy import zeros
from qiime.format import format_otu_table


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
    try:
        table = zeros((len(all_otus), len(all_libs)), int)
    except MemoryError, e:
        stderr.write('memory error, check format of input otu file\n')
        stderr.write('are there really %s otus and %s samples?\n' %
            (len(all_otus), len(all_libs)))
        stderr.write('traceback follows:\n')
        raise(e)
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
