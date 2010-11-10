#!/usr/bin/env python

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Rob Knight", "Catherine Lozupone", "Justin Kuczynski","Julia Goodrich"]
__license__ = "GPL"
__version__ = "1.2.0"
__maintainer__ = "Daniel McDonald"
__email__ = "wasade@gmail.com"
__status__ = "Release"

"""Contains code for summarizing OTU table with taxa in last field.
"""
from collections import defaultdict
from sys import stdout, stderr
from optparse import OptionParser
from string import strip
from numpy import array
from qiime.parse import parse_otu_table, parse_mapping_file

def make_summary(otu_table, level): 
    """Returns taxonomy summary data

    header is a list of:
    [(Taxon),sample1,sample2,...]

    taxonomy_summary is a list of lists of:
    [[(taxon1),count,count,...],[(taxon2),count,count,...]...]
    """
    header = ['Taxon']
    header.extend(otu_table[0]) # sample ids

    counts_by_consensus, sample_map = sum_counts_by_consensus(otu_table, level)

    taxonomy_summary = []
    for consensus, otu_counts in sorted(counts_by_consensus.items()):
        new_row = [(consensus)]
        new_row.extend(otu_counts)
        taxonomy_summary.append(new_row)

    return taxonomy_summary, header

def sum_counts_by_consensus(otu_table, level, missing_name='Other'):
    """Returns a dict keyed by consensus, valued by otu counts

    otu counts are summed together if they have the same consensus

    if the consensus string doesn't reach to level, missing_name is appended on
    until the taxonomy string is of length level
    """
    result = {}
    sample_map = dict([(s,i) for i,s in enumerate(otu_table[0])])

    for counts, consensus in zip(otu_table[2], otu_table[3]):
        n_ranks = len(consensus)
        if n_ranks > level:
            consensus = consensus[:level]
        elif n_ranks < level:
            consensus.extend([missing_name for i in range(level - n_ranks)])
        else:
            # consensus is the correct number of levels
            pass

        consensus = tuple(consensus)
        if consensus in result:
            result[consensus] += counts
        else:
            result[consensus] = counts.copy()

    return result, sample_map

def add_summary_mapping(otu_table, mapping, level): 
    """Returns sample summary of sample counts by taxon
    
    Summary is keyed by sample_id, valued by otu counts for each taxon
    Taxon order is a list of taxons where idx n corresponds to otu count idx n
    """
    counts_by_consensus, sample_map = sum_counts_by_consensus(otu_table, level)
    
    summary = defaultdict(list)
    for row in mapping:
        # grab otu idx if the sample exists, otherwise ignore it
        sample_id = row[0]
        if sample_id not in sample_map:
            continue
        otu_idx = sample_map[sample_id]

        for consensus, counts in sorted(counts_by_consensus.items()):
            summary[sample_id].append(counts[otu_idx])

    taxon_order = sorted(counts_by_consensus.keys())

    return summary, taxon_order

