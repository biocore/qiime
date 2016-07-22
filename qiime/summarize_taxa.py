#!/usr/bin/env python

__author__ = "Rob Knight"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = [
    "Rob Knight", "Catherine Lozupone", "Justin Kuczynski", "Julia Goodrich",
    "Antonio Gonzalez Pena", "Jose Carlos Clemente Litran"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "wasade@gmail.com"

"""Contains code for summarizing OTU table with taxa in last field.
"""
from collections import defaultdict
from sys import stdout, stderr
from optparse import OptionParser
from string import strip
from numpy import array


def make_summary(otu_table,
                 level,
                 upper_percentage,
                 lower_percentage,
                 md_as_string=False,
                 md_identifier="taxonomy"):
    """Returns taxonomy summary data

    header is a list of:
    [(Taxon),sample1,sample2,...]

    taxonomy_summary is a list of lists of:
    [[(taxon1),count,count,...],[(taxon2),count,count,...]...]
    """
    header = ['Taxon']
    header.extend(otu_table.ids())

    counts_by_consensus, sample_map = sum_counts_by_consensus(otu_table,
                                                              level,
                                                              "Other",
                                                              md_as_string,
                                                              md_identifier)

    total_counts = float(sum([sum(i) for i in counts_by_consensus.values()]))
    taxonomy_summary = []
    for consensus, otu_counts in sorted(counts_by_consensus.items()):
        if lower_percentage is not None and \
                otu_counts.sum() > lower_percentage * total_counts:
            continue
        elif upper_percentage is not None and \
                otu_counts.sum() < upper_percentage * total_counts:
            continue
        new_row = [(consensus)]
        new_row.extend(otu_counts)
        taxonomy_summary.append(new_row)

    # We need level-1 cause biom starts from 0
    collapse_f = lambda id_, md: md[md_identifier][level-1]
    collapsed = otu_table.collapse(collapse_f, norm=False, min_group_size=0,
                                   axis='observation')

    return taxonomy_summary, header


def sum_counts_by_consensus(otu_table,
                            level,
                            missing_name='Other',
                            md_as_string=False,
                            md_identifier='taxonomy'):
    """Returns a dict keyed by consensus, valued by otu counts

    otu counts are summed together if they have the same consensus

    if the consensus string doesn't reach to level, missing_name is appended on
    until the taxonomy string is of length level
    """
    if otu_table.metadata(axis='observation') is None:
        raise ValueError("BIOM table does not contain any "
                         "observation metadata (e.g., taxonomy)."
                         " You can add metadata to it using the "
                         "'biom add-metadata' command.")

    result = {}
    sample_map = otu_table._index()

    # Define a function to process the metadata prior to summarizing - this
    # is more convenient than having to check md_as_string on every iteration
    # in the for loop below
    if md_as_string:
        def process_md(v):
            return v.split(';')
    else:
        def process_md(v):
            return v

    for (otu_val, otu_id, otu_metadata) in otu_table.iter(axis='observation'):
        if md_identifier not in otu_metadata:
            raise KeyError("Metadata category '%s' not in OTU %s. Can't "
                           "continue. Did you pass the correct metadata "
                           "identifier?" % (md_identifier, otu_id))

        consensus = process_md(otu_metadata[md_identifier])
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
            result[consensus] += otu_val
        else:
            result[consensus] = otu_val.copy()

    return result, sample_map


def add_summary_mapping(otu_table,
                        mapping,
                        level,
                        md_as_string=False,
                        md_identifier='taxonomy'):
    """Returns sample summary of sample counts by taxon

    Summary is keyed by sample_id, valued by otu counts for each taxon
    Taxon order is a list of taxons where idx n corresponds to otu count idx n
    """
    counts_by_consensus, sample_map = sum_counts_by_consensus(otu_table,
                                                              level,
                                                              "Other",
                                                              md_as_string,
                                                              md_identifier)

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
