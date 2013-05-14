#!/usr/bin/env python
"""Wrapper for Tax2Tree version 1.0

Includes wrapper for Tax2Tree.
"""

__author__ = "Kyle Patnode"
__copyright__ = "Copyright 2012, The QIIME Project"
__credits__ = ["Kyle Patnode", "Daniel McDonald", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Kyle Patnode"
__email__ = "kpatnode1@gmail.com"
__status__ = "Development"

from t2t import nlevel
from cogent.parse.fasta import MinimalFastaParser

def clean_output(assigned_constrings, seq_names):
    """Extracts seq_names consensus strings from assigned_constrings"""
    results = {}
    assigned_constrings = dict([x.split('\t') for x in assigned_constrings])

    for name in seq_names:
        results[name] = (assigned_constrings[name.split()[0]], 0)

    return results

def expand_constrings(sample_id_map):
    """Adds spaces to the constrings to fix incompatibility with nlevel."""
    new_id_map = dict([x.split('\t') for x in sample_id_map])

    for seq_id, constring in new_id_map.iteritems():
        #nlevel looks for '; ' when separating constrings.
        #Unfortunately, not all id maps(green genes) have a space so this fixes that.
        if len(constring.split('; ')) == 1:
            new_id_map[seq_id] = '; '.join(new_id_map[seq_id].split(';'))

    return new_id_map

def generate_constrings(tree, tipname_map, verbose=False):
    """Assigns taxonomy to unidentified sequences in tree.

    Returns all sequence IDs on tree."""
    counts = nlevel.collect_names_at_ranks_counts(tree)
    min_count = 2
    nlevel.decorate_ntips(tree)
    nlevel.decorate_name_relative_freqs(tree, counts, min_count)
    nlevel.set_ranksafe(tree)
    nlevel.pick_names(tree)
    nlevel.name_node_score_fold(tree)

    if verbose:
        print "Tree score: ", nlevel.score_tree(tree)

    nlevel.set_preliminary_name_and_rank(tree)
    contree, contree_lookup = nlevel.make_consensus_tree(tipname_map.values())
    nlevel.backfill_names_gap(tree, contree_lookup)
    nlevel.commonname_promotion(tree)
    nlevel.make_names_unique(tree, append_suffix=False)

    constrings = nlevel.pull_consensus_strings(tree)

    return constrings

def prep_consensus(sample_id_map, seq_names):
    """Adds representative sequences to the consensus map and marks them to be identified."""
    new_id_map = expand_constrings(sample_id_map)

    for duplicate in set(seq_names).intersection(set(new_id_map.keys())):
        #Resolve duplicate names
        new_id_map[duplicate + '_R'] = new_id_map.pop(duplicate)

    unclassified_map = dict()
    #Creates a map of placeholders for the sequences to be identified.
    for name in seq_names:
        unclassified_map[name] = 'Unclassified'

    new_id_map.update(unclassified_map)

    result = []
    for seq_id, constring in new_id_map.iteritems():
        result.append(seq_id+'\t'+constring)

    return result
