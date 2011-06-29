#!/usr/bin/env python
# File created on 19 Mar 2011
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Release"

import numpy
import random

from qiime.util import make_option

from qiime.parse import parse_otu_table
from qiime.util import parse_command_line_parameters
from qiime.format import format_otu_table
from qiime.sort import natsort

def sim_otu_table(sample_ids, otu_ids, otu_mtx, otu_metadata, tree, 
    num_replicates, dissimilarity):
    """ make n samples related to each sample in an input otu table

    input:
    the constituents of an otu table
     * sample_ids: list
     * otu_ids: list
     * otu_mtx: numpy2d otus as rows
     * otu_metadata: list, either empty or of length len(otu_ids)
    tree: a PhyloNode tree or compatible,
    num_replicates: how many replicates for each input sample
    dissimilarity: how phylogenetically dissimilar each replicate should be 
    relative to the original

    output is a tuple with the constituents of an otu table, possibly missing 
    some otu metadata:
    (res_sam_names, res_otus, res_otu_mtx, res_otu_metadata)
    """

    tree_tips = [tip.Name for tip in tree.tips()]

    # hold sample abundance vector in a dict (sample_dict) temporarily
    # format of sample_dict: otu_id: num_seqs
    sample_dicts = []
    res_sam_names = []
    for i in range(len(sample_ids)):
        sample_vector = otu_mtx[:,i]
        for j in range(num_replicates):
            sample_dict = {}
            for k in range(len(otu_ids)):
                otu_abundance = sample_vector[k]
                if otu_abundance == 0: continue
                new_otu_id = get_new_otu_id(otu_ids[k], tree, dissimilarity)
                # beware, get_new_otu_id may return something we already 
                # have in the table
                if sample_dict.has_key(new_otu_id):
                    sample_dict[new_otu_id] += otu_abundance
                else:
                    sample_dict[new_otu_id] = otu_abundance
            sample_dicts.append(sample_dict)
            res_sam_names.append(sample_ids[i] + '.' + str(j))

    res_otu_mtx, res_otus = combine_sample_dicts(sample_dicts)

    res_otu_metadata = []
    if otu_metadata == None or otu_metadata == []:
        pass #leave res_otu_metadata as []
    else:
        for otu_id in res_otus:
            # if otu was in original table, just copy it's metadata
            try:
                res_otu_metadata.append(otu_metadata[otu_ids.index(otu_id)])
            except ValueError:
            # else just append the empty string
                res_otu_metadata.append('')
    return res_sam_names, res_otus, res_otu_mtx, res_otu_metadata

def get_new_otu_id(old_otu_id, tree, dissim):
    """ simulates an otu switching to related one

    input a tipname, a tree, and a distance to walk up the tree
    ouputs the name of the new, randomly chosen, tree tip
    output tip name may be the same as 
    """
    node = tree.getNodeMatchingName(old_otu_id) #starts at tip
    distance_up_tree = 0
    while (not node.isRoot()) and (distance_up_tree + node.Length < dissim):
        distance_up_tree += node.Length
        node = node.Parent

    # another option is to do 50-50 descending each node, 
    # so we don't bias for overrepresented clades
    if node.isTip():
        return node.Name
    else:
        return random.choice([tip.Name for tip in node.tips()])

def combine_sample_dicts(sample_dicts):
    """ combines a list of sample_dicts into one otu table

    sample dicts is a list of dicts, each one {otu_id:num_seqs}

    output is a tuple:
    (otu_mtx (rows are otus), otu_ids (list))
    * otu_mtx has samples in order of dicts, otus sorted with natsort 
    / human sort
    * otu_mtx will have all otus mentioned as keys in sample_dicts, even if
    they are abundance 0  ({otu_id:0,...})
    such otus will simply be rows of zeros
    """
    all_otu_ids = []
    for s in sample_dicts:
        all_otu_ids.extend(s.keys())

    all_otu_ids = list(set(all_otu_ids))
    all_otu_ids = natsort(all_otu_ids)

    # get index once now, for all samples, instead of all_otu_ids.index()
    indices = {}
    for i in range(len(all_otu_ids)):
        indices[all_otu_ids[i]] = i

    otu_mtx = numpy.zeros((len(all_otu_ids),len(sample_dicts)),int) 
    # otus (rows) by samples (cols)
    for i, sample_dict in enumerate(sample_dicts):
        for otu, abund in sample_dict.items():
            otu_mtx[indices[otu],i] = abund

    return otu_mtx, all_otu_ids