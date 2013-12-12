#!/usr/bin/env python
# File created on 11 Oct 2011
from __future__ import division

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"

from cogent.core.tree import PhyloNode
from cogent.parse.tree import DndParser
import re

def convert_tree_tips(align_map,tree_fp):
    """ rename the starting tree to correspond to the new phylip names, 
        which are assigned to each sequence """
    
    # flip key value pairs
    tree_tip_to_seq_name={}
    for i in align_map:
        tree_tip_to_seq_name[align_map[i]] = i

    # change the tip labels to phylip labels
    open_tree=open(tree_fp)
    tree=DndParser(open_tree, constructor=PhyloNode)
    for node in tree.tips():
        node.Name = tree_tip_to_seq_name[node.Name]
    
    return tree
    
def write_updated_tree_file(updated_tree_fp,tree):
    """ write the tree """
    
    open_tree_fp=open(updated_tree_fp,'w')
    open_tree_fp.write(tree.getNewick(with_distances=True))
    open_tree_fp.close()

    return

def strip_and_rename_unwanted_labels_from_tree(align_map,tree):
    """ rename tree tips to match the input sequence names """
    
    # iterate over tips and strip unwanted text
    for node in tree.tips():
        removed_query_str=re.sub('QUERY___','',str(node.Name))
        new_node_name=re.sub('___\d+','',str(removed_query_str))
        if new_node_name in align_map:
            node.Name = align_map[new_node_name]
            
    return tree
    