#!/usr/bin/env python

"""Takes a tree and cuts off specified children"""

from sys import argv
from cogent.parse.tree import DndParser
from cogent.seqsim.tree import RangeNode

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2010, The QIIME Project" #consider project name
__credits__ = ["Daniel McDonald"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"
__status__ = "Pre-release"

def parse_otu_file(lines):
    """otu\ttaxonid\tthreshold"""
    result = {}
    for l in lines:
        otu, taxonid, threshold = l.strip().split()
        if float(threshold) > 0.01:
            break
        result[int(otu)] = taxonid
    return result

def load_tree(tree_str):
    """Loads tree and assigns ids. Returns tree and mapping: id -> node"""
    tree = DndParser(tree_str, constructor=RangeNode)
    tree.assignIds()
    d = tree.indexByAttr('Id')
    return tree, d

def cutoff_children(id_index, otu_map):
    """Cut off the children of a given node"""
    for otu in otu_map:
        node = id_index[otu]
        node.Name = str(otu)
        node.Children = []

def main():
    tree, id_index = load_tree(open(argv[1]))
    otu_map = parse_otu_file(open(argv[2]))
    cutoff_children(id_index, otu_map)

    new_tree = open(argv[1] + '-otus_at_tips.ntree','w')
    new_tree.write(tree.getNewick(with_distances=True))
    new_tree.close()

if __name__ == '__main__':
    main()  
