#!/usr/bin/env python
# File created on 16 Nov 2011
from __future__ import division

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"

from random import choice


def decorate_numtips(tree):
    """ This function maps the number of children (tips) for each node """

    # iterate over tree and give the score for each node
    for n in tree.postorder(include_self=True):
        if n.istip():
            # if the node is a tip then the number of children is 1
            n.Score = 1
        else:
            # if the node is not a tip then get the number of children
            n.Score = len(n.tips())

    return tree


def decorate_depth(tree):
    """ This function maps the depth of each node on the tree """

    # iterate over the tree and assign depth scores
    for node in tree.levelorder(include_self=True):

        if node.Parent is None:
            # if the node is the root, then the score should start at 0
            node.Score = 0
        else:
            # if the node is not the root, then add it to the parents depth
            node.Score = node.Parent.Score + 1

    return tree


def get_insert_dict(tree, names):
    """ This function returns the nodes labeled as inserted (names) """

    # iterate over the tips and determine if the tip is one of the ones to be
    # removed
    d = {}
    for n in tree.tips():
        if n.Name and n.Name in names:
            if n.Name not in d:
                d[n.Name] = []
            d[n.Name].append(n)

    return d


def drop_duplicate_nodes(tree, inserted_nodes):
    """ remove nodes from tree """

    # iterate over the nodes to be removed and see which one has the highest
    # score which means it will be kept
    for name, node_list in inserted_nodes.items():
        if len(node_list) == 1:
            continue

        # sort the scores so the highest score is kept
        node_list.sort(key=lambda x: x.Score, reverse=True)
        node_to_keep = node_list[0]

        # this will be useful if we want to randomly remove without basing
        # on score
        #node_to_keep = choice(node_list)

        # iterate over the list of nodes to remove and if the node is not
        # designated as a keeper, then remove it from the parent
        for n in node_list:
            if n is node_to_keep:
                continue
            else:
                if n.Parent is not None:
                    n.Parent.remove(n)

    # prune the tree
    tree.prune()

    return tree
