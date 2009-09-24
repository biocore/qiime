#!/usr/bin/env python
from __future__ import division
from optparse import OptionParser
from cogent.core.tree import PhyloNode
from cogent.parse.tree import DndParser
import os.path

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Justin Kuczynski"] #remember to add yourself
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Prototype"

"""compares jackknife/bootstrapped trees with a master, and outputs support for
nodes.  main() is called when run from command line, primary function of this
module is bootstrap_support()

command line usage example:
python tree_compare.py -m TEST/sample_cluster.tre -s TEST/rare_unifrac_upgma -o TEST/unifrac_jackknife
"""


def main(options):
    tree_file_names = os.listdir(options.support_path)
    if not os.path.exists(options.output_path):
        os.mkdir(options.output_path)
    support_trees = []
    for fname in tree_file_names:
        f = open(os.path.join(options.support_path, fname))
        tree = DndParser(f, PhyloNode)
        support_trees.append(tree)
        f.close()
    master_tree = DndParser(open(options.master_tree_path), PhyloNode)

    # get support of each node in master
    new_master, bootstraps = bootstrap_support(master_tree, support_trees)

    # write out modified master tree with internal nodes, and bootstrap support
    # values to 2 files
    fname = os.path.join(options.output_path, "master_tree.tre")
    f = open(fname, 'w')
    f.write(new_master.getNewick(with_distances=True))
    f.close()
    
    f = open(os.path.join(options.output_path, 'jackknife_support.txt'), 'w')
    for key, val in bootstraps.items():
        f.write("\t".join([str(key),str(val)]) + "\n")
    f.close()
    


def bootstrap_support(master_tree, trees):
    """ typically this calculates bootstrap support of each node in master_tree
    a tree supports a given master_tree_node if there exists a node in tree
    where node.tips == master_tree_node.tips (by name)
    not specific to bootstrap, does node support for trees generated in any
    manner

    bootstrap support of .5 => 50% of trees support that node

    input:
    PhyloNode objects, trees is a list of PhyloNode objects
    output: (modified master, bootstrap_supports)
    * modified master_tree (with internally named nodes)
    * bootstrap_supports: list of (node name, bootstrap support)
    """
    setup_master_tree(master_tree)
    for sub_tree in trees:
        compare(master_tree, sub_tree)
    num_trees = len(trees)
    bootstrap_supports = {}
    for node in master_tree.iterNontips(include_self=True):
        node.bootstrap_support = node.bootstrap_support/num_trees
        bootstrap_supports[node.Name] = node.bootstrap_support

    return master_tree, bootstrap_supports

def setup_master_tree(master):
    """ inits bootstrap_support on all nontip nodes, and ensures unique names
    modifies master in place
    """
    i = 0
    for node in master.iterNontips(include_self=True):
        node.bootstrap_support = 0
        if getattr(node, 'Name', None) == None or \
            getattr(node, 'Name', None) == "":
            node.Name = "node" + str(i)
            i += 1
    if len(set(master.getNodeNames())) != len(master.getNodeNames()):
        raise ValueError("can't setup tree, nonunique node names")


def compare(master, subsampled_tree):
    """ compares master tree to subsampled_tree tree, modifies master in place
    """
    if set(master.getTipNames()) != set(subsampled_tree.getTipNames()):
        raise ValueError(\
            "different number of tips in subsampled_tree and master.\n"+\
            "typically some sapmples with few sequences were skipped in"+\
            " the support trees.  remove offending samples from the master"+\
            " tree, and try again")
    # list of lists.  each elem is list of tip names for a node
    subsampled_tree_nodes_names = [node.getTipNames() for node in \
        subsampled_tree.iterNontips(include_self=True)]
    for master_node in master.iterNontips(include_self=True):
        if set(master_node.getTipNames()) in \
            map(set, subsampled_tree_nodes_names):
            master_node.bootstrap_support += 1

def make_cmd_parser():
    """returns command-line options"""
    usage="""
example: python tree_compare.py -m TEST/sample_cluster.tre -s TEST/rare_unifrac_upgma/ -o TEST/unifrac_jackknife/
this makes the folder TEST/unifrac_jackknife.  In that is the master tree,
with internal nodes named, and a separate bootstrap/jackknife support file

master tree must have the same tips as support trees.  if your support trees
omit (e.g.:) samples with few sequences, make a new master tree with those tips
omitted
"""

    parser = OptionParser(usage=usage)
    parser.add_option('-m', '--master', dest='master_tree_path',
        help='master tree path')
    parser.add_option('-s', '--support', dest='support_path',
        help='path to dir containing support trees')
    parser.add_option('-o', '--output', dest='output_path',
        help='output directory, writes two files here')
    options, args = parser.parse_args()
    return options


if __name__ == '__main__':
    options = make_cmd_parser()
    main(options)
