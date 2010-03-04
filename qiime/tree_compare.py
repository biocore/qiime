#!/usr/bin/env python
from __future__ import division
from optparse import OptionParser
from cogent.core.tree import PhyloNode
from cogent.parse.tree import DndParser
import qiime.parse
import os.path
import sys
from warnings import warn
import copy

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "0.92-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Pre-release"

def load_tree_files(master_tree_file, support_dir):
    """Load trees from filepaths
    
    checks if support filenames indicate that support are from different 
    distance methods.  If so, warns user.
    loads trees into phylonode objects
    returns master_tree, [support_trees]
    """
    tree_file_names = os.listdir(support_dir)
    # ignore invisible files like .DS_Store
    tree_file_names = [fname for fname in tree_file_names if not \
        fname.startswith('.')]

    #try to warn user if using multiple types of trees    
    try:
        base_names = []
        for fname in tree_file_names:
            base_names.append(qiime.parse.parse_rarefaction_fname(fname)[0])
    except ValueError:
        pass
    else:
        if len(set(base_names)) > 1:
            warnstr = """
warning: support trees are named differently, please be sure you're not 
comparing trees generated in different manners, unless you're quite sure 
that's what you intend to do.  types: """ + str(set(base_names)) + """
continuing anyway..."""
            warn(warnstr)

    master_tree = DndParser(open(master_tree_file, 'U'), PhyloNode)
    support_trees = []
    for fname in tree_file_names:
        try:
            f = open(os.path.join(support_dir, fname), 'U')
            tree = DndParser(f, PhyloNode)
            tree.filepath = fname
            support_trees.append(tree)
            f.close()
        except IOError, err:
            sys.stderr.write('error loading support tree ' + fname + '\n')
            exit(1)
    return master_tree, support_trees
    
def write_bootstrap_support_files(master_tree, bootstraps, output_dir,
    num_support_trees):
    """write bootsrap support data to 3 files in output_dir
    
    assumes input master tree has unique node names
    writes that tree to file using newick format, along with a tab delimited
    text file listing the bootstrap support of each internal node in that tree.
    Also writes a pseudo-newick file with internal nodes named by their
    bootstrap values"""
    
    

    # master tree as passed
    fname = os.path.join(output_dir, "master_tree.tre")
    f = open(fname, 'w')
    f.write(master_tree.getNewick(with_distances=True))
    f.close()
    
    # support of nodes in tab delimited text
    f = open(os.path.join(output_dir, 'jackknife_support.txt'), 'w')
    f.write('#total support trees considered: ' + str(num_support_trees) +'\n')
    f.write('#node support is fractional - in range [0,1]\n')
    for key, val in bootstraps.items():
        f.write("\t".join([str(key),str(val)]) + "\n")
    f.close()
    
    # tree with nodes named by support values
    pseudo_newick_master = copy.deepcopy(master_tree)
    f = open(os.path.join(output_dir, 'jackknife_named_nodes.tre'), 'w')
    for name, support_val in bootstraps.items():
        node = pseudo_newick_master.getNodeMatchingName(name)
        node.Name = str(support_val)
    f.write(pseudo_newick_master.getNewick(with_distances=True))
    f.close()
    
def bootstrap_support(master_tree, trees):
    """ calculate bootstrap/jackknife support of master, by trees

    this calculates bootstrap support of each nontip node in master_tree
    a tree supports a given master_tree_node if there exists a node in tree
    where node.tips == master_tree_node.tips (by name)
    not specific to bootstrap, does node support for trees generated in any
    manner (e.g.: jackknifing)

    bootstrap support of .5 => 50% of trees support that node

    input:
    PhyloNode objects, trees is a list of PhyloNode objects
    output: (modified master, bootstrap_supports)
    * new master_tree, modified with internally named nodes
    * bootstrap_supports: list of (node name, bootstrap support)
    """
    new_master = setup_master_tree(master_tree)
    for sub_tree in trees:
        # modifies new_master in place
        tree_support(new_master, sub_tree)
    num_trees = len(trees)
    bootstrap_supports = {}
    for node in new_master.iterNontips(include_self=True):
        node.bootstrap_support = node.bootstrap_support/num_trees
        bootstrap_supports[node.Name] = node.bootstrap_support

    return new_master, bootstrap_supports

def setup_master_tree(master):
    """ inits bootstrap_support on all nontip nodes, and ensures unique names
    
    returns a nearly identical copy, with uniquely named internal nodes,
    and node.bootstrap_support set to 0
    """
    new_master = copy.deepcopy(master)
    i = 0
    for node in new_master.iterNontips(include_self=True):
        node.bootstrap_support = 0
        if getattr(node, 'Name', None) == None or \
            getattr(node, 'Name', None) == "":
            node.Name = "node" + str(i)
            i += 1
    if len(set(new_master.getNodeNames())) != len(new_master.getNodeNames()):
        node_names = master.getNodeNames()
        
        nonuniques = []
        for name in set(node_names):
            if node_names.count(name) > 1:
                nonuniques.append(name)
        raise ValueError("can't setup master tree, nonunique node names" +
            str(nonuniques))

    return new_master

def tree_support(master, subsampled_tree):
    """ compares master tree to subsampled_tree, modifies master in place
    
    this calculates bootstrap support of each nontip node in the master tree
    a given master_tree_node is supported if there exists a node in subsampled
    tree where sub_tree_node.tips == master_tree_node.tips (by name)
    
    not specific to bootstrap, does node support for trees generated in any
    manner (e.g.: jackknifing)
    master is modified to have node.bootstrap_support incremented by 1 if
    subsampled tree has support for that node
    """
    if set(master.getTipNames()) != set(subsampled_tree.getTipNames()):
        raise ValueError("""
problem with subsampled tree """ + subsampled_tree.filepath + """
different number of tips in subsampled_tree and master.\n
typically some samples with few sequences were skipped in
the support trees.  be sure to keep all samples when doing rarefaction of
otu tables (see --small_included option), or
remove offending samples from the master tree, and try again)
""")
    # subsampled_tree_nodes_names is a list of lists.
    # each elem is list of tip names for a node
    subsampled_tree_nodes_names = []
    for node in subsampled_tree.iterNontips(include_self=True):
        subsampled_tree_nodes_names.append(node.getTipNames())
    
    #now a list of sets, each set is tip names for a specific node
    subsampled_tree_nodes_names = map(set, subsampled_tree_nodes_names)
        
    for master_node in master.iterNontips(include_self=True):
        if set(master_node.getTipNames()) in subsampled_tree_nodes_names:
            master_node.bootstrap_support += 1
