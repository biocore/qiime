#!/usr/bin/env python
from __future__ import division
from optparse import OptionParser
from cogent.core.tree import PhyloNode
from cogent.parse.tree import DndParser
import qiime.parse
import os.path
import sys
from warnings import warn

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2009, the PyCogent Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Prototype"

"""compares jackknife/bootstrapped trees with a master, and outputs support for
nodes.  main() is called when run from command line, primary function of this
module is bootstrap_support()

see command line usage example called with -h
"""


def main(options):
    tree_file_names = os.listdir(options.support_dir)
    tree_file_names = [fname for fname in tree_file_names if not \
        fname.startswith('.')]
    try:
        #try to warn user if using multiple types of trees
        base_names = [qiime.parse.parse_rarefaction_fname(fname)[0]
            for fname in tree_file_names]
        if len(set(base_names)) > 1:
            warnstr = """
warning: support trees are named differently, please be sure you're not 
comparing trees generated in different manners, unless you're quite sure 
that's what you intend to do.  types: """ + str(set(base_names))
            warn(warnstr)
    except:
        pass

    support_trees = []
    for fname in tree_file_names:
        try:
            f = open(os.path.join(options.support_dir, fname))
            tree = DndParser(f, PhyloNode)
            tree.filepath = fname
            support_trees.append(tree)
            f.close()
        except IOError, err:
            sys.stderr.write('error loading support tree ' + fname + '\n')
            exit(1)
    master_tree = DndParser(open(options.master_tree), PhyloNode)

    # get support of each node in master
    new_master, bootstraps = bootstrap_support(master_tree, support_trees)

    # write out modified master tree with internal nodes, and bootstrap support
    # values to 2 files
    fname = os.path.join(options.output_dir, "master_tree.tre")
    f = open(fname, 'w')
    f.write(new_master.getNewick(with_distances=True))
    f.close()
    
    f = open(os.path.join(options.output_dir, 'jackknife_support.txt'), 'w')
    f.write('total support trees considered: ' + str(len(support_trees)) + '\n')
    for key, val in bootstraps.items():
        f.write("\t".join([str(key),str(val)]) + "\n")
    f.close()
    
    # we ruin master_tree here, this better be the end
    f = open(os.path.join(options.output_dir, 'jackknife_named_nodes.tre'), 'w')
    for name, support_val in bootstraps.items():
        node = new_master.getNodeMatchingName(name)
        node.Name = str(support_val)
    f.write(new_master.getNewick(with_distances=True))
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
        raise ValueError("""
didn't like subsampled tree """ + subsampled_tree.filepath + """
different number of tips in subsampled_tree and master.\n
typically some samples with few sequences were skipped in
the support trees.  be sure to keep all samples when doing rarefaction of
otu tables (see --small_included option), or
remove offending samples from the master tree, and try again)
""")
    # list of lists.  each elem is list of tip names for a node
    subsampled_tree_nodes_names = [node.getTipNames() for node in \
        subsampled_tree.iterNontips(include_self=True)]
    for master_node in master.iterNontips(include_self=True):
        if set(master_node.getTipNames()) in \
            map(set, subsampled_tree_nodes_names):
            master_node.bootstrap_support += 1
            
            
usage_str = """usage: %prog [options] {-i INPUT_PATH -o OUTPUT_DIR -m MASTER_TREE -s SUPPORT_DIR}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:
python %prog -m sample_cluster.tre -s rare_unifrac_upgma/ -o unifrac_jackknife/
this makes the folder unifrac_jackknife.  In that is the master tree,
with internal nodes named uniquely, a separate bootstrap/jackknife support file,
and a jackknife_named_nodes.tre tree, for use with e.g.: figtree
output jackknife support values 

master tree must have the same tips as support trees.  if your support trees
omit some tips (e.g.: samples with few sequences),
make a new master tree with those tips omitted
"""
def parse_command_line_parameters():
    """returns command-line options"""

    if len(sys.argv) == 1:
        sys.argv.append('-h')
    usage = usage_str
    version = '%prog ' + str(__version__)
    parser = OptionParser(usage=usage, version=version)

    parser.add_option('-m', '--master_tree',
        help='master tree filepath [REQUIRED]')

    parser.add_option('-s', '--support_dir',
        help='path to dir containing support trees [REQUIRED]')

    parser.add_option('-o', '--output_dir',
        help='output directory, writes two files here'+\
        "makes dir if it doesn't exist [REQUIRED]")

    opts, args = parser.parse_args()
    if len(args) != 0:
        parser.error("positional argument detected.  make sure all"+\
         ' parameters are identified.' +\
         '\ne.g.: include the \"-m\" in \"-m MINIMUM_LENGTH\"')
         
    required_options = ['master_tree','support_dir','output_dir']
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 
    return opts, args


if __name__ == '__main__':
    opts,args = parse_command_line_parameters()
    if not os.path.exists(opts.output_dir):
    	os.makedirs(opts.output_dir)
    main(opts)
