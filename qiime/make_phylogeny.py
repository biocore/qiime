#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Rob Knight", "Justin Kuczynski", "Daniel McDonald",
               "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"


"""Contains code for making a phylogeny from aligned sequences, using several techniques.

This module has the responsibility for taking a set of sequences and
returning an alignment. Mostly, it will be thin wrappers for code
already in cogent.app.*, to which wrappers for e.g. NAST need to be
added..
"""

from skbio.parse.sequences import parse_fasta
from skbio.alignment import Alignment
from skbio.sequence import DNA

from qiime.util import FunctionWithParams
# app controllers that implement align_unaligned_seqs
import bfillings.muscle_v38
import bfillings.clustalw
import bfillings.mafft
import bfillings.raxml_v730
import bfillings.fasttree
import bfillings.clearcut


class TreeBuilder(FunctionWithParams):

    """A TreeBuilder takes a aligned set of sequences and returns a tree.

    This is an abstract class: subclasses should implement the __call__
    method.

    Note: sequence ids should be preserved during this process, i.e. the
    description lines should be saved/restored if the phylogeny app is
    destructive to them.

    Specific wrappers need to be written for nontraditional approaches,
    including:

    (a) RAxML assignment, where sequences are assigned to internal nodes of
        an existing tree
    (b) BLAST-based assignment, where sequences are assigned to existing
        nodes of a tree based on best blast hit (assigned only to terminal
        nodes if using a single hit, but could easily imagine assigning to
        an internal node based on a set of indistinguishably good hits).
    (c) BLAST-like approach using VMATCH or other suffix array library, or
        using oligonucleotide freqs like the RDP classifier does to assign
        to an arbitrary internal node in an existing tree, etc.
    """
    Name = 'TreeBuilder'

    def __init__(self, params):
        """Return new TreeBuilder object with specified params.

        Note: expect params to contain both generic and per-method (e.g. for
        raxml vs. fasttree vs. whatever) params, so leaving it as a dict
        rather than setting attributes. Some standard entries in params are:

        Application: 3rd-party application used, if any, e.g. raxml
        """
        self.Params = params

    def __call__(self, aln_path, result_path=None, log_path=None):
        """Returns tree from alignment.

        Parameters:
        aln_path: path to file of aligned sequences
        result_path: path to file of results. If specified, should
        dump the result to the desired path as fasta, otherwise should
        return skbio.core.tree.TreeNode object.
        log_path: path to log, which should include dump of params.
        """
        raise NotImplementedError("TreeBuilder is an abstract class")


class CogentTreeBuilder(TreeBuilder):

    """Generic tree builder using Cogent tree methods."""

    Name = 'CogentTreeBuilder'

    def getResult(self, aln_path, *args, **kwargs):
        """Returns alignment from sequences.

        Currently does not allow parameter tuning of program and uses
        default parameters -- this is bad and should be fixed.

        #TODO: allow command-line access to important aln params.
        """
        module = self.Params['Module']
        # standard qiime says we just consider the first word as the unique ID
        # the rest of the defline of the fasta alignment often doesn't match
        # the otu names in the otu table
        with open(aln_path) as aln_f:
            seqs = Alignment.from_fasta_records(
                parse_fasta(aln_f, label_to_name=lambda x: x.split()[0]),
                DNA)
        # This ugly little line of code lets us pass a skbio Alignment when a
        # a cogent alignment is expected.
        seqs.getIntMap = seqs.int_map
        result = module.build_tree_from_alignment(seqs)

        try:
            root_method = kwargs['root_method']
            if root_method == 'midpoint':
                result = root_midpt(result)
            elif root_method == 'tree_method_default':
                pass
        except KeyError:
            pass
        return result

    def __call__(self, result_path=None, log_path=None, *args, **kwargs):
        """Calls superclass method to align seqs"""
        return FunctionWithParams.__call__(self, result_path=result_path,
                                           log_path=log_path, *args, **kwargs)

tree_method_constructors = {}
tree_module_names = {'muscle': bfillings.muscle_v38,
                     'clustalw': bfillings.clustalw,
                     #'mafft':bfillings.mafft,
                     # current version of Mafft does not support tree building
                     'fasttree': bfillings.fasttree,
                     'raxml_v730': bfillings.raxml_v730,
                     'clearcut': bfillings.clearcut
                     }

# def maxTipTipDistance(tree):
#        """returns the max distance between any pair of tips
#
#        Also returns the tip names  that it is between as a tuple"""
#        distmtx, tip_order = tree.tipToTipDistances()
#        idx_max = divmod(distmtx.argmax(),distmtx.shape[1])
#        max_pair = (tip_order[idx_max[0]].Name, tip_order[idx_max[1]].Name)
#        return distmtx[idx_max], max_pair
#
# def decorate_max_tip_to_tip_distance(self):
#    """Propagate tip distance information up the tree
#
#    This method was originally implemented by Julia Goodrich with the intent
#    of being able to determine max tip to tip distances between nodes on large
#    trees efficiently. The code has been modified to track the specific tips
#    the distance is between
#    """
#    for n in self.postorder():
#        if n.isTip():
#            n.MaxDistTips = [(0.0, n.Name), (0.0, n.Name)]
#        else:
#            n.MaxDistTips = [max(c.MaxDistTips) for c in n.Children][:2]
#
#    max_dist = 0.0
#    max_names = [None, None]
#    for n in tree.nontips():
#        tip_a, tip_b = n.MaxTipsTips
#        dist = tip_a[0] + tip_b[0]
#        if dist > max_dist:
#            max_dist = dist
#            max_names = [tip_a[1], tip_b[1]]
#
#    return max_dist, max_names


def root_midpt(tree):
    """ this was instead of PhyloNode.rootAtMidpoint(), which is slow and broke

    this should be deprecated in a future release once the release version
    of PyCogent's tree.rootAtMidpoint() is identical to this function

    this fn doesn't preserve the internal node naming or structure,
    but does keep tip to tip distances correct.  uses unrootedDeepcopy()
    """
    #max_dist, tip_names, int_node = getMaxTipTipDistance(tree)
    max_dist, tip_names, int_node = tree.getMaxTipTipDistance()

    half_max_dist = max_dist / 2.0
    if max_dist == 0.0:  # only pathological cases with no lengths
        return tree.unrootedDeepcopy()
    tip1 = tree.getNodeMatchingName(tip_names[0])
    tip2 = tree.getNodeMatchingName(tip_names[1])
    # last comm ancestor
    lca = tree.getConnectingNode(tip_names[0], tip_names[1])
    if tip1.distance(lca) > half_max_dist:
        climb_node = tip1
    else:
        climb_node = tip2

    dist_climbed = 0.0
    while dist_climbed + climb_node.Length < half_max_dist:
        dist_climbed += climb_node.Length
        climb_node = climb_node.Parent

    # now midpt is either at on the branch to climb_node's  parent
    # or midpt is at climb_node's parent
    # print dist_climbed, half_max_dist, 'dists cl hamax'
    if dist_climbed + climb_node.Length == half_max_dist:
        # climb to midpoint spot
        climb_node = climb_node.Parent
        if climb_node.isTip():
            raise RuntimeError('error trying to root tree at tip')
        else:
            # print climb_node.Name, 'clmb node'
            return climb_node.unrootedDeepcopy()

    else:
        # make a new node on climb_node's branch to its parent
        tmp_node_name = "TMP_ROOT_NODE_NAME"
        parent = climb_node.Parent
        parent.removeNode(climb_node)
        climb_node.Parent = None
        new_node = parent.__class__()
        new_node.Name = tmp_node_name

        # adjust branch lengths
        old_br_len = climb_node.Length
        climb_node.Length = half_max_dist - dist_climbed
        new_node.Length = old_br_len - climb_node.Length

        if climb_node.Length < 0.0 or new_node.Length < 0.0:
            raise RuntimeError(
                'attempting to create a negative branch length!')

        # re-attach tree
        parent.append(new_node)
        new_node.append(climb_node)

        # reroot and remove the temporary node name
        new_tree = tree.rootedAt(tmp_node_name)
        new_root = new_tree.getNodeMatchingName(tmp_node_name)
        new_root.Name = None

        return new_tree
