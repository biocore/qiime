#!/usr/bin/env python

"""additions to cogent.core.tree. 

Note: these changes are in cogent svn but correctly under the PhyloNode
object. These cannot be added as member methods to the PhyloNode in here
because it would require massive awkward changes and adding in a _lot_ of
modified cogent files into pycogent_backports
"""

from numpy import argsort

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Release"

def setMaxTipTipDistance(tree):
    """Propagate tip distance information up the tree

    This method was originally implemented by Julia Goodrich with the intent
    of being able to determine max tip to tip distances between nodes on 
    large trees efficiently. The code has been modified to track the 
    specific tips the distance is between
    """
    for n in tree.postorder():
        if n.isTip():
            n.MaxDistTips = [[0.0, n.Name], [0.0, n.Name]]
        else:
            if len(n.Children) == 1:
                tip_a, tip_b = n.Children[0].MaxDistTips
                tip_a[0] += n.Children[0].Length or 0.0
                tip_b[0] += n.Children[0].Length or 0.0
            else:
                tip_info = [(max(c.MaxDistTips), c) for c in n.Children]
                dists = [i[0][0] for i in tip_info]
                best_idx = argsort(dists)[-2:]
                tip_a, child_a = tip_info[best_idx[0]]
                tip_b, child_b = tip_info[best_idx[1]]
                tip_a[0] += child_a.Length or 0.0
                tip_b[0] += child_b.Length or 0.0
            n.MaxDistTips = [tip_a, tip_b]

def getMaxTipTipDistance(tree):
    """Returns the max tip tip distance between any pair of tips
        
    Returns (dist, tip_names, internal_node)
    """
    if not hasattr(tree, 'MaxDistTips'):
        setMaxTipTipDistance(tree)

    longest = 0.0
    names = [None,None]
    best_node = None
    for n in tree.nontips(include_self=True):
        tip_a, tip_b = n.MaxDistTips
        dist = (tip_a[0] + tip_b[0])

        if dist > longest:
            longest = dist
            best_node = n
            names = [tip_a[1], tip_b[1]]
    return longest, names, best_node

