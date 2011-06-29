#!/usr/bin/env python

from cogent.core.tree import PhyloNode
from cogent.util.unit_test import TestCase, main
from qiime.pycogent_backports.tree import getMaxTipTipDistance, \
        setMaxTipTipDistance

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Daniel McDonald"
__email__ = "mcdonadt@colorado.edu"
__status__ = "Release"

class TreeTests(TestCase):
    def setUp(self):
        # setup tree here copied from cogent test_core.test_tree for PhyloNode
        # most likely written by Rob Knight
        nodes = dict([(x, PhyloNode(x)) for x in 'abcdefgh'])
        nodes['a'].append(nodes['b'])
        nodes['b'].append(nodes['c'])
        nodes['c'].append(nodes['d'])
        nodes['c'].append(nodes['e'])
        nodes['c'].append(nodes['f'])
        nodes['f'].append(nodes['g'])
        nodes['a'].append(nodes['h'])
        self.TreeNode = nodes
        self.TreeRoot = nodes['a']
        nodes['a'].Length = None
        nodes['b'].Length = 0
        nodes['c'].Length = 3
        nodes['d'].Length = 1
        nodes['e'].Length = 4
        nodes['f'].Length = 2
        nodes['g'].Length = 3
        nodes['h'].Length = 2 

    def test_getMaxTipTipDistance(self):
        """getMaxTipTipDistance should get max tip distance across tree"""
        nodes, tree = self.TreeNode, self.TreeRoot
        dist, names, node = getMaxTipTipDistance(tree)
        self.assertEqual(dist, 15.0) # due to nodes with single descendents!!
        self.assertEqual(sorted(names), ['e','g'])
        self.assertEqual(node.Name, 'b')

    def test_setMaxTipTipDistance(self):
        """setMaxTipTipDistance sets MaxDistTips across tree"""
        nodes, tree = self.TreeNode, self.TreeRoot
        setMaxTipTipDistance(tree)
        tip_a, tip_b = tree.MaxDistTips
        self.assertEqual(tip_a[0] + tip_b[0], 10)
        self.assertEqual(sorted([tip_a[1],tip_b[1]]), ['g','h'])

if __name__ == '__main__':
    main()
