#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Release"

"""tests the tree_compare.py module."""

from cogent.util.unit_test import TestCase, main
from qiime.parse import parse_newick
import qiime.tree_compare as tc

class TreeCompareTests(TestCase):
    """ tests only top level functions
    """
    def test_bootstrap_support(self):
        """ bootstrap_support should have correct bootstrap for a tree with 

        unlabeled internal nodes 
        """
        master_tree = parse_newick('((a:2,b:3):2,(c:1,d:2):7);')
        """
             /-------.5 /-a
        ---1|           \-b
             \------.5 /-c
                       \-d
        """
        t1 = parse_newick('((a:6,b:8.2):2,(c:1,d:2):7);') # same structure
        t2 = parse_newick('((a:2,b:3,c:33):2,d:7);') # abc are siblings
        new_master, bootstraps = tc.bootstrap_support(master_tree, [t1, t2])
        self.assertFloatEqual(sorted(bootstraps.values()),sorted([1.0, .5, .5]))
        
    def test_bootstrap_support_labeled(self):
        """ bootstrap_support should have correct bootstrap on a tree 

        with labeled internal nodes
        """
        master_tree = parse_newick('((a:2,b:3)ab:2,(c:1,d:2)cd:7)rt;')
        """
             /-------.5 /-a
        ---1|           \-b
             \------.5 /-c
                       \-d
        """
        t1 = parse_newick('((a:6,b:8.2)hi:2,(c:1,d:2):7);') # same structure
        t2 = parse_newick('((a:2,b:3,c:33)ho:2,d:7);') # abc are siblings
        new_master, bootstraps = tc.bootstrap_support(master_tree, [t1, t2])
        expected = dict([('ab', .5),('cd',.5),('rt',1.0)])
        self.assertFloatEqual(bootstraps, expected)
        
    def test_bootstrap_support_subset(self):
        """ bootstrap_support should have correct bootstrap on a tree 

        when one support tree is missing a tip
        """
        master_tree = parse_newick('((a:2,b:3)ab:2,(c:1,d:2)cd:7)rt;')
        """
             /-------.5 /-a
        ---1|           \-b
             \------.5 /-c
                       \-d
        """
        t1 = parse_newick('((a:6,b:8.2)hi:2,(c:1,d:2):7);') # same structure
        t2 = parse_newick('((a:2,b:3,c:33)ho:2,d:7);') # abc are siblings
        t3 = parse_newick('((a:6)hi:2,(c:1,d:2):7);') # b missing
        t4 = parse_newick('(a:8,(c:1,d:2):7);') # b missing, and pruned
        new_master, bootstraps = tc.bootstrap_support(master_tree, 
            [t1, t2,t3,t4])
        expected = dict([('cd',.75),('rt',1.0)])
        self.assertFloatEqual(bootstraps, expected)
        
    def test_tree_support(self):
        """ tree_support should correctly modify node.bootstrap_support
        """
        master_tree = parse_newick('((a:2,b:3)ab:2,(c:1,d:2)cd:7)rt;')
        """
             /-------.5 /-a
        ---1|           \-b
             \------.5 /-c
                       \-d
        """
        t2 = parse_newick('((a:2,b:3,c:33)ho:2,d:7);') # abc are siblings
        
        tc.tree_support(master_tree, t2)
        self.assertFloatEqual(\
            master_tree.getNodeMatchingName('rt').bootstrap_support,1.0) 
        
if __name__ =='__main__':
    main()
