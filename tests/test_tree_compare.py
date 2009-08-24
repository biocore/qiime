#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Justin Kuczynski"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Prototype"

"""tests the tree_compare.py module."""

from cogent.util.unit_test import TestCase, main
from cogent.maths.unifrac.fast_unifrac import DndParser
from pipe454.tree_compare import bootstrap_support

class TreeCompareTests(TestCase):
    """ tests only top level functions
    """
    def test_bootstrap_support(self):
        """ tree with unlabeled internal nodes should have correct bootstrap
        values"""
        master_tree = DndParser('((a:2,b:3):2,(c:1,d:2):7);')
        """
             /-------.5 /-a
        ---1|           \-b
             \------.5 /-c
                       \-d
        """
        t1 = DndParser('((a:6,b:8.2):2,(c:1,d:2):7);') # same structure
        t2 = DndParser('((a:2,b:3,c:33):2,d:7);') # abc are siblings
        new_master, bootstraps = bootstrap_support(master_tree, [t1, t2])
        self.assertFloatEqual(sorted(bootstraps.values()),sorted([1.0, .5, .5]))
        
    def test_bootstrap_support_labeled(self):
        """ tree with labeled internal nodes should have correct bootstrap
        values"""
        master_tree = DndParser('((a:2,b:3)ab:2,(c:1,d:2)cd:7)rt;')
        """
             /-------.5 /-a
        ---1|           \-b
             \------.5 /-c
                       \-d
        """
        t1 = DndParser('((a:6,b:8.2)hi:2,(c:1,d:2):7);') # same structure
        t2 = DndParser('((a:2,b:3,c:33)ho:2,d:7);') # abc are siblings
        new_master, bootstraps = bootstrap_support(master_tree, [t1, t2])
        expected = dict([('ab', .5),('cd',.5),('rt',1.0)])
        self.assertFloatEqual(bootstraps, expected)
        
if __name__ =='__main__':
    main()
