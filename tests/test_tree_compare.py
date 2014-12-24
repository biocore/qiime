#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.9.0-rc1"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

"""tests the tree_compare.py module."""

from unittest import TestCase, main
from numpy.testing import assert_almost_equal
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
        t1 = parse_newick('((a:6,b:8.2):2,(c:1,d:2):7);')  # same structure
        t2 = parse_newick('((a:2,b:3,c:33):2,d:7);')  # abc are siblings
        new_master, bootstraps = tc.bootstrap_support(master_tree, [t1, t2])
        assert_almost_equal(
            sorted(bootstraps.values()), sorted([1.0, .5, .5]))

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
        t1 = parse_newick('((a:6,b:8.2)hi:2,(c:1,d:2):7);')  # same structure
        t2 = parse_newick('((a:2,b:3,c:33)ho:2,d:7);')  # abc are siblings
        new_master, bootstraps = tc.bootstrap_support(master_tree, [t1, t2])
        expected = dict([('ab', .5), ('cd', .5), ('rt', 1.0)])
        self.assertDictEqual(bootstraps, expected)

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
        t1 = parse_newick('((a:6,b:8.2)hi:2,(c:1,d:2):7);')  # same structure
        t2 = parse_newick('((a:2,b:3,c:33)ho:2,d:7);')  # abc are siblings
        t3 = parse_newick('((a:6)hi:2,(c:1,d:2):7);')  # b missing
        t4 = parse_newick('(a:8,(c:1,d:2):7);')  # b missing, and pruned
        new_master, bootstraps = tc.bootstrap_support(master_tree,
                                                      [t1, t2, t3, t4])
        expected = dict([('cd', .75), ('rt', 1.0)])
        self.assertDictEqual(bootstraps, expected)

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
        t2 = parse_newick('((a:2,b:3,c:33)ho:2,d:7);')  # abc are siblings

        tc.tree_support(master_tree, t2)
        assert_almost_equal(
            master_tree.getNodeMatchingName('rt').bootstrap_support, 1.0)

    def test_setup_master_tree_alltips(self):
        """tests setup_master_tree"""
        master_tree = parse_newick('((a:2,b:3):2,(c:1,d:2)foo:7);')
        t1 = parse_newick('((a:6,b:8.2):2,(c:1,d:2):7);')  # same structure
        t2 = parse_newick('((a:2,b:3,c:33):2,d:7);')  # abc are siblings
        support_tree_tipnames = ['a', 'b', 'c', 'd']
        exp = "((a:2.0,b:3.0)node1:2.0,(c:1.0,d:2.0)foo:7.0)node0;"
        new_master = tc.setup_master_tree(master_tree, [t1, t2])
        self.assertEqual(new_master.getNewick(with_distances=True), exp)

        desc_tips_root = frozenset(['a', 'b', 'c', 'd'])
        desc_tips_c0 = frozenset(['a', 'b'])
        desc_tips_c1 = frozenset(['c', 'd'])
        desc_tips_c0c0 = frozenset(['a'])
        desc_tips_c0c1 = frozenset(['b'])
        desc_tips_c1c0 = frozenset(['c'])
        desc_tips_c1c1 = frozenset(['d'])

        self.assertEqual(frozenset(new_master.getTipNames()), desc_tips_root)
        self.assertEqual(
            frozenset(new_master.Children[0].getTipNames()),
            desc_tips_c0)
        self.assertEqual(
            frozenset(new_master.Children[1].getTipNames()),
            desc_tips_c1)
        c0 = master_tree.Children[0]
        c1 = master_tree.Children[1]
        self.assertEqual(
            frozenset(c0.Children[0].getTipNames()),
            desc_tips_c0c0)
        self.assertEqual(
            frozenset(c0.Children[1].getTipNames()),
            desc_tips_c0c1)
        self.assertEqual(
            frozenset(c1.Children[0].getTipNames()),
            desc_tips_c1c0)
        self.assertEqual(
            frozenset(c1.Children[1].getTipNames()),
            desc_tips_c1c1)

    def test_setup_master_tree_missingtips(self):
        """tests setup_master_tree"""
        master_tree = parse_newick('((a:2,b:3):2,(c:1,d:2)foo:7);')
        t1 = parse_newick('((a:6,b:8.2):2,(c:1,d:2):7);')  # same structure
        t2 = parse_newick('((a:2,c:33):2,d:7);')  # where's b?
        support_tree_tipnames = ['a', 'c', 'd']
        exp = "((c:1.0,d:2.0)foo:7.0,a:4.0)node0;"
        new_master = tc.setup_master_tree(master_tree, [t1, t2])
        self.assertEqual(new_master.getNewick(with_distances=True), exp)

        desc_tips_root = frozenset(['a', 'c', 'd'])
        desc_tips_c1 = frozenset(['a'])
        desc_tips_c0 = frozenset(['c', 'd'])
        desc_tips_c0c0 = frozenset(['c'])
        desc_tips_c0c1 = frozenset(['d'])

        self.assertEqual(frozenset(new_master.getTipNames()), desc_tips_root)
        self.assertEqual(
            frozenset(new_master.Children[0].getTipNames()),
            desc_tips_c0)
        self.assertEqual(
            frozenset(new_master.Children[1].getTipNames()),
            desc_tips_c1)
        c0 = new_master.Children[0]
        self.assertEqual(
            frozenset(c0.Children[0].getTipNames()),
            desc_tips_c0c0)
        self.assertEqual(
            frozenset(c0.Children[1].getTipNames()),
            desc_tips_c0c1)

if __name__ == '__main__':
    main()
