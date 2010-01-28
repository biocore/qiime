#!/usr/bin/env python

from cogent.util.unit_test import TestCase, main
from cogent.parse.tree import DndParser
from collapse_taxontree_nodes import CollapseNode, ncbi_taxonid_filter, \
        pick_tree_OTUs, get_tip_relatives

__author__ = "Daniel McDonald"
__copyright__ = "Copyright 2010, The QIIME Project" #consider project name
__credits__ = ["Daniel McDonald"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Daniel McDonald"
__email__ = "daniel.mcdonald@colorado.edu"
__status__ = "Pre-release"

class CollapseNodeTests(TestCase):
    def setUp(self):
        t_str = "(((a:0.1,b:0.2)c:0.1,d:0.2)e:0.2,(f:0.1,g:0.1)h:0.4,i:0.2)r;"
        self.tree = DndParser(t_str, constructor=CollapseNode)

    def test_init(self):
        self.assertEqual(self.tree.AtThreshold, 0.0)
        self.assertEqual(self.tree.TaxonId, None)
        self.assertEqual(self.tree.TaxonomyInformation, {})

    def test_setRepresentativeTaxonId(self):
        idx = self.tree.indexByAttr('Name')
        idx['a'].TaxonId = 1
        idx['b'].TaxonId = 1
        idx['d'].TaxonId = 2

        idx['f'].TaxonId = None
        idx['g'].TaxonId = None
        idx['i'].TaxonId = 1

        idx['c'].setRepresentativeTaxonId()
        self.assertEqual(idx['c'].TaxonId, 1)
        idx['e'].setRepresentativeTaxonId()
        self.assertEqual(idx['e'].TaxonId, 1)
        idx['h'].setRepresentativeTaxonId()
        self.assertEqual(idx['h'].TaxonId, None)
        idx['r'].setRepresentativeTaxonId()
        self.assertEqual(idx['r'].TaxonId, 1)

    def test_setTaxonInfoAtTips(self):
        tax_filter = lambda x: x
        tax_lookup = {'a':1,'d':2,'g':3,'i':4}
        self.tree.setTaxonInfoAtTips(tax_lookup, tax_filter)
        self.assertEqual(self.tree.getNodeMatchingName('a').TaxonId, 1)
        self.assertEqual(self.tree.getNodeMatchingName('b').TaxonId, None)
        self.assertEqual(self.tree.getNodeMatchingName('d').TaxonId, 2)
        self.assertEqual(self.tree.getNodeMatchingName('f').TaxonId, None)
        self.assertEqual(self.tree.getNodeMatchingName('g').TaxonId, 3)
        self.assertEqual(self.tree.getNodeMatchingName('i').TaxonId, 4)

class CollapseTaxonTreeNodesTests(TestCase):
    def setUp(self):
        t_str = "(((a:0.1,b:0.2)c:0.1,d:0.2)e:0.2,(f:0.1,g:0.1)h:0.4,i:0.2)r;"
        self.tree = DndParser(t_str, constructor=CollapseNode)
    def test_ncbi_taxonid_filter(self):
        good = {'a':1,'b':2,'ncbi_tax_string_format_2':'foo','ncbi_tax_id':1}
        bad = {'a':1,'b':2,'ncbi_tax_string_format_2':'UnKNoWn','ncbi_tax_id':2}
        exp_good = 1
        exp_bad = None
        obs_good = ncbi_taxonid_filter(good)
        obs_bad = ncbi_taxonid_filter(bad)
        self.assertEqual(obs_good, exp_good)
        self.assertEqual(obs_bad, exp_bad)

    def test_get_tip_relatives(self):
        """Test if we are getting the deepest node below threshold"""
        idx = self.tree.indexByAttr('Name')
        idx['a'].TaxonId = 1
        idx['b'].TaxonId = None
        idx['d'].TaxonId = 1
        tip_a = self.tree.getNodeMatchingName('a')
        exp_at_005 = self.tree.getNodeMatchingName('a')
        exp_at_011 = self.tree.getNodeMatchingName('a') 
        exp_at_016 = self.tree.getNodeMatchingName('c')
        exp_at_024 = self.tree.getNodeMatchingName('e')

        obs_005 = get_tip_relatives(tip_a, 0.05)
        obs_011 = get_tip_relatives(tip_a, 0.11)
        obs_016 = get_tip_relatives(tip_a, 0.16)
        obs_024 = get_tip_relatives(tip_a, 0.24)

        self.assertEqual(obs_005, exp_at_005)
        self.assertEqual(obs_011, exp_at_011)
        self.assertEqual(obs_016, exp_at_016)
        self.assertEqual(obs_024, exp_at_024)

    def test_pick_tree_OTUs(self):
        """Test if we are getting the correct OTUs on a given tree"""
        exp_nodes = [self.tree.getNodeMatchingName('c'),
                     self.tree.getNodeMatchingName('d'),
                     self.tree.getNodeMatchingName('h'),
                     self.tree.getNodeMatchingName('i')]
        for n in self.tree.tips():
            n.TaxonId = 1
        pick_tree_OTUs(self.tree, 0.16)
        self.assertEqualItems(exp_nodes, self.tree.tips())

if __name__ == '__main__':
    main()
