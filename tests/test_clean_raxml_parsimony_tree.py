#!/usr/bin/env python
# File created on 16 Nov 2011
from __future__ import division

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.9.0"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"


from unittest import TestCase, main
from StringIO import StringIO
from cogent.parse.tree import DndParser
from cogent.core.tree import PhyloNode, TreeNode
from qiime.clean_raxml_parsimony_tree import decorate_numtips, decorate_depth,\
    get_insert_dict, drop_duplicate_nodes


class TopLevelTests(TestCase):

    """Tests of top-level functions"""

    def setUp(self):
        """define some top-level data"""

        self.tree = DndParser(StringIO(TEST_TREE), constructor=PhyloNode)

    def test_decorate_numtips(self):
        """decorate_numtips: decorate the number of tips below each node."""

        obs = decorate_numtips(self.tree)

        # make sure each tip only has 1 tip
        tips = obs.tips()
        for obs_tips in tips:
            exp = 1
            self.assertEqual(obs_tips.Score, exp)

        # check that the node scores are correct
        node_numtips = [1, 2, 8]
        for nodes in obs:
            self.assertTrue(nodes.Score in node_numtips)

    def test_decorate_depth(self):
        """decorate_depth: decorate the depth from the root each node is."""

        # decorate the depth of each node on the tree
        obs = decorate_depth(self.tree)

        # make sure each tip depth is between 1 and 5
        tips = obs.tips()
        tip_depth = [1, 2, 3, 4, 5]
        for obs_tips in tips:
            self.assertTrue(obs_tips.Score in tip_depth)

        # check that the node depth is 1
        for nodes in obs:
            exp = 1
            self.assertEqual(nodes.Score, exp)

    def test_get_insert_dict(self):
        """get_insert_dict: get the location of each inserted tip."""

        # pull out the tips for Species006
        obs = get_insert_dict(self.tree, 'Species006')

        # verify the dict is properly refrenced
        self.assertTrue('Species006' in obs)

        # if Species006 in dict, verify there are 3 tips
        if 'Species006' in obs:
            exp_len = 3
            obs_len = len(obs['Species006'])
            self.assertEqual(obs_len, exp_len)

    def test_drop_duplicate_nodes(self):
        """drop_duplicate_nodes: remove duplicate tips from tree based on either the number of tips or depth."""

        # pull out the tips for Species006
        inserted_nodes = get_insert_dict(self.tree, 'Species006')

        # decorate the depth of each node on the tree
        decorated_tree = decorate_depth(self.tree)

        # drop duplicate nodes based on depth (deterministic when equal depths)
        obs = drop_duplicate_nodes(decorated_tree, inserted_nodes)

        # verify the resulting tree is correct
        self.assertEqual(obs.getNewick(with_distances=True), EXPECTED_TREE)


TEST_TREE = """(Species003:0.1,(Species001:0.1,Species002:0.1):0.1,((Species006,Species007):0.0,(((Species006,Species007):0.0,Species004:0.1):0.1,((Species006,Species007):0.0,Species005:0.1):0.1):0.1):0.1);"""

EXPECTED_TREE = """(Species003:0.1,(Species001:0.1,Species002:0.1):0.1,((((Species006,Species007):0.0,Species004:0.1):0.1,(Species005:0.1,Species007:0.0):0.1):0.1,Species007:0.0):0.1);"""

if __name__ == "__main__":
    main()
