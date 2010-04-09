#!/usr/bin/env python

"""Tests of code for building tree from aligned sequences"""

__author__ = "Rob Knight"
__copyright__ = "Copyright 2010, The QIIME Project" 
#remember to add yourself if you make changes
__credits__ = ["Rob Knight"] 
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Rob Knight"
__email__ = "rob.knight@colorado.edu"
__status__ = "Release"

from os import remove
from cogent import LoadSeqs, DNA
from cogent.util.unit_test import TestCase, main
from cogent.app.util import get_tmp_filename
import cogent.app.fasttree
from qiime.make_phylogeny import TreeBuilder, CogentTreeBuilder

def safe_remove(f):
    try:
        remove(f)
    except OSError:
        pass

class TreeBuilderTests(TestCase):
    """Tests of the abstract TreeBuilder class"""

    def test_init(self):
        """Abstract TreeBuilder __init__ should store name, params"""
        p = TreeBuilder({})
        self.assertEqual(p.Name, 'TreeBuilder')
        self.assertEqual(p.Params, {})

    def test_call(self):
        """Abstract TreeBuilder __call__ should raise NotImplementedError"""
        p = TreeBuilder({})
        self.assertRaises(NotImplementedError, p, '/path/to/seqs')

class SharedSetupTestCase(TestCase):
    """Shared setup for aligner tests"""
       
    def tearDown(self):
        map(safe_remove,self._paths_to_clean_up)
 
class CogentTreeBuilderTests(SharedSetupTestCase):
    """Tests of the CogentTreeBuilder class"""
    def setUp(self):
        self.input_fp = get_tmp_filename(\
         prefix='CogentTreeBuilderTests_',suffix='.fasta')
        self._paths_to_clean_up =\
         [self.input_fp] 
        open(self.input_fp,'w').write(aln_for_tree)

    def test_call_correct_alignment(self):
        """CogentTreeBuilder: output expected alignment file
        """
        p = CogentTreeBuilder({'Module': cogent.app.fasttree})
        log_fp = get_tmp_filename(\
         prefix='CogentTreeBuilderTests_',suffix='.log')
        self._paths_to_clean_up.append(log_fp)
         
        actual = p(result_path=None, aln_path=self.input_fp,
            log_path=log_fp)
        expected = tree
        #note: lines in diff order w/ diff versions
        self.assertEqual(str(actual),expected)

    # this test assumes a certain newick format.  always works for me, but isn't
    # a good test
    # def test_midpoint_rooting(self):
    #     """CogentTreeBuilder: midpoint rooting should work"""
    #     p = CogentTreeBuilder({'Module': cogent.app.fasttree})
    #     log_fp = get_tmp_filename(\
    #      prefix='CogentTreeBuilderTests_',suffix='.log')
    #     self._paths_to_clean_up.append(log_fp)
    #      
    #     actual = p(result_path=None, aln_path=self.input_fp,
    #         log_path=log_fp,root_method='midpoint')
    #     expected = midpoint_tree
    #     #note: lines in diff order w/ diff versions
    #     self.assertEqual(str(actual),expected)

    def test_root_midpt(self):
        """midpoint should be selected correctly when it is an internal node
        """
        from qiime.make_phylogeny import root_midpt
        from cogent.core.tree import PhyloNode
        from cogent.parse.tree import DndParser
        tree = DndParser('(a:1,((c:1,d:3)n3:1,b:1)n2:1)rt;')
        tmid = root_midpt(tree)
        self.assertEqual(tmid.getDistances(),tree.getDistances())
        tipnames = tree.getTipNames()
        nontipnames = [t.Name for t in tree.nontips()]
        # for tipname in tipnames:
        #     tmid_tip = tmid.getNodeMatchingName(tipname)
        #     orig_tip = tree.getNodeMatchingName(tipname)
        #     for nontipname in nontipnames:
        #         tmid_dist=\
        #           tmid.getNodeMatchingName(nontipname).distance(tmid_tip)
        #         orig_dist=\
        #           tree.getNodeMatchingName(nontipname).distance(orig_tip)
        #         print nontipname, tipname, 'assert'
                # self.assertEqual(tmid_dist, orig_dist)
        self.assertTrue(tmid.isRoot())
        self.assertEqual(tmid.distance(\
            tmid.getNodeMatchingName('d')), 3)
                    
                    
                
                    
    def test_root_midpt2(self):
        """ midpoint between nodes should behave correctly"""
        from qiime.make_phylogeny import root_midpt
        from cogent.core.tree import PhyloNode
        from cogent.parse.tree import DndParser
        tree = DndParser('(a:1,((c:1,d:2.5)n3:1,b:1)n2:1)rt;')
        tmid = root_midpt(tree)
        self.assertEqual(tmid.getDistances(),tree.getDistances())
        tipnames = tree.getTipNames()
        nontipnames = [t.Name for t in tree.nontips()]
        self.assertTrue(tmid.isRoot())
        self.assertEqual(tmid.distance(tmid.getNodeMatchingName('d')), 2.75)


from copy import copy, deepcopy
from cogent import LoadTree
from cogent.core.tree import TreeNode, PhyloNode, TreeError
from cogent.parse.tree import DndParser
from cogent.maths.stats.test import correlation
from cogent.util.unit_test import TestCase, main
from numpy import array, arange        
class PhyloNodeTests(TestCase):
    """this tests only maxtiptipdistance"""

    def setUp(self):
        """Creates a standard tree"""
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

        self.s = '((H:1,G:1):2,(R:0.5,M:0.7):3);'
        self.t = DndParser(self.s,PhyloNode)
        self.s3 = '(((H:1,G:1,O:1):2,R:3):1,X:4);'
        self.t3 = DndParser(self.s3, PhyloNode)
        
    def test_maxTipTipDistance(self):
        """maxTipTipDistance returns the max dist between any pair of tips"""
        nodes, tree = self.TreeNode, self.TreeRoot
        max_dist, tip_pair = tree.maxTipTipDistance()
        self.assertEqual(max_dist, 10)
        try:
            self.assertEqual(tip_pair, ('h', 'g'))
        except AssertionError:
            self.assertEqual(tip_pair, ('g', 'h'))
aln_for_tree = """>jkl\n--TTACAC--\n>abc\nACACACAC--\n>ghi\nACAGACACTT\n>def\nACAGACAC--\n"""

tree = '(def:0.00014,ghi:0.00014,(abc:0.07248,jkl:0.40293)0.742:0.07156);'
midpoint_tree = '(jkl:0.237705,(abc:0.07248,(def:0.00014,ghi:0.00014)0.742:0.07156):0.165225);'

#run unit tests if run from command-line
if __name__ == '__main__':
    main()
