#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2010, The QIIME Project" #consider project name
__credits__ = ["Justin Kuczynski", "Rob Knight"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Release"

"""Contains tests for performing alpha diversity analyses within each sample."""

from numpy import array
import numpy
from cogent.util.unit_test import TestCase, main
from cogent.maths.unifrac.fast_unifrac import PD_whole_tree
from cogent.maths.stats.alpha_diversity import (observed_species, osd)
from cogent.util.misc import remove_files
from cogent.app.util import get_tmp_filename
from qiime.alpha_diversity import AlphaDiversityCalc, AlphaDiversityCalcs
import qiime.alpha_diversity as alph454
from qiime.parse import parse_newick

class AlphaDiversityCalcTests(TestCase):
    """Tests of the AlphaDiversityCalc class"""

    def setUp(self):
        """Define some test data."""
        self.otu_table = array([[2,0,0,1],
                                [1,1,1,1],
                                [0,0,0,0]])
        self.sample_names = list('XYZ')
        self.otu_names = list('abcd')
        self.otu_tuple = (self.sample_names, self.otu_names, self.otu_table.T,
        None)
        self.tree = parse_newick('((a:2,b:3):2,(c:1,d:2):7);')
        
        self.files_to_remove = []
        
    def tearDown(self):
        remove_files(self.files_to_remove)

    def test_init(self):
        """AlphaDiversity __init__ should store metric, name, params"""
        c = AlphaDiversityCalc(observed_species)
        self.assertEqual(c.Metric, observed_species)
        self.assertEqual(c.Params, {})

    def test_call(self):
        """AlphaDiversityCalc __call__ should call metric on data
        and return correct result"""
        c = AlphaDiversityCalc(observed_species)
        self.assertEqual(c(data_path=self.otu_table), [2,4,0])
    
    def test_multi_return(self):
        """AlphaDiversityCalc __call__ should call metric on data
        and return correct result for metric fn that returns a len 3 tuple
        """
        c = AlphaDiversityCalc(osd)
        res = c(data_path=self.otu_table)
        self.assertEqual(res, array([[2,1,1],
                                    [4,4,0],
                                    [0,0,0]]))

    def test_1sample(self):
        """ should work if only testing one sample as well"""
        otu_table = array([[2,0,0,1]])
        c = AlphaDiversityCalc(alph454.alph.observed_species)
        self.assertEqual(c(data_path=otu_table), [2])
 
    def test_call_phylogenetic(self):
        """AlphaDiversityCalc __call__ should call metric on phylo data
        and return correct values"""
        c = AlphaDiversityCalc(metric=PD_whole_tree,
            is_phylogenetic=True)
        self.assertEqual(c(data_path=self.otu_table, tree_path=self.tree, \
            taxon_names = self.otu_names, sample_names=self.sample_names), 
            [13, 17, 0])
            
    def test_call_phylogenetic_escaped_names(self):
        """AlphaDiversityCalc __call__ should call metric on phylo data
        and return correct values"""
        c = AlphaDiversityCalc(metric=PD_whole_tree,
            is_phylogenetic=True)
        expected = [13, 17, 0]
        non_escaped_result = c(data_path=self.otu_table, tree_path=self.tree, \
            taxon_names = self.otu_names, sample_names=self.sample_names)
        
        
        otu_table = array([[2,0,0,1],
                           [1,1,1,1],
                           [0,0,0,0]])
        sample_names = list('XYZ')
        otu_names = ['a','b','c','d_']
        tree_str = "((a:2,'b':3):2,(c:1,'d_':2):7);"
        tree_fp = get_tmp_filename(prefix='Alpha_div_tests',suffix='.tre')
        open(tree_fp,'w').write(tree_str)
        self.files_to_remove.append(tree_fp)
        escaped_result  = c(data_path=otu_table, tree_path=tree_fp, \
            taxon_names = otu_names, sample_names=sample_names)
            
        self.assertEqual(non_escaped_result,escaped_result)

class AlphaDiversityCalcsTests(TestCase):
    """Tests of the AlphaDiversityCalcs class"""

    def setUp(self):
        """Define some test data."""
        self.otu_table = array([[2,0,0,1],
                                [1,1,1,1],
                                [0,0,0,0]])
        self.sample_names = list('XYZ')
        self.otu_names = list('abcd')
        self.otu_tuple = (self.sample_names, self.otu_names, self.otu_table.T,
        None)
        self.tree = parse_newick('((a:2,b:3):2,(c:1,d:2):7);')

    def test1(self):
        """ checks that output from AlphaDiversityCalcs is the right shape 
        when run on phylo, multiple return value nonphylo, and another nonphylo
        """
        calc1 = AlphaDiversityCalc(metric=observed_species)
        calc2 = AlphaDiversityCalc(metric=PD_whole_tree,
            is_phylogenetic=True)
        calc3 = AlphaDiversityCalc(metric=osd)
        calcs = AlphaDiversityCalcs([calc1, calc2, calc3])
        results = calcs(data_path=self.otu_tuple, 
            tree_path=self.tree,
            result_path=None, log_path=None)
        self.assertEqual(results[0].shape, (3,5))
        self.assertEqual(len(results[1]), 3)
        self.assertEqual(len(results[2]), 5)

#run tests if called from command line
if __name__ == '__main__':
    main()
