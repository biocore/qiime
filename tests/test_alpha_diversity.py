#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project" #consider project name
__credits__ = ["Justin Kuczynski", "Rob Knight"] #remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.8.0"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Development"

"""Contains tests for performing alpha diversity analyses within each sample."""

from numpy import array
import numpy
from shutil import rmtree
from os import makedirs
from os.path import exists
from cogent.util.unit_test import TestCase, main
from cogent.maths.unifrac.fast_unifrac import PD_whole_tree
#from cogent.maths.stats.alpha_diversity import (observed_species, osd)
from qiime.pycogent_backports.alpha_diversity import (observed_species, osd)

from cogent.util.misc import remove_files
from qiime.util import get_tmp_filename, load_qiime_config
from qiime.alpha_diversity import AlphaDiversityCalc, AlphaDiversityCalcs
import qiime.alpha_diversity
from qiime.parse import parse_newick
from qiime.format import format_biom_table
from biom.table import table_factory, DenseOTUTable

class AlphaDiversitySharedSetUpTests(TestCase):

    def setUp(self):
        """Define some test data."""
        self.qiime_config = load_qiime_config()
        self.dirs_to_remove = []
        
        self.tmp_dir = self.qiime_config['temp_dir'] or '/tmp/'
        if not exists(self.tmp_dir):
            makedirs(self.tmp_dir)
            # if test creates the temp dir, also remove it
            self.dirs_to_remove.append(self.tmp_dir)
        
        self.otu_table1 = table_factory(data=array([[2,0,0,1],
                                                   [1,1,1,1],
                                                   [0,0,0,0]]).T,
                                       sample_ids=list('XYZ'),
                                       observation_ids=list('abcd'),
                                       constructor=DenseOTUTable)
        self.otu_table1_fp = get_tmp_filename(tmp_dir=self.tmp_dir,
                                             prefix='alpha_diversity_tests',
                                             suffix='.biom',
                                             result_constructor=str)
        open(self.otu_table1_fp,'w').write(\
         format_biom_table(self.otu_table1))

        self.otu_table2 = table_factory(data=array([[2,0,0,1],
                                                   [1,1,1,1],
                                                   [0,0,0,0]]).T,
                                       sample_ids=list('XYZ'),
                                       observation_ids=['a','b','c','d_'],
                                       constructor=DenseOTUTable)
        self.otu_table2_fp = get_tmp_filename(tmp_dir=self.tmp_dir,
                                             prefix='alpha_diversity_tests',
                                             suffix='.biom',
                                             result_constructor=str)
        open(self.otu_table2_fp,'w').write(\
         format_biom_table(self.otu_table2))
        
        self.single_sample_otu_table = table_factory(data=array([[2,0,0,1]]).T,
                                                     sample_ids=list('X'),
                                                     observation_ids=list('abcd'),
                                                     constructor=DenseOTUTable)
        self.single_sample_otu_table_fp = get_tmp_filename(tmp_dir=self.tmp_dir,
                                             prefix='alpha_diversity_tests',
                                             suffix='.biom',
                                             result_constructor=str)
        open(self.single_sample_otu_table_fp,'w').write(\
         format_biom_table(self.single_sample_otu_table))
        
        
        self.tree1 = parse_newick('((a:2,b:3):2,(c:1,d:2):7);')
        self.tree2 = parse_newick("((a:2,'b':3):2,(c:1,'d_':2):7);")
        
        self.files_to_remove = [self.otu_table1_fp,self.otu_table2_fp,
                                self.single_sample_otu_table_fp]
        
    def tearDown(self):
        """ """
        remove_files(self.files_to_remove)
        # remove directories last, so we don't get errors
        # trying to remove files which may be in the directories
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)

class AlphaDiversityCalcTests(AlphaDiversitySharedSetUpTests):
    """Tests of the AlphaDiversityCalc class"""

    def test_init(self):
        """AlphaDiversity __init__ should store metric, name, params"""
        c = AlphaDiversityCalc(observed_species)
        self.assertEqual(c.Metric, observed_species)
        self.assertEqual(c.Params, {})

    def test_call(self):
        """AlphaDiversityCalc __call__ should call metric on data
        and return correct result"""
        c = AlphaDiversityCalc(observed_species)
        self.assertEqual(c(data_path=self.otu_table1_fp), [2,4,0])
    
    def test_multi_return(self):
        """AlphaDiversityCalc __call__ should call metric on data
        and return correct result for metric fn that returns a len 3 tuple
        """
        c = AlphaDiversityCalc(osd)
        res = c(data_path=self.otu_table1_fp)
        self.assertEqual(res, array([[2,1,1],
                                    [4,4,0],
                                    [0,0,0]]))

    def test_1sample(self):
        """ should work if only testing one sample as well"""
        c = AlphaDiversityCalc(qiime.alpha_diversity.alph.observed_species)
        self.assertEqual(c(data_path=self.single_sample_otu_table_fp), [2])
 
    def test_call_phylogenetic(self):
        """AlphaDiversityCalc __call__ should call metric on phylo data
        and return correct values"""
        c = AlphaDiversityCalc(metric=PD_whole_tree,
            is_phylogenetic=True)
        self.assertEqual(c(data_path=self.otu_table1_fp, tree_path=self.tree1, \
            taxon_names=self.otu_table1.ObservationIds, 
            sample_names=self.otu_table1.SampleIds), 
            [13, 17, 0])
            
    def test_call_phylogenetic_escaped_names(self):
        """AlphaDiversityCalc __call__ should call metric on phylo data
        and return correct values"""
        
        c = AlphaDiversityCalc(metric=PD_whole_tree,is_phylogenetic=True)
        expected = [13., 17., 0.]
        
        non_escaped_result = c(data_path=self.otu_table1_fp, 
                               tree_path=self.tree1,
                               taxon_names = self.otu_table1.ObservationIds,
                               sample_names=self.otu_table1.SampleIds)
        
        escaped_result  = c(data_path=self.otu_table2_fp,
                            tree_path=self.tree2,
                            taxon_names = self.otu_table2.ObservationIds,
                            sample_names=self.otu_table2.SampleIds)
        
        self.assertEqual(non_escaped_result,expected)
        self.assertEqual(escaped_result,expected)
        self.assertEqual(non_escaped_result,escaped_result)

class AlphaDiversityCalcsTests(AlphaDiversitySharedSetUpTests):
    """Tests of the AlphaDiversityCalcs class"""

    def test1(self):
        """ checks that output from AlphaDiversityCalcs is the right shape 
        when run on phylo, multiple return value nonphylo, and another nonphylo
        """
        calc1 = AlphaDiversityCalc(metric=observed_species)
        calc2 = AlphaDiversityCalc(metric=PD_whole_tree,
            is_phylogenetic=True)
        calc3 = AlphaDiversityCalc(metric=osd)
        calcs = AlphaDiversityCalcs([calc1, calc2, calc3])
        results = calcs(data_path=self.otu_table1_fp, 
            tree_path=self.tree1,
            result_path=None, log_path=None)
        self.assertEqual(results[0].shape, (3,5))
        self.assertEqual(len(results[1]), 3)
        self.assertEqual(len(results[2]), 5)

#run tests if called from command line
if __name__ == '__main__':
    main()
