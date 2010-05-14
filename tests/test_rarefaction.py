#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["justin kuczynski", "Rob Knight"]
__license__ = "GPL"
__version__ = "1.1.0"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Release"

"""Contains tests for producing rarefied OTU tables."""

from qiime.rarefaction import RarefactionMaker, get_rare_data
from cogent.util.unit_test import TestCase, main
import numpy

class FunctionTests(TestCase):
    def setUp(self):
        self.otu_table_transpose = numpy.array([
            [2,0,0,1],
            [1,5,3,2],
            [0,0,0,0],
            ])
        self.otu_table = self.otu_table_transpose.T
        self.sample_names = list('YXZ')
        self.taxon_names = list('bacd')
        self.otu_tuple = (self.sample_names, self.taxon_names, 
            self.otu_table_transpose.T, None)
    
    def test_rarefy_to_list(self):
        """rarefy_to_list should rarefy correctly, same names, rm empty samples
        
        """
        maker = RarefactionMaker(self.otu_tuple, 0, 1, 1, 1)
        res = maker.rarefy_to_list(include_full=True)
        self.assertFloatEqual(res[-1][2], self.sample_names)
        self.assertFloatEqual(res[-1][3], self.taxon_names)
        self.assertFloatEqual(res[-1][4], self.otu_table_transpose.T)
        
        # each sample should have 1 seq, sample z should be removed
        self.assertFloatEqual((res[1][4]).sum(0),[1.0,1.0] )

    def test_get_empty_rare(self):
        """get_rare_data should be empty when depth > # seqs in any sample"""
        rare_sample_ids, rare_otu_table = get_rare_data(
            self.sample_names, self.otu_table, \
            50, include_small_samples=False)
        self.assertEqual(len(rare_sample_ids), 0)
        self.assertEqual(rare_otu_table.size, 0)    

    def test_get_overfull_rare(self):
        """get_rare_data should be identical to given in this case

        here, rare depth > any sample, and include_small... = True"""
        rare_sample_ids, rare_otu_table = get_rare_data(
            self.sample_names, self.otu_table, \
            50, include_small_samples=True)
        self.assertEqual(len(rare_sample_ids), 3)
        self.assertEqual(rare_otu_table.size, 12)
        for i, sam in enumerate(self.sample_names):
            for j, otu in enumerate(self.taxon_names):
                rare_val = rare_otu_table[self.taxon_names.index(otu),
                    rare_sample_ids.index(sam)]
                self.assertEqual(rare_val, self.otu_table[j,i]) 

    def test_get_11depth_rare(self):
        """get_rare_data should get only sample X

        """
        rare_sample_ids, rare_otu_table = get_rare_data(
            self.sample_names, self.otu_table, \
            11, include_small_samples=False)
        self.assertEqual(rare_sample_ids, ['X'])
        #rare_otu_table[numpy.argsort(rare_otu_ids)]
        self.assertEqual(rare_otu_table[numpy.argsort(self.taxon_names)][:,0], 
            numpy.array([5,1,3,2]))

        
#run tests if called from command line
if __name__ == '__main__':
    main()
