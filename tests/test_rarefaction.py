#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["justin kuczynski", "Rob Knight", "Jose Carlos Clemente Litran"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Development"

"""Contains tests for producing rarefied OTU tables."""

from qiime.rarefaction import RarefactionMaker, get_rare_data, remove_empty_otus
from cogent.util.unit_test import TestCase, main
import numpy
from qiime.pycogent_backports.rich_otu_table import SparseOTUTable, to_ll_mat, TableException
from qiime.util import get_tmp_filename, load_qiime_config

class FunctionTests(TestCase):
    def setUp(self):
        # Since Justin's original OTU table had one column with all
        #  zeros, I have to specificy here one entry in the sparse
        #  representation in that column even if it's zero, just
        #  so to_ll_mat knows of the existence of that column

        self.qiime_config = load_qiime_config()
        self.tmp_dir = self.qiime_config['temp_dir'] or '/tmp/'

        self.otu_table_values = {(0, 0):2.0, (0, 1):1.0, (0,2):0.0,
                                 (1, 1):5.0,
                                 (2, 1):3.0,
                                 (3, 0):1.0, (3, 1):2.0}

        self.otu_table = SparseOTUTable(to_ll_mat(self.otu_table_values),
                                        ['Y', 'X', 'Z'],
                                        ['b', 'a', 'c', 'd'],
                                        [None, None, None],
                                        [None, None, None, None])

        self.otu_table_str = self.otu_table.getBiomFormatJsonString()

        self.otu_table_fp = get_tmp_filename(tmp_dir=self.tmp_dir,
                                             prefix='test_rarefaction',suffix='.biom')
        open(self.otu_table_fp,'w').write(self.otu_table_str)

        #self.otu_table_transpose = numpy.array([
        #    [2,0,0,1],
        #    [1,5,3,2],
        #    [0,0,0,0],
        #    ])
        #self.otu_table = self.otu_table_transpose.T
        #self.sample_names = list('YXZ')
        #self.taxon_names = list('bacd')
        #self.otu_tuple = (self.sample_names, self.taxon_names, 
        #    self.otu_table_transpose.T, None)
    
    def test_rarefy_to_list(self):
        """rarefy_to_list should rarefy correctly, same names, rm empty samples
        
        """
        #maker = RarefactionMaker(self.otu_tuple, 0, 1, 1, 1)
        maker = RarefactionMaker(self.otu_table_fp, 0, 1, 1, 1)
        res = maker.rarefy_to_list(include_full=True)
        self.assertFloatEqual(res[-1][2].SampleIds, self.otu_table.SampleIds)
        self.assertFloatEqual(res[-1][2].ObservationIds, self.otu_table.ObservationIds)
        #self.assertFloatEqual(res[-1][2], self.otu_table_transpose.T)
        self.assertEqual(res[-1][2], self.otu_table)
        
        # each sample should have 1 seq, sample z should be removed
        #self.assertFloatEqual((res[1][4]).sum(0),[1.0,1.0] )

        sample_value_sum = []
        for val in res[1][2].iterSampleData():
            sample_value_sum.append(val.sum())
        self.assertFloatEqual(sample_value_sum, [1.0, 1.0])

    def test_get_empty_rare(self):
        """get_rare_data should be empty when depth > # seqs in any sample"""
        #rare_sample_ids, rare_otu_table = get_rare_data(
        #    self.sample_names, self.otu_table, \
        #    50, include_small_samples=False)
        self.assertRaises(TableException, get_rare_data,self.otu_table,
                                           50, include_small_samples=False)

        #self.assertEqual(len(rare_otu_table.SampleIds), 0)
        #self.assertEqual(rare_otu_table.size, 0)

        #self.assertEqual(rare_otu_table.isEmpty,True)
        

    def test_get_overfull_rare(self):
        """get_rare_data should be identical to given in this case

        here, rare depth > any sample, and include_small... = True"""
        #rare_sample_ids, rare_otu_table = get_rare_data(
        #    self.sample_names, self.otu_table, \
        #    50, include_small_samples=True)
        rare_otu_table = get_rare_data(self.otu_table,
                                       50, include_small_samples=True)
        #self.assertEqual(len(rare_sample_ids), 3)
        #self.assertEqual(rare_otu_table.size, 12)
        self.assertEqual(len(rare_otu_table.SampleIds), 3)
        # 4 observations times 3 samples = size 12 before
        self.assertEqual(len(rare_otu_table.ObservationIds), 4)
        #for i, sam in enumerate(self.sample_names):
        #    for j, otu in enumerate(self.taxon_names):
        #        rare_val = rare_otu_table[self.taxon_names.index(otu),
        #            rare_sample_ids.index(sam)]
        #        self.assertEqual(rare_val, self.otu_table[j,i]) 
        for sam in self.otu_table.SampleIds:
            for otu in self.otu_table.ObservationIds:
                rare_val = rare_otu_table.getValueByIds(otu, sam)
                self.assertEqual(rare_otu_table.getValueByIds(otu, sam),
                                 self.otu_table.getValueByIds(otu, sam))

    def test_get_11depth_rare(self):
        """get_rare_data should get only sample X

        """
        #rare_sample_ids, rare_otu_table = get_rare_data(
        #    self.sample_names, self.otu_table, \
        #    11, include_small_samples=False)
        rare_otu_table = get_rare_data(self.otu_table,
                                       11, include_small_samples=False)
        #self.assertEqual(rare_sample_ids, ['X'])
        self.assertEqual(rare_otu_table.SampleIds, ['X'])
        #rare_otu_table[numpy.argsort(rare_otu_ids)]
        #self.assertEqual(rare_otu_table[numpy.argsort(self.taxon_names)][:,0], 
        #    numpy.array([5,1,3,2]))

        # a very complicated way to test things
        #self.assertEqual(rare_otu_table[numpy.argsort(self.taxon_names)][:,0], 
        #                 numpy.array([5,1,3,2]))
        rare_values = [val[0] for (val, otu_id, meta) in rare_otu_table.iterObservations()]
        self.assertEqual(rare_values,[1.0, 5.0, 3.0, 2.0])

if __name__ == '__main__':
    main()
