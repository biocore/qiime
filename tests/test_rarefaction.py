#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["justin kuczynski", "Rob Knight", "Jose Carlos Clemente Litran"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Development"

"""Contains tests for producing rarefied OTU tables."""

from cogent.util.unit_test import TestCase, main
import numpy
from qiime.util import get_tmp_filename, load_qiime_config
from qiime.rarefaction import (RarefactionMaker,
                               get_rare_data,
                               remove_empty_otus)
from qiime.format import format_biom_table
from biom.table import table_factory, TableException
from os import remove

class FunctionTests(TestCase):
    def setUp(self):

        self.qiime_config = load_qiime_config()
        self.tmp_dir = self.qiime_config['temp_dir'] or '/tmp/'

        self.otu_table = numpy.array([[2,1,0],
                                      [0,5,0],
                                      [0,3,0],
                                      [1,2,0]])
        self.sample_names = list('YXZ')
        self.taxon_names = list('bacd')

        self.otu_table = table_factory(self.otu_table,
                                        self.sample_names,
                                        self.taxon_names)

        self.otu_table_str = format_biom_table(self.otu_table)

        self.otu_table_fp = get_tmp_filename(tmp_dir=self.tmp_dir,
                                             prefix='test_rarefaction',suffix='.biom')
        open(self.otu_table_fp,'w').write(self.otu_table_str)
        
        self._paths_to_clean_up=[self.otu_table_fp]
    
    def tearDown(self):
        """ cleanup temporary files """
        map(remove,self._paths_to_clean_up)
    
    
    def test_rarefy_to_list(self):
        """rarefy_to_list should rarefy correctly, same names
        
        """
        maker = RarefactionMaker(self.otu_table_fp, 0, 1, 1, 1)
        res = maker.rarefy_to_list(include_full=True)
        self.assertFloatEqual(res[-1][2].SampleIds, self.otu_table.SampleIds)
        self.assertFloatEqual(res[-1][2].ObservationIds, self.otu_table.ObservationIds)
        self.assertEqual(res[-1][2], self.otu_table)

        sample_value_sum = []
        for val in res[1][2].iterSampleData():
            sample_value_sum.append(val.sum())
        self.assertFloatEqual(sample_value_sum, [1.0, 1.0])

    def test_get_empty_rare(self):
        """get_rare_data should be empty when depth > # seqs in any sample"""
        self.assertRaises(TableException, get_rare_data,self.otu_table,
                                           50, include_small_samples=False)
        

    def test_get_overfull_rare(self):
        """get_rare_data should be identical to given in this case

        here, rare depth > any sample, and include_small... = True"""
        rare_otu_table = get_rare_data(self.otu_table,
                                       50, include_small_samples=True)
        self.assertEqual(len(rare_otu_table.SampleIds), 3)
        # 4 observations times 3 samples = size 12 before
        self.assertEqual(len(rare_otu_table.ObservationIds), 4)
        for sam in self.otu_table.SampleIds:
            for otu in self.otu_table.ObservationIds:
                rare_val = rare_otu_table.getValueByIds(otu, sam)
                self.assertEqual(rare_otu_table.getValueByIds(otu, sam),
                                 self.otu_table.getValueByIds(otu, sam))

    def test_get_11depth_rare(self):
        """get_rare_data should get only sample X

        """
        rare_otu_table = get_rare_data(self.otu_table,
                                       11, include_small_samples=False)
        self.assertEqual(rare_otu_table.SampleIds, ['X'])

        # a very complicated way to test things
        rare_values = [val[0] for (val, otu_id, meta) in rare_otu_table.iterObservations()]
        self.assertEqual(rare_values,[1.0, 5.0, 3.0, 2.0])

if __name__ == '__main__':
    main()
