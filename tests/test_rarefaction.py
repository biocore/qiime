#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["justin kuczynski", "Rob Knight",
               "Jose Carlos Clemente Litran", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

"""Contains tests for producing rarefied OTU tables."""

import os
from os import remove, rmdir, close
from shutil import rmtree
from tempfile import mkstemp, mkdtemp
from unittest import TestCase, main

import numpy.testing as npt
import numpy as np
from biom.table import Table, TableException
from biom.parse import parse_biom_table
from biom.util import biom_open

from qiime.rarefaction import RarefactionMaker, get_rare_data
from qiime.util import load_qiime_config, write_biom_table


class FunctionTests(TestCase):

    def setUp(self):

        self.qiime_config = load_qiime_config()
        self.tmp_dir = self.qiime_config['temp_dir'] or '/tmp/'

        self.otu_table_data = np.array([[2, 1, 0],
                                        [0, 5, 0],
                                        [0, 3, 0],
                                        [1, 2, 0]])
        self.sample_names = list('YXZ')
        self.taxon_names = list('bacd')
        self.otu_metadata = [{'domain': 'Archaea'},
                             {'domain': 'Bacteria'},
                             {'domain': 'Bacteria'},
                             {'domain': 'Bacteria'}]

        self.otu_table = Table(self.otu_table_data,
                               self.taxon_names,
                               self.sample_names)

        self.otu_table_meta = Table(self.otu_table_data,
                                    self.taxon_names, self.sample_names,
                                    observation_metadata=self.otu_metadata)

        fd, self.otu_table_fp = mkstemp(dir=self.tmp_dir,
                                        prefix='test_rarefaction',
                                        suffix='.biom')
        close(fd)
        fd, self.otu_table_meta_fp = mkstemp(dir=self.tmp_dir,
                                             prefix='test_rarefaction',
                                             suffix='.biom')
        close(fd)

        self.rare_dir = mkdtemp(dir=self.tmp_dir,
                                prefix='test_rarefaction_dir', suffix='')

        write_biom_table(self.otu_table, self.otu_table_fp)
        write_biom_table(self.otu_table_meta, self.otu_table_meta_fp)

        self._paths_to_clean_up = [self.otu_table_fp, self.otu_table_meta_fp]
        self._dirs_to_clean_up = [self.rare_dir]

    def tearDown(self):
        """ cleanup temporary files """
        map(remove, self._paths_to_clean_up)
        for d in self._dirs_to_clean_up:
            if os.path.exists(d):
                rmtree(d)

    def test_rarefy_to_list(self):
        """rarefy_to_list should rarefy correctly, same names

        """
        maker = RarefactionMaker(self.otu_table_fp, 0, 1, 1, 1)
        res = maker.rarefy_to_list(include_full=True)
        self.assertItemsEqual(res[-1][2].sample_ids, self.otu_table.sample_ids)
        self.assertItemsEqual(
            res[-1][2].observation_ids,
            self.otu_table.observation_ids)
        self.assertEqual(res[-1][2], self.otu_table)

        sample_value_sum = []
        for val in res[1][2].iter_data(axis='sample'):
            sample_value_sum.append(val.sum())
        npt.assert_almost_equal(sample_value_sum, [1.0, 1.0])

    def test_rarefy_to_files(self):
        """rarefy_to_files should write valid files

        """
        maker = RarefactionMaker(self.otu_table_fp, 0, 1, 1, 1)
        maker.rarefy_to_files(
            self.rare_dir,
            include_full=True,
            include_lineages=False)

        fname = os.path.join(self.rare_dir, "rarefaction_1_0.biom")
        with biom_open(fname, 'U') as biom_file:
            otu_table = Table.from_hdf5(biom_file)

        self.assertItemsEqual(
            otu_table.sample_ids,
            self.otu_table.sample_ids[:2])
        # third sample had 0 seqs, so it's gone

    def test_rarefy_to_files2(self):
        """rarefy_to_files should write valid files with some metadata on otus

        """
        maker = RarefactionMaker(self.otu_table_meta_fp, 0, 1, 1, 1)
        maker.rarefy_to_files(
            self.rare_dir,
            include_full=True,
            include_lineages=False)

        fname = os.path.join(self.rare_dir, "rarefaction_1_0.biom")
        with biom_open(fname, 'U') as biom_file:
            otu_table = Table.from_hdf5(biom_file)

        self.assertItemsEqual(
            otu_table.sample_ids,
            self.otu_table.sample_ids[:2])
        # third sample had 0 seqs, so it's gone

    def test_get_empty_rare(self):
        """get_rare_data should be empty when depth > # seqs in any sample"""
        self.assertRaises(TableException, get_rare_data, self.otu_table,
                          50, include_small_samples=False)

    def test_get_overfull_rare(self):
        """get_rare_data should be identical to given in this case

        here, rare depth > any sample, and include_small... = True"""
        rare_otu_table = get_rare_data(self.otu_table,
                                       50, include_small_samples=True)
        self.assertEqual(len(rare_otu_table.sample_ids), 3)
        # 4 observations times 3 samples = size 12 before
        self.assertEqual(len(rare_otu_table.observation_ids), 4)
        for sam in self.otu_table.sample_ids:
            for otu in self.otu_table.observation_ids:
                rare_val = rare_otu_table.get_value_by_ids(otu, sam)
                self.assertEqual(rare_otu_table.get_value_by_ids(otu, sam),
                                 self.otu_table.get_value_by_ids(otu, sam))

    def test_get_11depth_rare(self):
        """get_rare_data should get only sample X

        """
        rare_otu_table = get_rare_data(self.otu_table,
                                       11, include_small_samples=False)
        self.assertEqual(rare_otu_table.sample_ids, ('X',))

        # a very complicated way to test things
        rare_values = [val[0]
                       for (val, otu_id, meta) in rare_otu_table.iter(axis='observation')]
        self.assertEqual(rare_values, [1.0, 5.0, 3.0, 2.0])

if __name__ == '__main__':
    main()
