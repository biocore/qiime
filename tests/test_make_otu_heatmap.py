#!/usr/bin/env python
# file test_make_otu_heatmap_html.py

__author__ = "Dan Knights"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Dan Knights"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Dan Knights"
__email__ = "daniel.knights@colorado.edu"

from os import close
from os.path import exists
from tempfile import mkstemp

from numpy import array, log10, asarray, float64, argwhere
from unittest import TestCase, main
from numpy.testing import assert_almost_equal
from skbio.util.misc import remove_files
from qiime.make_otu_heatmap import (extract_metadata_column,
                                    get_order_from_categories, get_order_from_tree, make_otu_labels,
                                    names_to_indices, get_log_transform, get_clusters,
                                    get_fontsize, plot_heatmap)

from biom.table import Table


class TopLevelTests(TestCase):

    """Tests of top-level functions"""

    def setUp(self):
        """define some top-level data"""

        self.otu_table_values = array([[0, 0, 9, 5, 3, 1],
                                       [1, 5, 4, 0, 3, 2],
                                       [2, 3, 1, 1, 2, 5]])
        {(0, 2): 9.0, (0, 3): 5.0, (0, 4): 3.0, (0, 5): 1.0,
         (1, 0): 1.0, (1, 1): 5.0, (1, 2): 4.0, (1, 4): 3.0, (1, 5): 2.0,
         (2, 0): 2.0, (2, 1): 3.0, (2, 2): 1.0, (2, 3): 1.0, (2, 4): 2.0, (2, 5): 5.0}
        self.otu_table = Table(self.otu_table_values,
                                       ['OTU1', 'OTU2', 'OTU3'],
                                       ['Sample1', 'Sample2', 'Sample3',
                                        'Sample4', 'Sample5', 'Sample6'],
                                       [{"taxonomy": ['Bacteria']},
                                        {"taxonomy": ['Archaea']},
                                        {"taxonomy": ['Streptococcus']}],
                                        [None, None, None, None, None, None])
        self.otu_table_f = Table(self.otu_table_values,
                                         ['OTU1', 'OTU2', 'OTU3'],
                                         ['Sample1', 'Sample2', 'Sample3',
                                          'Sample4', 'Sample5', 'Sample6'],
                                         [{"taxonomy": ['1A', '1B', '1C', 'Bacteria']},
                                          {"taxonomy":
                                           ['2A', '2B', '2C', 'Archaea']},
                                          {"taxonomy": ['3A', '3B', '3C', 'Streptococcus']}],
                                          [None, None, None, None, None, None])

        self.full_lineages = [['1A', '1B', '1C', 'Bacteria'],
                              ['2A', '2B', '2C', 'Archaea'],
                              ['3A', '3B', '3C', 'Streptococcus']]
        self.metadata = [[['Sample1', 'NA', 'A'],
                          ['Sample2', 'NA', 'B'],
                          ['Sample3', 'NA', 'A'],
                          ['Sample4', 'NA', 'B'],
                          ['Sample5', 'NA', 'A'],
                          ['Sample6', 'NA', 'B']],
                         ['SampleID', 'CAT1', 'CAT2'], []]
        self.tree_text = ["('OTU3',('OTU1','OTU2'))"]
        fh, self.tmp_heatmap_fpath = mkstemp(prefix='test_heatmap_',
                                            suffix='.pdf')
        close(fh)

    def test_extract_metadata_column(self):
        """Extracts correct column from mapping file"""
        obs = extract_metadata_column(self.otu_table.sample_ids,
                                      self.metadata, category='CAT2')
        exp = ['A', 'B', 'A', 'B', 'A', 'B']
        self.assertEqual(obs, exp)

    def test_get_order_from_categories(self):
        """Sample indices should be clustered within each category"""
        category_labels = ['A', 'B', 'A', 'B', 'A', 'B']
        obs = get_order_from_categories(self.otu_table, category_labels)
        group_string = "".join([category_labels[i] for i in obs])
        self.assertTrue("AAABBB" == group_string or group_string == "BBBAAA")

    def test_get_order_from_tree(self):
        obs = get_order_from_tree(
            self.otu_table.observation_ids,
            self.tree_text)
        exp = [2, 0, 1]
        assert_almost_equal(obs, exp)

    def test_make_otu_labels(self):
        lineages = []
        for val, id, meta in self.otu_table.iter(axis='observation'):
            lineages.append([v for v in meta['taxonomy']])
        obs = make_otu_labels(self.otu_table.observation_ids,
                              lineages, n_levels=1)
        exp = ['Bacteria (OTU1)', 'Archaea (OTU2)', 'Streptococcus (OTU3)']
        self.assertEqual(obs, exp)

        full_lineages = []
        for val, id, meta in self.otu_table_f.iter(axis='observation'):
            full_lineages.append([v for v in meta['taxonomy']])
        obs = make_otu_labels(self.otu_table_f.observation_ids,
                              full_lineages, n_levels=3)
        exp = ['1B;1C;Bacteria (OTU1)',
               '2B;2C;Archaea (OTU2)',
               '3B;3C;Streptococcus (OTU3)']
        self.assertEqual(obs, exp)

    def test_names_to_indices(self):
        new_order = ['Sample4', 'Sample2', 'Sample3',
                     'Sample6', 'Sample5', 'Sample1']
        obs = names_to_indices(self.otu_table.sample_ids, new_order)
        exp = [3, 1, 2, 5, 4, 0]
        assert_almost_equal(obs, exp)

    def test_get_log_transform(self):
        obs = get_log_transform(self.otu_table)

        data = [val for val in self.otu_table.iter_data(axis='observation')]
        xform = asarray(data, dtype=float64)

        for (i, val) in enumerate(obs.iter_data(axis='observation')):
            non_zeros = argwhere(xform[i] != 0)
            xform[i, non_zeros] = log10(xform[i, non_zeros])
            assert_almost_equal(val, xform[i])

    def test_get_clusters(self):
        data = asarray([val for val in self.otu_table.iter_data(axis='observation')])
        obs = get_clusters(data, axis='row')
        self.assertTrue([0, 1, 2] == obs or obs == [1, 2, 0])
        obs = get_clusters(data, axis='column')
        exp = [2, 3, 1, 4, 0, 5]
        self.assertEqual(obs, exp)

    def test_plot_heatmap(self):
        plot_heatmap(
            self.otu_table, self.otu_table.observation_ids, self.otu_table.sample_ids,
            filename=self.tmp_heatmap_fpath)
        self.assertEqual(exists(self.tmp_heatmap_fpath), True)
        remove_files(set([self.tmp_heatmap_fpath]))


# run tests if called from command line
if __name__ == "__main__":
    main()
