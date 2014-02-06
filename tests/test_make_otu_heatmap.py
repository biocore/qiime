#!/usr/bin/env python
# file test_make_otu_heatmap_html.py

__author__ = "Dan Knights"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Dan Knights"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Dan Knights"
__email__ = "daniel.knights@colorado.edu"

from os.path import exists
from numpy import array, log10, asarray, float64
from cogent.util.unit_test import TestCase, main
from cogent.util.misc import remove_files
from qiime.make_otu_heatmap import (extract_metadata_column,
                                    get_order_from_categories, get_order_from_tree, make_otu_labels,
                                    names_to_indices, get_log_transform, get_clusters,
                                    get_fontsize, plot_heatmap)
from qiime.util import get_tmp_filename
from biom.table import table_factory


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
        self.otu_table = table_factory(self.otu_table_values,
                                       ['Sample1', 'Sample2', 'Sample3',
                                        'Sample4', 'Sample5', 'Sample6'],
                                       ['OTU1', 'OTU2', 'OTU3'],
                                       [None, None, None, None, None, None],
                                       [{"taxonomy": ['Bacteria']},
                                        {"taxonomy": ['Archaea']},
                                        {"taxonomy": ['Streptococcus']}])
        self.otu_table_f = table_factory(self.otu_table_values,
                                         ['Sample1', 'Sample2', 'Sample3',
                                          'Sample4', 'Sample5', 'Sample6'],
                                         ['OTU1', 'OTU2', 'OTU3'],
                                         [None, None, None, None, None, None],
                                         [{"taxonomy": ['1A', '1B', '1C', 'Bacteria']},
                                          {"taxonomy":
                                           ['2A', '2B', '2C', 'Archaea']},
                                          {"taxonomy": ['3A', '3B', '3C', 'Streptococcus']}])

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
        self.tmp_heatmap_fpath = get_tmp_filename(
            prefix='test_heatmap_',
            suffix='.pdf'
        )

    def test_extract_metadata_column(self):
        """Extracts correct column from mapping file"""
        obs = extract_metadata_column(self.otu_table.SampleIds,
                                      self.metadata, category='CAT2')
        exp = ['A', 'B', 'A', 'B', 'A', 'B']
        self.assertEqual(obs, exp)

    def test_get_order_from_categories(self):
        """Sample indices should be clustered within each category"""
        category_labels = ['A', 'B', 'A', 'B', 'A', 'B']
        obs = get_order_from_categories(self.otu_table, category_labels)
        exp = [0, 4, 2, 1, 5, 3]
        self.assertEqual(obs, exp)

    def test_get_order_from_tree(self):
        obs = get_order_from_tree(
            self.otu_table.ObservationIds,
            self.tree_text)
        exp = [2, 0, 1]
        self.assertEqual(obs, exp)

    def test_make_otu_labels(self):
        lineages = []
        for val, id, meta in self.otu_table.iterObservations():
            lineages.append([v for v in meta['taxonomy']])
        obs = make_otu_labels(self.otu_table.ObservationIds,
                              lineages, n_levels=1)
        exp = ['Bacteria (OTU1)', 'Archaea (OTU2)', 'Streptococcus (OTU3)']
        self.assertEqual(obs, exp)

        full_lineages = []
        for val, id, meta in self.otu_table_f.iterObservations():
            full_lineages.append([v for v in meta['taxonomy']])
        obs = make_otu_labels(self.otu_table_f.ObservationIds,
                              full_lineages, n_levels=3)
        exp = ['1B;1C;Bacteria (OTU1)',
               '2B;2C;Archaea (OTU2)',
               '3B;3C;Streptococcus (OTU3)']
        self.assertEqual(obs, exp)

    def test_names_to_indices(self):
        new_order = ['Sample4', 'Sample2', 'Sample3',
                     'Sample6', 'Sample5', 'Sample1']
        obs = names_to_indices(self.otu_table.SampleIds, new_order)
        exp = [3, 1, 2, 5, 4, 0]
        self.assertEqual(obs, exp)

    def test_get_log_transform(self):
        eps = .01
        obs = get_log_transform(self.otu_table, eps=eps)

        data = [val for val in self.otu_table.iterObservationData()]
        xform = asarray(data, dtype=float64)
        xform[xform == 0] = eps

        for (i, val) in enumerate(obs.iterObservationData()):
            self.assertEqual(val, log10(xform[i]))

    def test_get_clusters(self):
        data = asarray([val for val in self.otu_table.iterObservationData()])
        obs = get_clusters(data, axis='row')
        exp = [0, 1, 2]
        self.assertEqual(obs, exp)
        obs = get_clusters(data, axis='column')
        exp = [0, 5, 4, 1, 2, 3]
        self.assertEqual(obs, exp)

    def test_plot_heatmap(self):
        plot_heatmap(
            self.otu_table, self.otu_table.ObservationIds, self.otu_table.SampleIds,
            filename=self.tmp_heatmap_fpath)
        self.assertEqual(exists(self.tmp_heatmap_fpath), True)
        remove_files(set([self.tmp_heatmap_fpath]))


# run tests if called from command line
if __name__ == "__main__":
    main()
