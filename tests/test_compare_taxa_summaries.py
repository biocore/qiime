#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"
__status__ = "Development"

"""Test suite for the compare_taxa_summaries.py module."""

from numpy import array
from cogent.util.unit_test import TestCase, main

from qiime.compare_taxa_summaries import (add_filename_suffix,
        compare_taxa_summaries, _compute_all_to_expected_correlations,
        _compute_paired_sample_correlations, _format_correlation_vector,
        _format_taxa_summary, _get_correlation_function, _get_rank,
        _make_compatible_taxa_summaries, parse_sample_id_map,
        _sort_and_fill_taxa_summaries, _pearson_correlation,
        _spearman_correlation)

class CompareTaxaSummariesTests(TestCase):
    """Tests for the compare_taxa_summaries.py module."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.taxa_summary1 = (['Even1','Even2','Even3'],
            ['Bacteria;Actinobacteria;Actinobacteria(class);Actinobacteridae',
             'Bacteria;Bacteroidetes/Chlorobigroup;Bacteroidetes;Bacteroidia',
             'Bacteria;Firmicutes;Bacilli;Lactobacillales',
             'Bacteria;Firmicutes;Clostridia;Clostridiales',
             'Bacteria;Firmicutes;Erysipelotrichi;Erysipelotrichales',
             'Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales',
             'No blast hit;Other'],
            array([[0.0880247251673, 0.0721968465746, 0.081371761759],
                   [0.192137761955, 0.191095101593, 0.188504131885],
                   [0.0264895739603, 0.0259942669171, 0.0318460745596],
                   [0.491800007824, 0.526186212556, 0.49911159984],
                   [0.0311411916592, 0.0184083913576, 0.0282325481054],
                   [0.166137214246, 0.163087129528, 0.168923372865],
                   [0.00426952518811, 0.00303205147361, 0.0020105109874]]))

        self.taxa_summary2 = (['Even4','Even5','Even6'],
            ['Bacteria;Actinobacteria;Actinobacteria(class);NotARealTaxa',
             'Bacteria;AnotherFakeTaxa',
             'Bacteria;Bacteroidetes/Chlorobigroup;Bacteroidetes;Bacteroidia',
             'Bacteria;Firmicutes;Bacilli;Lactobacillales',
             'Bacteria;Firmicutes;Clostridia;Clostridiales',
             'Bacteria;Firmicutes;Erysipelotrichi;Erysipelotrichales',
             'Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales',
             'No blast hit;Other'],
            array([[0.99, 0.11, 0.075],
                   [0.1921, 0.19109, 0.18],
                   [0.192137761955, 0.191095101593, 0.188504131885],
                   [0.0264895739603, 0.0259942669171, 0.0318460745596],
                   [0.491800007824, 0.526186212556, 0.49911159984],
                   [0.0311411916592, 0.0184083913576, 0.0282325481054],
                   [0.166137214246, 0.163087129528, 0.168923372865],
                   [0.00426952518811, 0.00303205147361, 0.0020105109874]]))

        self.taxa_summary3 = (['Even7','Even8'], ['Eukarya'],
                                   array([[1.0, 1.0]]))

        # Has the same taxa as self.taxa_summary3_data (above).
        self.taxa_summary4 = (['Even1','Even2'], ['Eukarya'],
                                   array([[0.5, 0.6]]))

        # A sample ID map for testing making compatible taxa summaries.
        self.sample_id_map1 = {'Even7':'Even1'}

        # A sample ID map for testing making compatible taxa summaries using
        # two substitutions.
        self.sample_id_map2 = {'Even7':'Even1', 'Even8':'Even2'}

        # A sample ID map for testing many-to-one mappings.
        self.sample_id_map3 = {'Even7':'Even1', 'Even8':'Even1'}

        # A sample ID map with a bad sample ID as the value.
        self.sample_id_map4 = {'Even7':'Even1', 'Even8':'foo'}

        # A sample ID map with bad sample IDs as the keys.
        self.sample_id_map5 = {'foo':'Even1', 'bar':'Even2'}

        # No intersection with self.taxa_summary4_data.
        self.taxa_summary5 = (['Even1','Even2'], ['foo'],
                                   array([[0.5, 0.6]]))

        # Different sample ID from self.taxa_summary5.
        self.taxa_summary6 = (['Even1','Even3'], ['foo'],
                                   array([[0.1, 0.6]]))

        # Samples are in different orders from self.taxa_summary6.
        self.taxa_summary7 = (['Even3','Even1'], ['foo'],
                                   array([[0.2, 0.77]]))

        # Samples are not in alphabetical order, and we have multiple taxa.
        self.taxa_summary8 = (['Even3','Even1', 'S7'], ['foo', 'bar'],
                                   array([[0.2, 0.77, 0.001],
                                          [0.45, 0.9, 0.0]]))

        # For testing expected comparison mode.
        self.taxa_summary_exp1 = (['Expected'], ['Eukarya', 'Bacteria'],
                                   array([[0.5], [0.6]]))
        self.taxa_summary_obs1 = (['S1', 'S2'], ['Eukarya', 'Bacteria'],
                                   array([[0.4, 0.5], [0.5, 0.7]]))
        self.taxa_summary_obs1_mismatch = (['S1', 'S2'],
                                           ['Eukarya', 'Bacteria', 'Archaea'],
                                           array([[0.4, 0.5], [0.5, 0.7],
                                                  [0.1, 0.2]]))
        # For testing expected comparison mode using spearman (contains
        # repeats).
        self.taxa_summary_exp2 = (['Expected'],
                                  ['Eukarya', 'Bacteria', 'Archaea'],
                                   array([[0.5], [0.6], [0.4]]))
        self.taxa_summary_obs2 = (['S1', 'S2'],
                                  ['Eukarya', 'Bacteria', 'Archaea'],
                                  array([[0.4, 0.5], [0.5, 0.7], [0.4, 0.4]]))

        # For testing paired comparison mode.
        self.taxa_summary_paired1 = (['S1', 'S2'],
                ['Eukarya', 'Bacteria', 'Archaea'],
                array([[0.4, 0.5], [0.5, 0.7], [0.4, 0.4]]))
        self.taxa_summary_paired2 = (['S1', 'S2'],
                ['Eukarya', 'Bacteria', 'Archaea'],
                array([[0.5, 0.6], [0.7, 0.8], [0.5, 0.6]]))
        self.taxa_summary_paired3 = (['E1', 'E1'],
                ['Eukarya', 'Bacteria', 'Archaea'],
                array([[0.5, 0.6], [0.7, 0.8], [0.5, 0.6]]))
        self.taxa_summary_paired4 = (['E1', 'E2'],
                ['Eukarya', 'Bacteria', 'Archaea'],
                array([[0.5, 0.6], [0.7, 0.8], [0.5, 0.6]]))
        self.taxa_summary_paired5 = (['E1', 'E2'],
                ['Eukarya', 'Bacteria', 'Archaea', 'Foobar'],
                array([[0.5, 0.6], [0.7, 0.8], [0.5, 0.6], [0.1, 0.9]]))
        self.taxa_summary_paired_samp_id_map1 = {'S1':'E1', 'S2':'E1'}
        self.taxa_summary_paired_samp_id_map2 = {'S1':'E2'}

        # For formatting as correlation vectors.
        self.corr_vec1 = [('S1', 'T1', 0.7777777777)]
        self.corr_vec2 = [('S1', 'T1', 0.7777777777),
                          ('S2', 'T2', 0.1),
                          ('S3', 'T3', 100.68)]

        # Sample ID maps for testing the parser.
        self.sample_id_map_lines1 = ['\t\t\n', '', ' ', '\n', 'S1\tT1',
                                    'S2\tT2', '\n \t']
        self.sample_id_map_lines2 = ['\t\t\n', '', ' ', '\n', 'S1\tT1',
                                    'S2\tT2', '\n \t', 'S1\tT4']

        # For testing spearman correlation, taken from test_stats.BioEnvTests.
        self.a = [1,2,4,3,1,6,7,8,10,4]
        self.b = [2,10,20,1,3,7,5,11,6,13]
        self.c = [7,1,20,13,3,57,5,121,2,9]
        self.r = (1.7,10,20,1.7,3,7,5,11,6.5,13)
        self.x = (1, 2, 4, 3, 1, 6, 7, 8, 10, 4, 100, 2, 3, 77)

    def test_compare_taxa_summaries_expected_pearson(self):
        """Functions correctly using 'expected' comparison mode and pearson."""
        exp = ('Taxon\tS1\tS2\nBacteria\t0.5\t0.7\nEukarya\t0.4\t0.5\n',
               'Taxon\tExpected\nBacteria\t0.6\nEukarya\t0.5\n',
               '# Correlation coefficient: pearson\nS1\tExpected\t1.0000\nS2\t'
               'Expected\t1.0000\n')
        obs = compare_taxa_summaries(self.taxa_summary_obs1,
                self.taxa_summary_exp1, 'expected', 'pearson')
        self.assertEqual(obs, exp)

    def test_compare_taxa_summaries_expected_spearman(self):
        """Functions correctly using 'expected' comparison mode w/ spearman."""
        exp = ('Taxon\tS1\tS2\nBacteria\t0.5\t0.7\nEukarya\t0.4\t0.5\n',
               'Taxon\tExpected\nBacteria\t0.6\nEukarya\t0.5\n',
               '# Correlation coefficient: spearman\nS1\tExpected\t1.0000\nS2'
               '\tExpected\t1.0000\n')
        obs = compare_taxa_summaries(self.taxa_summary_obs1,
                self.taxa_summary_exp1, 'expected', 'spearman')
        self.assertEqual(obs, exp)

    def test_compare_taxa_summaries_paired_pearson(self):
        """Functions correctly using 'paired' comparison mode and pearson."""
        exp = ('Taxon\tS1\tS2\nArchaea\t0.4\t0.4\nBacteria\t0.5\t0.7\nEukarya'
               '\t0.4\t0.5\n',
               'Taxon\tS1\tS2\nArchaea\t0.5\t0.6\nBacteria\t0.7\t0.8\nEukarya'
               '\t0.5\t0.6\n',
               '# Correlation coefficient: pearson\nS1\tS1\t1.0000\nS2\tS2\t'
               '0.9449\n')
        obs = compare_taxa_summaries(self.taxa_summary_paired1,
                self.taxa_summary_paired2, 'paired', 'pearson')
        self.assertEqual(obs, exp)

    def test_compare_taxa_summaries_paired_sample_id_map(self):
        """Functions correctly using 'paired' comp mode and sample id map."""
        exp = ('Taxon\tS1\tS2\nArchaea\t0.4\t0.4\nBacteria\t0.5\t0.7\nEukarya'
               '\t0.4\t0.5\n', 'Taxon\tE1\tE2\nArchaea\t0.5\t0.6\nBacteria\t'
               '0.7\t0.8\nEukarya\t0.5\t0.6\n', '# Correlation coefficient: '
               'pearson\nS1\tE1\t1.0000\nS2\tE1\t0.9449\n')
        obs = compare_taxa_summaries(self.taxa_summary_paired1,
                self.taxa_summary_paired4, 'paired', 'pearson',
                self.taxa_summary_paired_samp_id_map1)
        self.assertEqual(obs, exp)

    def test_compare_taxa_summaries_paired_sample_id_map_partial(self):
        """Functions correctly using a partial sample id map."""
        exp = ('Taxon\tS1\tS2\nArchaea\t0.4\t0.4\nBacteria\t0.5\t0.7\nEukarya'
               '\t0.4\t0.5\n', 'Taxon\tE1\tE2\nArchaea\t0.5\t0.6\nBacteria\t'
               '0.7\t0.8\nEukarya\t0.5\t0.6\n', '# Correlation coefficient: '
               'pearson\nS1\tE2\t1.0000\n')
        obs = compare_taxa_summaries(self.taxa_summary_paired1,
                self.taxa_summary_paired4, 'paired', 'pearson',
                self.taxa_summary_paired_samp_id_map2)
        self.assertEqual(obs, exp)

    def test_compare_taxa_summaries_paired_sample_id_map_mismatched_taxa(self):
        """Functions correctly using a sample id map and mismatched taxa."""
        exp = ('Taxon\tS1\tS2\nArchaea\t0.4\t0.4\nBacteria\t0.5\t0.7\nEukarya'
               '\t0.4\t0.5\nFoobar\t0.0\t0.0\n', 'Taxon\tE1\tE2\nArchaea\t0.5'
               '\t0.6\nBacteria\t0.7\t0.8\nEukarya\t0.5\t0.6\nFoobar\t0.1\t0.9'
               '\n', '# Correlation coefficient: pearson\nS1\tE2\t-0.6264\n')
        obs = compare_taxa_summaries(self.taxa_summary_paired1,
                self.taxa_summary_paired5, 'paired', 'pearson',
                self.taxa_summary_paired_samp_id_map2)
        self.assertEqual(obs, exp)

    def test_compare_taxa_summaries_invalid_comparison_mode(self):
        """Throws error on unrecognized comparison mode."""
        self.assertRaises(ValueError, compare_taxa_summaries,
                self.taxa_summary_obs1, self.taxa_summary_exp1, 'foo',
                'pearson')

    def test_parse_sample_id_map(self):
        """Test parsing a sample id map functions correctly."""
        exp = {'S1':'T1', 'S2':'T2'}
        obs = parse_sample_id_map(self.sample_id_map_lines1)
        self.assertEqual(obs, exp)

    def test_parse_sample_id_map_repeat_sample_ids(self):
        """Test parsing a sample id map with non-unique first column fails."""
        self.assertRaises(ValueError, parse_sample_id_map,
                          self.sample_id_map_lines2)

    def test_add_filename_suffix(self):
        """Test adding a suffix to a filename works correctly."""
        self.assertEqual(add_filename_suffix('/foo/bar/baz.txt', 'z'),
                                             'bazz.txt')
        self.assertEqual(add_filename_suffix('baz.txt', 'z'),
                                             'bazz.txt')
        self.assertEqual(add_filename_suffix('/foo/bar/baz', 'z'),
                                             'bazz')
        self.assertEqual(add_filename_suffix('baz', 'z'),
                                             'bazz')
        self.assertEqual(add_filename_suffix('/baz.fasta.txt', 'z'),
                                             'baz.fastaz.txt')
        self.assertEqual(add_filename_suffix('baz.fasta.txt', 'z'),
                                             'baz.fastaz.txt')
        self.assertEqual(add_filename_suffix('/foo/', 'z'), 'z')

    def test_format_correlation_vector(self):
        """Test formatting correlations works correctly."""
        # One row.
        exp = 'S1\tT1\t0.7778\n'
        obs = _format_correlation_vector(self.corr_vec1)
        self.assertEqual(obs, exp)

        # Multiple rows.
        exp = 'S1\tT1\t0.7778\nS2\tT2\t0.1000\nS3\tT3\t100.6800\n'
        obs = _format_correlation_vector(self.corr_vec2)
        self.assertEqual(obs, exp)

    def test_format_correlation_vector_with_header(self):
        """Test formatting correlations with a header works correctly."""
        # One row.
        exp = 'foo\nS1\tT1\t0.7778\n'
        obs = _format_correlation_vector(self.corr_vec1, 'foo')
        self.assertEqual(obs, exp)

        # Multiple rows.
        exp = '#foobar\nS1\tT1\t0.7778\nS2\tT2\t0.1000\nS3\tT3\t100.6800\n'
        obs = _format_correlation_vector(self.corr_vec2, '#foobar')
        self.assertEqual(obs, exp)

    def test_format_taxa_summary(self):
        """Test formatting a taxa summary works correctly."""
        # More than one sample.
        exp = 'Taxon\tEven7\tEven8\nEukarya\t1.0\t1.0\n'
        obs = _format_taxa_summary(self.taxa_summary3)
        self.assertEqual(obs, exp)

        # More than one taxon.
        exp = 'Taxon\tExpected\nEukarya\t0.5\nBacteria\t0.6\nArchaea\t0.4\n'
        obs = _format_taxa_summary(self.taxa_summary_exp2)
        self.assertEqual(obs, exp)

    def test_get_correlation_function(self):
        """Test returns correct correlation function."""
        self.assertEqual(_get_correlation_function('pearson'),
                         _pearson_correlation)
        self.assertEqual(_get_correlation_function('spearman'),
                         _spearman_correlation)

    def test_get_correlation_function_invalid_correlation_type(self):
        """Test with an invalid (unrecognized) correlation type."""
        self.assertRaises(ValueError, _get_correlation_function, 'foo')

    def test_make_compatible_taxa_summaries(self):
        """Test making compatible taxa summaries works correctly on two ts."""
        exp = [(['Even1'], ['foo'], array([[ 0.5]])), (['Even1'], ['foo'],
                array([[ 0.1]]))]
        obs = _make_compatible_taxa_summaries(self.taxa_summary5,
                                              self.taxa_summary6)
        self.assertFloatEqual(obs, exp)

    def test_make_compatible_taxa_summaries_already_compatible(self):
        """Test on taxa summaries that are already compatible."""
        # Using two compatible ts.
        obs = _make_compatible_taxa_summaries(self.taxa_summary4,
                                              self.taxa_summary5)
        self.assertFloatEqual(obs, (self.taxa_summary4, self.taxa_summary5))

        # Using the same ts twice.
        obs = _make_compatible_taxa_summaries(self.taxa_summary4,
                                              self.taxa_summary4)
        self.assertFloatEqual(obs, (self.taxa_summary4, self.taxa_summary4))

        obs = _make_compatible_taxa_summaries(self.taxa_summary5,
                                              self.taxa_summary5)
        self.assertFloatEqual(obs, (self.taxa_summary5, self.taxa_summary5))

    def test_make_compatible_taxa_summaries_unordered(self):
        """Test on taxa summaries whose sample IDs are in different orders."""
        # Using two compatible ts.
        exp = [(['Even1', 'Even3'], ['foo'], array([[ 0.1,  0.6]])),
               (['Even1', 'Even3'], ['foo'], array([[ 0.77,  0.2 ]]))]
        obs = _make_compatible_taxa_summaries(self.taxa_summary6,
                                              self.taxa_summary7)
        self.assertFloatEqual(obs, exp)

    def test_make_compatible_taxa_summaries_sample_id_map_one_sub(self):
        """Test making compatible ts using a sample ID map with one sub."""
        # A simple 1-1 substitution.
        exp = ((['Even7'], ['Eukarya'], array([[ 1.]])),
               (['Even1'], ['Eukarya'], array([[ 0.5]])))
        obs = _make_compatible_taxa_summaries(self.taxa_summary3,
                self.taxa_summary4, self.sample_id_map1)
        self.assertFloatEqual(obs, exp)

    def test_make_compatible_taxa_summaries_sample_id_map_two_sub(self):
        """Test making compatible ts using a sample ID map with two subs."""
        # Two simple 1-1 substitutions.
        exp = ((['Even7', 'Even8'], ['Eukarya'], array([[ 1.,  1.]])),
               (['Even1', 'Even2'], ['Eukarya'], array([[ 0.5,  0.6]])))
        obs = _make_compatible_taxa_summaries(self.taxa_summary3,
                self.taxa_summary4, self.sample_id_map2)
        self.assertFloatEqual(obs, exp)

    def test_make_compatible_taxa_summaries_sample_id_map_many_to_one(self):
        """Test making compatible ts using sample ID map with M-1 sub."""
        # A many-to-1 substitution.
        exp = ((['Even7', 'Even8'], ['Eukarya'], array([[ 1.,  1.]])),
               (['Even1', 'Even1'], ['Eukarya'], array([[ 0.5,  0.5]])))
        obs = _make_compatible_taxa_summaries(self.taxa_summary3,
                self.taxa_summary4, self.sample_id_map3)
        self.assertFloatEqual(obs, exp)

    def test_make_compatible_taxa_summaries_incompatible(self):
        """Test on taxa summaries that have no common sample IDs."""
        self.assertRaises(ValueError, _make_compatible_taxa_summaries,
                          self.taxa_summary3, self.taxa_summary4)
        self.assertRaises(ValueError, _make_compatible_taxa_summaries,
            self.taxa_summary1, self.taxa_summary2)

    def test_make_compatible_taxa_summaries_bad_sample_id_map(self):
        """Test using a bad sample ID map."""
        self.assertRaises(ValueError, _make_compatible_taxa_summaries,
                          self.taxa_summary3, self.taxa_summary4,
                          self.sample_id_map4)
        self.assertRaises(ValueError, _make_compatible_taxa_summaries,
            self.taxa_summary3, self.taxa_summary4, self.sample_id_map5)

    def test_sort_and_fill_taxa_summaries(self):
        """Test _sort_and_fill_taxa_summaries functions as expected."""
        exp = [
            (['Even1','Even2','Even3'],
             ['Bacteria;Actinobacteria;Actinobacteria(class);Actinobacteridae',
              'Bacteria;Actinobacteria;Actinobacteria(class);NotARealTaxa',
              'Bacteria;AnotherFakeTaxa',
              'Bacteria;Bacteroidetes/Chlorobigroup;Bacteroidetes;Bacteroidia',
              'Bacteria;Firmicutes;Bacilli;Lactobacillales',
              'Bacteria;Firmicutes;Clostridia;Clostridiales',
              'Bacteria;Firmicutes;Erysipelotrichi;Erysipelotrichales',
              'Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales',
              'Eukarya',
              'No blast hit;Other'],
             array([[0.0880247251673, 0.0721968465746, 0.081371761759],
                    [0.,0.,0.],
                    [0.,0.,0.],
                    [0.192137761955, 0.191095101593, 0.188504131885],
                    [0.0264895739603, 0.0259942669171, 0.0318460745596],
                    [0.491800007824, 0.526186212556, 0.49911159984],
                    [0.0311411916592, 0.0184083913576, 0.0282325481054],
                    [0.166137214246, 0.163087129528, 0.168923372865],
                    [0.,0.,0.],
                    [0.00426952518811, 0.00303205147361, 0.0020105109874]])),
            (['Even4','Even5','Even6'],
             ['Bacteria;Actinobacteria;Actinobacteria(class);Actinobacteridae',
              'Bacteria;Actinobacteria;Actinobacteria(class);NotARealTaxa',
              'Bacteria;AnotherFakeTaxa',
              'Bacteria;Bacteroidetes/Chlorobigroup;Bacteroidetes;Bacteroidia',
              'Bacteria;Firmicutes;Bacilli;Lactobacillales',
              'Bacteria;Firmicutes;Clostridia;Clostridiales',
              'Bacteria;Firmicutes;Erysipelotrichi;Erysipelotrichales',
              'Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales',
              'Eukarya',
              'No blast hit;Other'],
             array([[0.,0.,0.],
                    [0.99, 0.11, 0.075],
                    [0.1921, 0.19109, 0.18],
                    [0.192137761955, 0.191095101593, 0.188504131885],
                    [0.0264895739603, 0.0259942669171, 0.0318460745596],
                    [0.491800007824, 0.526186212556, 0.49911159984],
                    [0.0311411916592, 0.0184083913576, 0.0282325481054],
                    [0.166137214246, 0.163087129528, 0.168923372865],
                    [0.,0.,0.],
                    [0.00426952518811, 0.00303205147361, 0.0020105109874]])),
            (['Even7','Even8'],
             ['Bacteria;Actinobacteria;Actinobacteria(class);Actinobacteridae',
              'Bacteria;Actinobacteria;Actinobacteria(class);NotARealTaxa',
              'Bacteria;AnotherFakeTaxa',
              'Bacteria;Bacteroidetes/Chlorobigroup;Bacteroidetes;Bacteroidia',
              'Bacteria;Firmicutes;Bacilli;Lactobacillales',
              'Bacteria;Firmicutes;Clostridia;Clostridiales',
              'Bacteria;Firmicutes;Erysipelotrichi;Erysipelotrichales',
              'Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales',
              'Eukarya',
              'No blast hit;Other'],
             array([[0.,0.],
                    [0.,0.],
                    [0.,0.],
                    [0.,0.],
                    [0.,0.],
                    [0.,0.],
                    [0.,0.],
                    [0.,0.],
                    [1.,1.],
                    [0.,0.]]))
        ]

        obs = _sort_and_fill_taxa_summaries([self.taxa_summary1,
                                             self.taxa_summary2,
                                             self.taxa_summary3])
        self.assertFloatEqual(obs, exp)

    def test_sort_and_fill_taxa_summaries_same(self):
        """Test when the taxa summaries are the same."""
        exp = [(['Even7','Even8'], ['Eukarya'], array([[1.0, 1.0]])),
               (['Even7','Even8'], ['Eukarya'], array([[1.0, 1.0]]))]
        obs = _sort_and_fill_taxa_summaries([self.taxa_summary3,
                                             self.taxa_summary3])
        self.assertFloatEqual(obs, exp)

        exp = [(['Even7','Even8'], ['Eukarya'], array([[1.0, 1.0]])),
               (['Even1','Even2'], ['Eukarya'], array([[0.5, 0.6]]))]
        obs = _sort_and_fill_taxa_summaries([self.taxa_summary3,
                                             self.taxa_summary4])
        self.assertFloatEqual(obs, exp)

        # Test the other direction.
        exp = [(['Even1','Even2'], ['Eukarya'], array([[0.5, 0.6]])),
               (['Even7','Even8'], ['Eukarya'], array([[1.0, 1.0]]))]
        obs = _sort_and_fill_taxa_summaries([self.taxa_summary4,
                                             self.taxa_summary3])
        self.assertFloatEqual(obs, exp)

    def test_sort_and_fill_taxa_summaries_no_intersection(self):
        """Test when the taxa summaries have no intersection."""
        exp = [(['Even1','Even2'], ['Eukarya', 'foo'], array([[0.5, 0.6],
                                                              [0.0, 0.0]])),
               (['Even1','Even2'], ['Eukarya', 'foo'], array([[0.0, 0.0],
                                                              [0.5, 0.6]]))]
        obs = _sort_and_fill_taxa_summaries([self.taxa_summary4,
                                             self.taxa_summary5])
        self.assertFloatEqual(obs, exp)

    def test_compute_all_to_expected_correlations_pearson(self):
        """Test functions correctly with pearson correlation."""
        exp = [('S1', 'Expected', 1.0), ('S2', 'Expected', 1.0)]
        obs = _compute_all_to_expected_correlations(self.taxa_summary_obs1,
                                       self.taxa_summary_exp1,
                                       _pearson_correlation)
        self.assertFloatEqual(obs, exp)

    def test_compute_all_to_expected_correlations_spearman(self):
        """Test functions correctly with spearman correlation."""
        # No repeats in ranks.
        exp = [('S1', 'Expected', 1.0), ('S2', 'Expected', 1.0)]
        obs = _compute_all_to_expected_correlations(self.taxa_summary_obs1,
                                       self.taxa_summary_exp1,
                                       _spearman_correlation)
        self.assertFloatEqual(obs, exp)

        # Repeats in ranks.
        exp = [('S1', 'Expected', 0.866025), ('S2', 'Expected', 1.0)]
        obs = _compute_all_to_expected_correlations(self.taxa_summary_obs2,
                                       self.taxa_summary_exp2,
                                       _spearman_correlation)
        self.assertFloatEqual(obs, exp)

    def test_compute_all_to_expected_correlations_invalid_sample_count(self):
        """Test running on expected taxa summary without exactly one sample."""
        self.assertRaises(ValueError, _compute_all_to_expected_correlations,
                          self.taxa_summary1, self.taxa_summary2,
                          _pearson_correlation)

    def test_compute_all_to_expected_correlations_incompatible_taxa(self):
        """Test running on taxa summaries that have mismatched taxa info."""
        self.assertRaises(ValueError, _compute_all_to_expected_correlations,
                          self.taxa_summary_obs1_mismatch,
                          self.taxa_summary_exp1,
                          _pearson_correlation)

    def test_compute_paired_sample_correlations(self):
        """Test runs correctly on standard input taxa summaries."""
        exp = [('S1', 'S1', 1.0), ('S2', 'S2', 0.94491118252306538)]
        obs = _compute_paired_sample_correlations(self.taxa_summary_paired1,
                self.taxa_summary_paired2, _pearson_correlation)
        self.assertFloatEqual(obs, exp)

    def test_compute_paired_sample_correlations_mismatched_ids(self):
        """Test runs correctly on taxa summaries with mismatched IDs."""
        exp = [('S1', 'E1', 1.0), ('S2', 'E1', 0.94491118252306538)]
        obs = _compute_paired_sample_correlations(self.taxa_summary_paired1,
                self.taxa_summary_paired3, _pearson_correlation)
        self.assertFloatEqual(obs, exp)

    def test_compute_paired_sample_correlations_incompatible_taxa(self):
        """Test running on taxa summaries that have mismatched taxa info."""
        self.assertRaises(ValueError, _compute_paired_sample_correlations,
                          self.taxa_summary4, self.taxa_summary5,
                          _pearson_correlation)

    def test_compute_paired_sample_correlations_incompatible_samples(self):
        """Test on incompatible taxa summaries (mismatched sample lengths)."""
        self.assertRaises(ValueError, _compute_paired_sample_correlations,
                          self.taxa_summary1, self.taxa_summary3,
                          _pearson_correlation)

    def test_pearson_correlation_invalid_input(self):
        """Test running pearson correlation on bad input."""
        self.assertRaises(ValueError, _pearson_correlation,
                          [1.4, 2.5], [5.6, 8.8, 9.0])
        self.assertRaises(ValueError, _pearson_correlation, [1.4], [5.6])

    # The following tests are taken from test_stats.BioEnvTests to test the
    # spearman function. They have been modified only to work with the function
    # instead of the BioEnv object containing that method.
    def test_get_rank(self):
        """Test the _get_rank method with valid input."""
        exp = ([1.5,3.5,7.5,5.5,1.5,9.0,10.0,11.0,12.0,7.5,14.0,3.5,5.5,13.0],
               4)
        obs = _get_rank(self.x)
        self.assertFloatEqual(exp,obs)

        exp = ([1.5,3.0,5.5,4.0,1.5,7.0,8.0,9.0,10.0,5.5],2)
        obs = _get_rank(self.a)
        self.assertFloatEqual(exp,obs)

        exp = ([2,7,10,1,3,6,4,8,5,9],0)
        obs = _get_rank(self.b)
        self.assertFloatEqual(exp,obs)

        exp = ([1.5,7.0,10.0,1.5,3.0,6.0,4.0,8.0,5.0,9.0], 1)
        obs = _get_rank(self.r)
        self.assertFloatEqual(exp,obs)

        exp = ([],0)
        obs = _get_rank([])
        self.assertEqual(exp,obs)

    def test_get_rank_invalid_input(self):
        """Test the _get_rank method with invalid input."""
        vec = [1, 'a', 3, 2.5, 3, 1]
        self.assertRaises(TypeError, _get_rank, vec)

        vec = [1, 2, {1:2}, 2.5, 3, 1]
        self.assertRaises(TypeError, _get_rank, vec)

        vec = [1, 2, [23,1], 2.5, 3, 1]
        self.assertRaises(TypeError, _get_rank, vec)

        vec = [1, 2, (1,), 2.5, 3, 1]
        self.assertRaises(TypeError, _get_rank, vec)

    def test_spearman_correlation(self):
        """Test the _spearman_correlation method."""
        # One vector has no ties.
        exp = 0.3719581
        obs = _spearman_correlation(self.a,self.b)
        self.assertFloatEqual(exp,obs)

        # Both vectors have no ties.
        exp = 0.2969697
        obs = _spearman_correlation(self.b,self.c)
        self.assertFloatEqual(exp,obs)

        # Both vectors have ties.
        exp = 0.388381
        obs = _spearman_correlation(self.a,self.r)
        self.assertFloatEqual(exp,obs)

    def test_spearman_correlation_invalid_input(self):
        """Test the _spearman_correlation method with invalid input."""
        self.assertRaises(ValueError, _spearman_correlation, [],[])
        self.assertRaises(ValueError, _spearman_correlation, self.a,[])
        self.assertRaises(ValueError, _spearman_correlation, {0:2}, [1,2,3])


if __name__ == "__main__":
    main()
