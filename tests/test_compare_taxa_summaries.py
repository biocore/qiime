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
        compare_taxa_summaries, _compute_correlation,
        _format_correlation_vector, _format_taxa_summary,
        _get_correlation_function, _get_rank, _make_compatible_taxa_summaries,
        parse_sample_id_map, _sort_and_fill_taxa_summaries,
        _pearson_correlation, _spearman, _spearman_correlation)

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
        self.sample_id_map1 = {'Even7':'1', 'Even1':'1', 'Even8':'2',
                               'Even2':'2'}

        # A sample ID map with an extra sample ID not found in the taxa
        # summaries.
        self.sample_id_map2 = {'Even7':'1', 'Even1':'1', 'Even8':'2',
                               'Even2':'2', 'Foo':'3'}

        # A sample ID map without all original sample IDs mapped.
        self.sample_id_map3 = {'Even7':'1', 'Even1':'1'}

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

        # For testing expected comparison mode with an expected sample ID (same
        # data as self.taxa_summary_exp2, but with extra sample).
        self.taxa_summary_exp3 = (['Expected', 'Expected2'],
                                  ['Eukarya', 'Bacteria', 'Archaea'],
                                   array([[0.5, 0.1], [0.6, 0.1], [0.4, 0.1]]))

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
        # Full mapping.
        self.taxa_summary_paired_samp_id_map1 = {'S1':'a', 'S2':'b',
                                                 'E1':'a','E2':'b'}
        # Partial mapping.
        self.taxa_summary_paired_samp_id_map2 = {'S1':'a', 'S2':'b',
                                                 'E1':'c', 'E2':'a'}

        # For formatting as correlation vectors.
        self.corr_vec1 = [('S1', 'T1', 0.7777777777, 0, 0)]
        self.corr_vec2 = [('S1', 'T1', 0.7777777777, 0, 0),
                          ('S2', 'T2', 0.1, 0.05, 0.15),
                          ('S3', 'T3', 100.68, 0.9, 1)]

        # Sample ID maps for testing the parser.
        self.sample_id_map_lines1 = ['\t\t\n', '', ' ', '\n', 'S1\ta',
                                    'S2\tb', '\n \t', 'T1\ta', 'T2\tb']
        self.sample_id_map_lines2 = ['\t\t\n', '', ' ', '\n', 'S1\ta',
                                    'S2\tb', '\n \t', 'S1\tc']
        self.sample_id_map_lines3 = ['S1\ta', 'T1\ta', 'S2\ta']

        # For testing _spearman_correlation, taken from the Spearman wikipedia
        # article.
        self.spearman_data1 = [106, 86, 100, 101, 99, 103, 97, 113, 112, 110]
        self.spearman_data2 = [7, 0, 27, 50, 28, 29, 20, 12, 6, 17]
        
        # For testing _spearman, taken from test_stats.BioEnvTests.
        self.a = [1,2,4,3,1,6,7,8,10,4]
        self.b = [2,10,20,1,3,7,5,11,6,13]
        self.c = [7,1,20,13,3,57,5,121,2,9]
        self.r = (1.7,10,20,1.7,3,7,5,11,6.5,13)
        self.x = (1, 2, 4, 3, 1, 6, 7, 8, 10, 4, 100, 2, 3, 77)

    def test_compare_taxa_summaries_expected_pearson(self):
        """Functions correctly using 'expected' comparison mode and pearson."""
        # Overall correlation.
        exp = ('Taxon\tS1\tS2\nBacteria\t0.5\t0.7\nEukarya\t0.4\t0.5\n',
               'Taxon\tExpected\nBacteria\t0.6\nEukarya\t0.5\n',
               '# Correlation coefficient: pearson. Performed a two-tailed '
               'test of significance using a t-distribution.\nCorrelation '
               'coefficient\tp-value\n0.6882\t0.3118\n', None)
        obs = compare_taxa_summaries(self.taxa_summary_obs1,
                self.taxa_summary_exp1, 'expected', 'pearson')
        self.assertEqual(obs, exp)

        # Broken down by sample pairs.
        exp = ('Taxon\tS1\tS2\nBacteria\t0.5\t0.7\nEukarya\t0.4\t0.5\n',
               'Taxon\tExpected\nBacteria\t0.6\nEukarya\t0.5\n',
               '# Correlation coefficient: pearson. Performed a two-tailed '
               'test of significance using a t-distribution.\nCorrelation '
               'coefficient\tp-value\n0.6882\t0.3118\n',
               '# Correlation coefficient: pearson. Performed a two-tailed '
               'test of significance using a t-distribution.\nSample ID\t'
               'Sample ID\tCorrelation coefficient\tp-value\tp-value '
               '(Bonferroni-corrected)\nS1\tExpected\t1.0000\t1.0000\t1.0000\n'
               'S2\tExpected\t1.0000\t1.0000\t1.0000\n')
        obs = compare_taxa_summaries(self.taxa_summary_obs1,
                self.taxa_summary_exp1, 'expected', 'pearson', True)
        self.assertEqual(obs, exp)

    def test_compare_taxa_summaries_spearman_warning(self):
        """Functions correctly using spearman with <= 10 obs."""
        # Verified with R's cor.test function.
        exp = ('Taxon\tS1\tS2\nBacteria\t0.5\t0.7\nEukarya\t0.4\t0.5\n',
            'Taxon\tExpected\nBacteria\t0.6\nEukarya\t0.5\n',
            "# Correlation coefficient: spearman. Performed a two-tailed "
            "test of significance using a t-distribution.\n# Since there "
            "were 10 or fewer observations when calculating Spearman's "
            "rank correlation coefficient, the p-value is not accurate "
            "when using the t-distribution. Please see Biometry (Sokal "
            "and Rohlf, 3rd edition) page 600 for more details.\nCorrelation "
            "coefficient\tp-value\n0.7071\t0.2929\n",
            '# Correlation coefficient: spearman. Performed a two-tailed '
            'test of significance using a t-distribution.\n# Since there '
            'were 10 or fewer taxa in the sorted and filled taxa '
            'summary files, the p-values and Bonferroni-corrected p-values '
            'are not accurate when using the t-distribution. Please see '
            'Biometry (Sokal and Rohlf, 3rd edition) page 600 for more '
            'details.\nSample ID\tSample ID\tCorrelation coefficient\t'
            'p-value\tp-value (Bonferroni-corrected)\nS1\tExpected\t1.0000\t'
            '1.0000\t1.0000\nS2\tExpected\t1.0000\t1.0000\t1.0000\n')
        obs = compare_taxa_summaries(self.taxa_summary_obs1,
                self.taxa_summary_exp1, 'expected', 'spearman', True)
        self.assertEqual(obs, exp)

    def test_compare_taxa_summaries_spearman_no_warning(self):
        """Functions correctly using spearman with > 10 obs."""
        # Verified with R's cor.test function.

        # We should receive a warning only for the individual tests, but not
        # for the overall test.
        exp = ('# Correlation coefficient: spearman. Performed a '
            'two-tailed test of significance using a t-distribution.\n'
            '# Number of samples that matched between the taxa summary '
            'files: 3\nCorrelation coefficient\tp-value\n0.2629\t0.1853\n',
            '# Correlation coefficient: spearman. Performed a two-tailed '
            'test of significance using a t-distribution.\n# Number of '
            'samples that matched between the taxa summary files: 3\n'
            '# Since there were 10 or fewer taxa in the sorted and filled '
            'taxa summary files, the p-values and Bonferroni-corrected '
            'p-values are not accurate when using the t-distribution. '
            'Please see Biometry (Sokal and Rohlf, 3rd edition) page 600 '
            'for more details.\nSample ID\tSample ID\t'
            'Correlation coefficient\tp-value\tp-value '
            '(Bonferroni-corrected)\nEven1\tEven4\t0.0753\t0.8473\t'
            '1.0000\nEven2\tEven5\t0.4017\t0.2839\t0.8517\nEven3\t'
            'Even6\t0.4017\t0.2839\t0.8517\n')
        sample_id_map = {'Even1':'a', 'Even2':'b', 'Even3':'c',
                         'Even4':'a', 'Even5':'b', 'Even6':'c'}
        obs = compare_taxa_summaries(self.taxa_summary1,
                self.taxa_summary2, 'paired', 'spearman', True, sample_id_map)
        self.assertEqual(obs[2:], exp)

    def test_compare_taxa_summaries_paired_pearson(self):
        """Functions correctly using 'paired' comparison mode and pearson."""
        exp = ('Taxon\tS1\tS2\nArchaea\t0.4\t0.4\nBacteria\t0.5\t'
            '0.7\nEukarya\t0.4\t0.5\n', 'Taxon\tS1\tS2\nArchaea\t0.5\t'
            '0.6\nBacteria\t0.7\t0.8\nEukarya\t0.5\t0.6\n',
            '# Correlation coefficient: pearson. Performed a two-tailed '
            'test of significance using a t-distribution.\n'
            '# Number of samples that matched between the taxa summary '
            'files: 2\nCorrelation coefficient\tp-value\n0.9024\t'
            '0.0138\n', None)
        obs = compare_taxa_summaries(self.taxa_summary_paired1,
                self.taxa_summary_paired2, 'paired', 'pearson')
        self.assertEqual(obs, exp)

    def test_compare_taxa_summaries_paired_sample_id_map(self):
        """Functions correctly using 'paired' comp mode and sample id map."""
        exp = ('Taxon\tS1\tS2\nArchaea\t0.4\t0.4\nBacteria\t0.5\t'
            '0.7\nEukarya\t0.4\t0.5\n', 'Taxon\tE1\tE2\nArchaea\t0.5\t'
            '0.6\nBacteria\t0.7\t0.8\nEukarya\t0.5\t0.6\n',
            '# Correlation coefficient: pearson. Performed a two-tailed '
            'test of significance using a t-distribution.\n'
            '# Number of samples that matched between the taxa summary '
            'files: 2\nCorrelation coefficient\tp-value\n0.9024\t'
            '0.0138\n', None)
        obs = compare_taxa_summaries(self.taxa_summary_paired1,
                self.taxa_summary_paired4, 'paired', 'pearson', False,
                self.taxa_summary_paired_samp_id_map1)
        self.assertEqual(obs, exp)

    def test_compare_taxa_summaries_paired_sample_id_map_partial(self):
        """Functions correctly using a partial sample id map."""
        # The sample ID map has some mappings that are not complete- i.e. a
        # sample from one file has a new sample ID that doesn't match any other
        # new sample IDs. In this case, the sample should be ignored.
        exp = ('Taxon\tS1\tS2\nArchaea\t0.4\t0.4\nBacteria\t0.5\t'
            '0.7\nEukarya\t0.4\t0.5\n', 'Taxon\tE1\tE2\nArchaea\t0.5'
            '\t0.6\nBacteria\t0.7\t0.8\nEukarya\t0.5\t0.6\n',
            '# Correlation coefficient: pearson. Performed a two-tailed '
            'test of significance using a t-distribution.\n# Number of '
            'samples that matched between the taxa summary '
            'files: 1\nCorrelation coefficient\tp-value\n1.0000\t0.0000\n',
            '# Correlation coefficient: pearson. Performed a two-tailed '
            'test of significance using a t-distribution.\n# Number of '
            'samples that matched between the taxa summary files: 1\n'
            'Sample ID\tSample ID\tCorrelation coefficient\tp-value\t'
            'p-value (Bonferroni-corrected)\nS1\tE2\t1.0000\t0.0000\t0.0000\n')
        obs = compare_taxa_summaries(self.taxa_summary_paired1,
                self.taxa_summary_paired4, 'paired', 'pearson', True,
                self.taxa_summary_paired_samp_id_map2)
        self.assertEqual(obs, exp)

    def test_compare_taxa_summaries_paired_sample_id_map_mismatched_taxa(self):
        """Functions correctly using a sample id map and mismatched taxa."""
        exp = ('Taxon\tS1\tS2\nArchaea\t0.4\t0.4\nBacteria\t0.5\t0.7\nEukarya'
               '\t0.4\t0.5\nFoobar\t0.0\t0.0\n', 'Taxon\tE1\tE2\nArchaea\t0.5'
               '\t0.6\nBacteria\t0.7\t0.8\nEukarya\t0.5\t0.6\nFoobar\t0.1\t0.9'
               '\n', '# Correlation coefficient: pearson. Performed a '
               'two-tailed test of significance using a t-distribution.\n'
               '# Number of samples that matched between the taxa '
               'summary files: 1\nCorrelation coefficient\tp-value\n'
               '-0.6264\t0.3736\n', None)
        obs = compare_taxa_summaries(self.taxa_summary_paired1,
                self.taxa_summary_paired5, 'paired', 'pearson', False,
                self.taxa_summary_paired_samp_id_map2)
        self.assertEqual(obs, exp)

    def test_compare_taxa_summaries_invalid_comparison_mode(self):
        """Throws error on unrecognized comparison mode."""
        self.assertRaises(ValueError, compare_taxa_summaries,
                self.taxa_summary_obs1, self.taxa_summary_exp1, 'foo',
                'pearson')

    def test_parse_sample_id_map(self):
        """Test parsing a sample id map functions correctly."""
        exp = {'S1':'a', 'S2':'b', 'T1':'a', 'T2':'b'}
        obs = parse_sample_id_map(self.sample_id_map_lines1)
        self.assertEqual(obs, exp)

    def test_parse_sample_id_map_repeat_sample_ids(self):
        """Test parsing a sample id map with non-unique first column fails."""
        self.assertRaises(ValueError, parse_sample_id_map,
                          self.sample_id_map_lines2)

    def test_parse_sample_id_map_many_to_one_mapping(self):
        """Test parsing a sample id map with many-to-one mapping fails."""
        self.assertRaises(ValueError, parse_sample_id_map,
                          self.sample_id_map_lines3)

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
        exp = 'Sample ID\tSample ID\tCorrelation coefficient\tp-value\t' + \
              'p-value (Bonferroni-corrected)\nS1\tT1\t0.7778\t0.0000\t' + \
              '0.0000\n'
        obs = _format_correlation_vector(self.corr_vec1)
        self.assertEqual(obs, exp)

        # Multiple rows.
        exp = 'Sample ID\tSample ID\tCorrelation coefficient\tp-value\t' + \
              'p-value (Bonferroni-corrected)\nS1\tT1\t0.7778\t0.0000\t' + \
              '0.0000\nS2\tT2\t0.1000\t0.0500\t0.1500\nS3\tT3\t100.6800\t' + \
              '0.9000\t1.0000\n'
        obs = _format_correlation_vector(self.corr_vec2)
        self.assertEqual(obs, exp)

    def test_format_correlation_vector_with_header(self):
        """Test formatting correlations with a header works correctly."""
        # One row.
        exp = 'foo\nSample ID\tSample ID\tCorrelation coefficient\t' + \
              'p-value\tp-value (Bonferroni-corrected)\nS1\tT1\t0.7778\t' + \
              '0.0000\t0.0000\n'
        obs = _format_correlation_vector(self.corr_vec1, 'foo')
        self.assertEqual(obs, exp)

        # Multiple rows.
        exp = '#foobar\nSample ID\tSample ID\tCorrelation coefficient\t' + \
              'p-value\tp-value (Bonferroni-corrected)\nS1\tT1\t0.7778\t' + \
              '0.0000\t0.0000\nS2\tT2\t0.1000\t0.0500\t0.1500\nS3\tT3\t' + \
              '100.6800\t0.9000\t1.0000\n'
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

    def test_make_compatible_taxa_summaries_sample_id_map(self):
        """Test making compatible ts using a sample ID map."""
        exp = ((['Even7', 'Even8'], ['Eukarya'], array([[ 1., 1.]])),
               (['Even1', 'Even2'], ['Eukarya'], array([[ 0.5, 0.6]])))
        obs = _make_compatible_taxa_summaries(self.taxa_summary3,
                self.taxa_summary4, self.sample_id_map1)
        self.assertFloatEqual(obs, exp)

    def test_make_compatible_taxa_summaries_sample_id_map_extra_sample(self):
        """Test using a sample ID map with an extra sample ID."""
        exp = ((['Even7', 'Even8'], ['Eukarya'], array([[ 1.,  1.]])),
               (['Even1', 'Even2'], ['Eukarya'], array([[ 0.5,  0.6]])))
        obs = _make_compatible_taxa_summaries(self.taxa_summary3,
                self.taxa_summary4, self.sample_id_map2)
        self.assertFloatEqual(obs, exp)

    def test_make_compatible_taxa_summaries_sample_id_map_incomplete_map(self):
        """Test using a sample ID map that is incomplete."""
        self.assertRaises(ValueError, _make_compatible_taxa_summaries,
                self.taxa_summary3, self.taxa_summary4, self.sample_id_map3)

    def test_make_compatible_taxa_summaries_incompatible(self):
        """Test on taxa summaries that have no common sample IDs."""
        self.assertRaises(ValueError, _make_compatible_taxa_summaries,
                          self.taxa_summary3, self.taxa_summary4)
        self.assertRaises(ValueError, _make_compatible_taxa_summaries,
            self.taxa_summary1, self.taxa_summary2)

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

    def test_compute_correlation_expected_pearson(self):
        """Test functions correctly with expected mode and pearson corr."""
        exp = ((0.68824720161169595, 0.31175279838830405), None)
        obs = _compute_correlation(self.taxa_summary_obs1,
                self.taxa_summary_exp1, 'expected', _pearson_correlation)
        self.assertFloatEqual(obs, exp)

    def test_compute_correlation_expected_pearson_detailed(self):
        """Test functions with expected mode, pearson corr, detailed tests."""
        exp = ((0.68824720161169595, 0.31175279838830405),
               [('S1', 'Expected', 1.0, 1, 1), ('S2', 'Expected', 1.0, 1, 1)])
        obs = _compute_correlation(self.taxa_summary_obs1,
                self.taxa_summary_exp1, 'expected', _pearson_correlation, True)
        self.assertFloatEqual(obs, exp)

    def test_compute_correlation_expected_spearman_detailed(self):
        """Test functions with expected mode, spearman corr, detailed tests."""
        # Verified with R's cor.test function.

        # No repeats in ranks.
        exp = ((0.707106781187, 0.292893218813),
               [('S1', 'Expected', 1.0, 1, 1), ('S2', 'Expected', 1.0, 1, 1)])
        obs = _compute_correlation(self.taxa_summary_obs1,
                                   self.taxa_summary_exp1,
                                   'expected',
                                   _spearman_correlation, True)
        self.assertFloatEqual(obs, exp)

        # Repeats in ranks.
        exp = ((0.839146391678, 0.0367298712193),
               [('S1', 'Expected', 0.866025, 0.33333333333333326,
                 0.66666666666666652),
                ('S2', 'Expected', 1.0, 0, 0)])
        obs = _compute_correlation(self.taxa_summary_obs2,
                                   self.taxa_summary_exp2,
                                   'expected',
                                   _spearman_correlation, True)
        self.assertFloatEqual(obs, exp)

    def test_compute_correlation_expected_expected_sample_id(self):
        """Test functions with expected mode using an expected sample ID."""
        exp = ((0.839146391678, 0.0367298712193),
               [('S1', 'Expected', 0.866025, 0.33333333333333326,
                 0.66666666666666652),
                ('S2', 'Expected', 1.0, 0, 0)])

        # Using a single-sample expected ts.
        obs = _compute_correlation(self.taxa_summary_obs2,
                                   self.taxa_summary_exp2,
                                   'expected',
                                   _spearman_correlation, True, 'Expected')
        self.assertFloatEqual(obs, exp)

        # Using a two-sample expected ts.
        obs = _compute_correlation(self.taxa_summary_obs2,
                                   self.taxa_summary_exp3,
                                   'expected',
                                   _spearman_correlation, True, 'Expected')
        self.assertFloatEqual(obs, exp)

    def test_compute_correlation_expected_invalid_sample_count(self):
        """Test running on expected taxa summary without exactly one sample."""
        self.assertRaises(ValueError, _compute_correlation,
                          self.taxa_summary1, self.taxa_summary2, 'expected',
                          _pearson_correlation)

    def test_compute_correlation_expected_invalid_expected_sample_id(self):
        """Test running in expected mode with invalid expected sample ID."""
        self.assertRaises(ValueError, _compute_correlation,
                self.taxa_summary_obs2, self.taxa_summary_exp3, 'expected',
                _pearson_correlation, expected_sample_id='foo')

    def test_compute_correlation_invalid_comparison_mode(self):
        """Test running using an invalid comparison mode."""
        self.assertRaises(ValueError, _compute_correlation,
                          self.taxa_summary1, self.taxa_summary2, 'foo',
                          _pearson_correlation)

    def test_compute_correlation_incompatible_taxa(self):
        """Test running on taxa summaries that have mismatched taxa info."""
        self.assertRaises(ValueError, _compute_correlation,
                          self.taxa_summary_obs1_mismatch,
                          self.taxa_summary_exp1,
                          'expected', _pearson_correlation)

    def test_compute_correlation_paired(self):
        """Test runs correctly on standard input taxa summaries."""
        # Verified using R's cor.test function.
        exp = ((0.90243902439024193, 0.013812916237431808), None)
        obs = _compute_correlation(self.taxa_summary_paired1,
                self.taxa_summary_paired2, 'paired', _pearson_correlation)
        self.assertFloatEqual(obs, exp)

    def test_compute_correlation_paired_detailed(self):
        """Test runs correctly on standard input ts with detailed tests."""
        # Verified using R's cor.test function.
        exp = ((0.90243902439024193, 0.013812916237431808),
               [('S1', 'S1', 1.0, 0, 0),
                ('S2', 'S2', 0.94491118252306538, 0.21229561501,
                 0.424591230019)])
        obs = _compute_correlation(self.taxa_summary_paired1,
                self.taxa_summary_paired2, 'paired', _pearson_correlation,
                True)
        self.assertFloatEqual(obs, exp)

    def test_compute_correlation_paired_mismatched_ids(self):
        """Test runs correctly on taxa summaries with mismatched IDs."""
        # Verified using R's cor.test function.
        exp = ((0.90243902439024193, 0.013812916237431808),
               [('S1', 'E1', 1.0, 0, 0),
                ('S2', 'E1', 0.94491118252306538, 0.21229561501,
                 0.424591230019)])
        obs = _compute_correlation(self.taxa_summary_paired1,
                self.taxa_summary_paired3, 'paired', _pearson_correlation,
                True)
        self.assertFloatEqual(obs, exp)

    def test_compute_correlation_paired_incompatible_samples(self):
        """Test on incompatible taxa summaries (mismatched sample lengths)."""
        self.assertRaises(ValueError, _compute_correlation,
                          self.taxa_summary1, self.taxa_summary3,
                          'paired', _spearman_correlation)

    def test_pearson_correlation_invalid_input(self):
        """Test running pearson correlation on bad input."""
        self.assertRaises(ValueError, _pearson_correlation,
                          [1.4, 2.5], [5.6, 8.8, 9.0])
        self.assertRaises(ValueError, _pearson_correlation, [1.4], [5.6])

    def test_spearman_correlation(self):
        """Test running spearman correlation on valid input."""
        # This example taken from Wikipedia page:
        # http://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient
        #
        # The p-value is off because the example uses a one-tailed test, while
        # we use a two-tailed test. Someone confirms the answer that we get
        # here for a two-tailed test:
        # http://stats.stackexchange.com/questions/22816/calculating-p-value-
        #     for-spearmans-rank-correlation-coefficient-example-on-wikip
        exp = (-0.17575757575757578, 0.62718834477648433)
        obs = _spearman_correlation(self.spearman_data1, self.spearman_data2)
        self.assertFloatEqual(obs, exp)

    def test_spearman_correlation_one_obs(self):
        """Test running spearman correlation on a single observation."""
        self.assertRaises(ValueError, _spearman_correlation, [1.0], [5.0])

    def test_spearman_no_variation(self):
        """Test the _spearman function with a vector having no variation."""
        exp = 0.0
        obs = _spearman([1, 1, 1], [1, 2, 3])
        self.assertEqual(obs, exp)

    # The following tests are taken from test_stats.BioEnvTests to test the
    # spearman function. They have been modified only to work with the function
    # instead of the BioEnv object containing that method.
    def test_spearman(self):
        """Test the _spearman function."""
        # One vector has no ties.
        exp = 0.3719581
        obs = _spearman(self.a,self.b)
        self.assertFloatEqual(obs,exp)

        # Both vectors have no ties.
        exp = 0.2969697
        obs = _spearman(self.b,self.c)
        self.assertFloatEqual(obs,exp)

        # Both vectors have ties.
        exp = 0.388381
        obs = _spearman(self.a,self.r)
        self.assertFloatEqual(obs,exp)

    def test_spearman_invalid_input(self):
        """Test the _spearman function with invalid input."""
        self.assertRaises(ValueError, _spearman, [],[])
        self.assertRaises(ValueError, _spearman, self.a,[])
        self.assertRaises(ValueError, _spearman, {0:2}, [1,2,3])

    def test_get_rank(self):
        """Test the _get_rank function with valid input."""
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
        """Test the _get_rank function with invalid input."""
        vec = [1, 'a', 3, 2.5, 3, 1]
        self.assertRaises(TypeError, _get_rank, vec)

        vec = [1, 2, {1:2}, 2.5, 3, 1]
        self.assertRaises(TypeError, _get_rank, vec)

        vec = [1, 2, [23,1], 2.5, 3, 1]
        self.assertRaises(TypeError, _get_rank, vec)

        vec = [1, 2, (1,), 2.5, 3, 1]
        self.assertRaises(TypeError, _get_rank, vec)


if __name__ == "__main__":
    main()
