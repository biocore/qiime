#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Test suite for the compare_taxa_summaries.py module."""

from string import digits
from numpy import array
from unittest import TestCase, main
from numpy.testing import assert_almost_equal
import numpy as np
from itertools import izip
from types import StringType, ListType, FloatType, TupleType

from qiime.compare_taxa_summaries import (compare_taxa_summaries,
                                          _compute_correlation, _make_compatible_taxa_summaries,
                                          _sort_and_fill_taxa_summaries)


class CompareTaxaSummariesTests(TestCase):

    """Tests for the compare_taxa_summaries.py module."""

    def compare_multiple_level_array(self, observed, expected):
        """ Compare multiple level arrays.

        It expecte observed and expected arrays, where each element is an
        array of elements.
        """
        if isinstance(observed, (TupleType, ListType)):
            for obs, exp in izip(observed, expected):
                self.compare_multiple_level_array(obs, exp)
        elif observed is not None and isinstance(observed, (np.number, np.ndarray, FloatType)):
            assert_almost_equal(observed, expected, decimal=5)
        else:
            self.assertEqual(observed, expected)

    def remove_nums(self, text):
        """Removes all digits from the given string.

        Returns the string will all digits removed. Useful for testing strings
        for equality in unit tests where you don't care about numeric values,
        or if some values are random.

        This code was taken from http://bytes.com/topic/python/answers/
            850562-finding-all-numbers-string-replacing

        Arguments:
            text - the string to remove digits from
        """
        return text.translate(None, digits)

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.value_for_seed = 20

        self.taxa_summary1 = (['Even1', 'Even2', 'Even3'],
                              ['Bacteria;Actinobacteria;Actinobacteria(class);Actinobacteridae',
                               'Bacteria;Bacteroidetes/Chlorobigroup;Bacteroidetes;Bacteroidia',
                               'Bacteria;Firmicutes;Bacilli;Lactobacillales',
                               'Bacteria;Firmicutes;Clostridia;Clostridiales',
                               'Bacteria;Firmicutes;Erysipelotrichi;Erysipelotrichales',
                               'Bacteria;Proteobacteria;Gammaproteobacteria;Enterobacteriales',
                               'No blast hit;Other'],
                              array(
                                  [[0.0880247251673, 0.0721968465746, 0.081371761759],
                                   [0.192137761955,
                                    0.191095101593,
                                    0.188504131885],
                                   [0.0264895739603,
                                    0.0259942669171,
                                    0.0318460745596],
                                   [0.491800007824,
                                    0.526186212556,
                                    0.49911159984],
                                   [0.0311411916592,
                                    0.0184083913576,
                                    0.0282325481054],
                                   [0.166137214246,
                                    0.163087129528,
                                    0.168923372865],
                                   [0.00426952518811, 0.00303205147361, 0.0020105109874]]))

        self.taxa_summary2 = (['Even4', 'Even5', 'Even6'],
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
                                     [0.192137761955,
                                      0.191095101593,
                                      0.188504131885],
                                     [0.0264895739603,
                                      0.0259942669171,
                                      0.0318460745596],
                                     [0.491800007824,
                                      0.526186212556,
                                      0.49911159984],
                                     [0.0311411916592,
                                      0.0184083913576,
                                      0.0282325481054],
                                     [0.166137214246,
                                      0.163087129528,
                                      0.168923372865],
                                     [0.00426952518811, 0.00303205147361, 0.0020105109874]]))

        self.taxa_summary3 = (['Even7', 'Even8'], ['Eukarya'],
                              array([[1.0, 1.0]]))

        # Has the same taxa as self.taxa_summary3_data (above).
        self.taxa_summary4 = (['Even1', 'Even2'], ['Eukarya'],
                              array([[0.5, 0.6]]))

        # A sample ID map for testing making compatible taxa summaries.
        self.sample_id_map1 = {'Even7': '1', 'Even1': '1', 'Even8': '2',
                               'Even2': '2'}

        # A sample ID map with an extra sample ID not found in the taxa
        # summaries.
        self.sample_id_map2 = {'Even7': '1', 'Even1': '1', 'Even8': '2',
                               'Even2': '2', 'Foo': '3'}

        # A sample ID map without all original sample IDs mapped.
        self.sample_id_map3 = {'Even7': '1', 'Even1': '1'}

        # No intersection with self.taxa_summary4_data.
        self.taxa_summary5 = (['Even1', 'Even2'], ['foo'],
                              array([[0.5, 0.6]]))

        # Different sample ID from self.taxa_summary5.
        self.taxa_summary6 = (['Even1', 'Even3'], ['foo'],
                              array([[0.1, 0.6]]))

        # Samples are in different orders from self.taxa_summary6.
        self.taxa_summary7 = (['Even3', 'Even1'], ['foo'],
                              array([[0.2, 0.77]]))

        # Samples are not in alphabetical order, and we have multiple taxa.
        self.taxa_summary8 = (['Even3', 'Even1', 'S7'], ['foo', 'bar'],
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
                                     ['Eukarya',
                                      'Bacteria',
                                      'Archaea',
                                      'Foobar'],
                                     array([[0.5, 0.6], [0.7, 0.8], [0.5, 0.6], [0.1, 0.9]]))
        # Full mapping.
        self.taxa_summary_paired_samp_id_map1 = {'S1': 'a', 'S2': 'b',
                                                 'E1': 'a', 'E2': 'b'}
        # Partial mapping.
        self.taxa_summary_paired_samp_id_map2 = {'S1': 'a', 'S2': 'b',
                                                 'E1': 'c', 'E2': 'a'}

    def test_compare_taxa_summaries_expected_pearson(self):
        """Functions correctly using 'expected' comparison mode and pearson."""
        # Overall correlation.
        exp = ('Taxon\tS1\tS2\nBacteria\t0.5\t0.7\nEukarya\t0.4\t0.5\n',
               'Taxon\tExpected\nBacteria\t0.6\nEukarya\t0.5\n',
               '# Correlation coefficient: pearson.\n# The parametric p-value(s) '
               'were calculated using a two-sided test of significance using a '
               't-distribution.\n# The nonparametric p-value(s) were calculated '
               'using a two-sided permutation test with  permutations.\n# The '
               'confidence interval(s) were constructed at a confidence level of '
               '.% using Fisher\'s z-transformation (see Sokal and Rohlf rd '
               'edition pg. ). The confidence interval(s) are two-sided.\n'
               'Correlation coefficient\tParametric p-value\tNonparametric '
               'p-value\tCI (lower)\tCI (upper)\n.\t.\t.\t-.\t.\n', None)

        obs = compare_taxa_summaries(self.taxa_summary_obs1,
                                     self.taxa_summary_exp1, 'expected', 'pearson')

        # We only test the structure of the returned string because of possible
        # roundoff error and stochastic p-values. The actual numbers will be
        # tested in the appropriate 'private' functions.
        obs = (obs[0], obs[1], self.remove_nums(obs[2]), obs[3])
        self.assertEqual(obs, exp)

        # Broken down by sample pairs, with an undefined confidence interval.
        exp = ('Taxon\tS1\tS2\nBacteria\t0.5\t0.7\nEukarya\t0.4\t0.5\n',
               'Taxon\tExpected\nBacteria\t0.6\nEukarya\t0.5\n',
               '# Correlation coefficient: pearson.\n# The parametric p-value(s) '
               'were calculated using a two-sided test of significance using a '
               't-distribution.\n# The nonparametric p-value(s) were calculated '
               'using a two-sided permutation test with  permutations.\n# The '
               'confidence interval(s) were constructed at a confidence level of '
               '.% using Fisher\'s z-transformation (see Sokal and Rohlf rd '
               'edition pg. ). The confidence interval(s) are two-sided.\n'
               'Correlation coefficient\tParametric p-value\tNonparametric '
               'p-value\tCI (lower)\tCI (upper)\n.\t.\t.\t-.\t.\n',
               '# Correlation coefficient: pearson.\n# The parametric p-value(s) '
               'were calculated using a two-sided test of significance using a '
               't-distribution.\n# The nonparametric p-value(s) were calculated '
               'using a two-sided permutation test with  permutations.\n# The '
               'confidence interval(s) were constructed at a confidence level of '
               '.% using Fisher\'s z-transformation (see Sokal and Rohlf rd '
               'edition pg. ). The confidence interval(s) are two-sided.\n'
               'Sample ID\tSample ID\tCorrelation coefficient\t'
               'Parametric p-value\tParametric p-value (Bonferroni-corrected)\t'
               'Nonparametric p-value\tNonparametric p-value '
               '(Bonferroni-corrected)\tCI (lower)\tCI (upper)\nS\tExpected\t'
               '.\t.\t.\t.\t.\tN/A\tN/A\nS\tExpected\t.\t.\t.\t.\t.\tN/A\tN/A\n')
        obs = compare_taxa_summaries(self.taxa_summary_obs1,
                                     self.taxa_summary_exp1, 'expected', 'pearson',
                                     perform_detailed_comparisons=True)
        obs = (obs[0], obs[1], self.remove_nums(obs[2]),
               self.remove_nums(obs[3]))
        self.assertEqual(obs, exp)

    def test_compare_taxa_summaries_spearman_warning(self):
        """Functions correctly using spearman with <= 10 obs."""
        exp = ('Taxon\tS1\tS2\nBacteria\t0.5\t0.7\nEukarya\t0.4\t0.5\n',
               'Taxon\tExpected\nBacteria\t0.6\nEukarya\t0.5\n',
               '# Correlation coefficient: spearman.\n# The parametric p-value(s) '
               'were calculated using a two-sided test of significance using a '
               't-distribution.\n# The nonparametric p-value(s) were calculated '
               'using a two-sided permutation test with  permutations.\n# The '
               'confidence interval(s) were constructed at a confidence level of '
               '.% using Fisher\'s z-transformation (see Sokal and Rohlf rd '
               'edition pg. ). The confidence interval(s) are two-sided.\n'
               '# Since there were  or fewer observations when calculating '
               'Spearman\'s rank correlation coefficient, the parametric p-value '
               'is not accurate when using the t-distribution. Please see '
               'Biometry (Sokal and Rohlf, rd edition) page  for more details.\n'
               'Correlation coefficient\tParametric p-value\tNonparametric '
               'p-value\tCI (lower)\tCI (upper)\n.\t.\t.\t-.\t.\n',
               '# Correlation coefficient: spearman.\n# The parametric p-value(s) '
               'were calculated using a two-sided test of significance using a '
               't-distribution.\n# The nonparametric p-value(s) were calculated '
               'using a two-sided permutation test with  permutations.\n# The '
               'confidence interval(s) were constructed at a confidence level of '
               '.% using Fisher\'s z-transformation (see Sokal and Rohlf rd '
               'edition pg. ). The confidence interval(s) are two-sided.\n'
               '# Since there were  or fewer taxa in the sorted and filled taxa '
               'summary files, the parametric p-values and Bonferroni-corrected '
               'parametric p-values are not accurate when using the '
               't-distribution. Please see Biometry (Sokal and Rohlf, '
               'rd edition) page  for more details.\nSample ID\tSample ID\t'
               'Correlation coefficient\tParametric p-value\tParametric p-value '
               '(Bonferroni-corrected)'
               '\tNonparametric p-value\tNonparametric p-value '
               '(Bonferroni-corrected)\tCI (lower)\tCI (upper)\nS\tExpected\t'
               '.\t.\t.\t.\t.\tN/A\tN/A\nS\tExpected\t.\t.\t.\t.\t.\tN/A\tN/A\n')
        obs = compare_taxa_summaries(self.taxa_summary_obs1,
                                     self.taxa_summary_exp1, 'expected', 'spearman',
                                     perform_detailed_comparisons=True)
        obs = (obs[0], obs[1], self.remove_nums(obs[2]),
               self.remove_nums(obs[3]))
        self.assertEqual(obs, exp)

    def test_compare_taxa_summaries_spearman_no_warning(self):
        """Functions correctly using spearman with > 10 obs."""
        # We should receive a warning only for the individual tests, but not
        # for the overall test.
        exp = ('# Correlation coefficient: spearman.\n# The parametric '
               'p-value(s) were calculated using a two-sided test of significance'
               ' using a '
               't-distribution.\n# The nonparametric p-value(s) were calculated '
               'using a two-sided permutation test with  permutations.\n# The '
               'confidence interval(s) were constructed at a confidence level of '
               '.% using Fisher\'s z-transformation (see Sokal and Rohlf rd '
               'edition pg. ). The confidence interval(s) are two-sided.\n# '
               'Number of samples that matched between the taxa summary files: '
               '\nCorrelation coefficient\tParametric p-value\tNonparametric '
               'p-value\tCI (lower)\tCI (upper)\n.\t.\t.\t-.\t.\n',
               '# Correlation coefficient: spearman.\n# The parametric '
               'p-value(s) were calculated using a two-sided test of significance'
               ' using a '
               't-distribution.\n# The nonparametric p-value(s) were calculated '
               'using a two-sided permutation test with  permutations.\n# The '
               'confidence interval(s) were constructed at a confidence level of '
               '.% using Fisher\'s z-transformation (see Sokal and Rohlf rd '
               'edition pg. ). The confidence interval(s) are two-sided.\n# '
               'Number of samples that matched between the taxa summary files: '
               '\n# Since there were  or fewer taxa in the sorted and filled '
               'taxa summary files, the parametric p-values and '
               'Bonferroni-corrected parametric p-values are not accurate when '
               'using the t-distribution. '
               'Please see Biometry (Sokal and Rohlf, rd edition) page  '
               'for more details.\n'
               'Sample ID\tSample ID\tCorrelation coefficient\t'
               'Parametric p-value\tParametric p-value (Bonferroni-corrected)\t'
               'Nonparametric p-value\tNonparametric p-value '
               '(Bonferroni-corrected)\tCI (lower)\tCI (upper)\n'
               'Even\tEven\t.\t.\t.\t.\t.\t-.\t.\n'
               'Even\tEven\t.\t.\t.\t.\t.\t-.\t.\n'
               'Even\tEven\t.\t.\t.\t.\t.\t-.\t.\n')

        sample_id_map = {'Even1': 'a', 'Even2': 'b', 'Even3': 'c',
                         'Even4': 'a', 'Even5': 'b', 'Even6': 'c'}
        obs = compare_taxa_summaries(self.taxa_summary1,
                                     self.taxa_summary2, 'paired', 'spearman',
                                     perform_detailed_comparisons=True, sample_id_map=sample_id_map)
        obs = (self.remove_nums(obs[2]), self.remove_nums(obs[3]))
        self.assertEqual(obs, exp)

    def test_compare_taxa_summaries_paired_pearson(self):
        """Functions correctly using 'paired' comparison mode and pearson."""
        exp = ('Taxon\tS1\tS2\nArchaea\t0.4\t0.4\nBacteria\t0.5\t'
               '0.7\nEukarya\t0.4\t0.5\n', 'Taxon\tS1\tS2\nArchaea\t0.5\t'
               '0.6\nBacteria\t0.7\t0.8\nEukarya\t0.5\t0.6\n',
               '# Correlation coefficient: pearson.\n# The parametric p-value(s) '
               'were calculated using a two-sided test of significance using a '
               't-distribution.\n# The nonparametric p-value(s) were calculated '
               'using a two-sided permutation test with  permutations.\n# The '
               'confidence interval(s) were constructed at a confidence level of '
               '.% using Fisher\'s z-transformation (see Sokal and Rohlf rd '
               'edition pg. ). The confidence interval(s) are two-sided.\n# '
               'Number of samples that matched between the taxa summary files: '
               '\nCorrelation coefficient\tParametric p-value\tNonparametric '
               'p-value\tCI (lower)\tCI (upper)\n.\t.\t.\t.\t.\n', None)

        obs = compare_taxa_summaries(self.taxa_summary_paired1,
                                     self.taxa_summary_paired2, 'paired', 'pearson')
        obs = (obs[0], obs[1], self.remove_nums(obs[2]),
               obs[3])
        self.assertEqual(obs, exp)

    def test_compare_taxa_summaries_paired_sample_id_map(self):
        """Functions correctly using 'paired' comp mode and sample id map."""
        exp = ('Taxon\tS1\tS2\nArchaea\t0.4\t0.4\nBacteria\t0.5\t'
               '0.7\nEukarya\t0.4\t0.5\n', 'Taxon\tE1\tE2\nArchaea\t0.5\t'
               '0.6\nBacteria\t0.7\t0.8\nEukarya\t0.5\t0.6\n',
               '# Correlation coefficient: pearson.\n# The parametric p-value(s) '
               'were calculated using a two-sided test of significance using a '
               't-distribution.\n# The nonparametric p-value(s) were calculated '
               'using a two-sided permutation test with  permutations.\n# The '
               'confidence interval(s) were constructed at a confidence level of '
               '.% using Fisher\'s z-transformation (see Sokal and Rohlf rd '
               'edition pg. ). The confidence interval(s) are two-sided.\n# '
               'Number of samples that matched between the taxa summary files: '
               '\nCorrelation coefficient\tParametric p-value\tNonparametric '
               'p-value\tCI (lower)\tCI (upper)\n.\t.\t.\t.\t.\n', None)
        obs = compare_taxa_summaries(self.taxa_summary_paired1,
                                     self.taxa_summary_paired4, 'paired', 'pearson',
                                     perform_detailed_comparisons=False,
                                     sample_id_map=self.taxa_summary_paired_samp_id_map1)
        self.assertEqual(obs[0], exp[0])
        self.assertEqual(obs[1], exp[1])
        self.assertEqual(self.remove_nums(obs[2]), exp[2])

    def test_compare_taxa_summaries_paired_sample_id_map_partial(self):
        """Functions correctly using a partial sample id map."""
        # The sample ID map has some mappings that are not complete- i.e. a
        # sample from one file has a new sample ID that doesn't match any other
        # new sample IDs. In this case, the sample should be ignored.
        exp = ('Taxon\tS1\tS2\nArchaea\t0.4\t0.4\nBacteria\t0.5\t'
               '0.7\nEukarya\t0.4\t0.5\n', 'Taxon\tE1\tE2\nArchaea\t0.5'
               '\t0.6\nBacteria\t0.7\t0.8\nEukarya\t0.5\t0.6\n',
               '# Correlation coefficient: pearson.\n# The parametric p-value(s) '
               'were calculated using a two-sided test of significance using a '
               't-distribution.\n# The '
               'confidence interval(s) were constructed at a confidence level of '
               '95.0% using Fisher\'s z-transformation (see Sokal and Rohlf 3rd '
               'edition pg. 575). The confidence interval(s) are two-sided.\n# '
               'Number of samples that matched between the taxa summary files: 1'
               '\nCorrelation coefficient\tParametric p-value\tNonparametric '
               'p-value\tCI (lower)\tCI (upper)\n1.0000\t0.0000\tN/A\tN/A\tN/A\n',
               '# Correlation coefficient: pearson.\n# The parametric p-value(s) '
               'were calculated using a two-sided test of significance using a '
               't-distribution.\n# The '
               'confidence interval(s) were constructed at a confidence level of '
               '95.0% using Fisher\'s z-transformation (see Sokal and Rohlf 3rd '
               'edition pg. 575). The confidence interval(s) are two-sided.\n# '
               'Number of samples that matched between the taxa summary files: 1'
               '\nSample ID\tSample ID\tCorrelation coefficient\tParametric '
               'p-value\tParametric p-value (Bonferroni-corrected)\t'
               'Nonparametric p-value\tNonparametric p-value '
               '(Bonferroni-corrected)\tCI (lower)\tCI (upper)\nS1\tE2\t1.0000\t'
               '0.0000\t0.0000\tN/A\tN/A\tN/A\tN/A\n')

        obs = compare_taxa_summaries(self.taxa_summary_paired1,
                                     self.taxa_summary_paired4, 'paired', 'pearson',
                                     num_permutations=0,
                                     perform_detailed_comparisons=True,
                                     sample_id_map=self.taxa_summary_paired_samp_id_map2)
        # We can test exactly because there aren't any stochastic p-values.
        self.assertEqual(obs, exp)

    def test_compare_taxa_summaries_paired_sample_id_map_mismatched_taxa(self):
        """Functions correctly using a sample id map and mismatched taxa."""
        exp = ('Taxon\tS1\tS2\nArchaea\t0.4\t0.4\nBacteria\t0.5\t0.7\nEukarya'
               '\t0.4\t0.5\nFoobar\t0.0\t0.0\n', 'Taxon\tE1\tE2\nArchaea\t0.5'
               '\t0.6\nBacteria\t0.7\t0.8\nEukarya\t0.5\t0.6\nFoobar\t0.1\t0.9'
               '\n',
               '# Correlation coefficient: pearson.\n# The parametric p-value(s) '
               'were calculated using a two-sided test of significance using a '
               't-distribution.\n# The nonparametric p-value(s) were calculated '
               'using a two-sided permutation test with  permutations.\n# The '
               'confidence interval(s) were constructed at a confidence level of '
               '.% using Fisher\'s z-transformation (see Sokal and Rohlf rd '
               'edition pg. ). The confidence interval(s) are two-sided.\n# '
               'Number of samples that matched between the taxa summary files: '
               '\nCorrelation coefficient\tParametric p-value\tNonparametric '
               'p-value\tCI (lower)\tCI (upper)\n-.\t.\t.\t-.\t.\n', None)

        obs = compare_taxa_summaries(self.taxa_summary_paired1,
                                     self.taxa_summary_paired5, 'paired', 'pearson',
                                     perform_detailed_comparisons=False,
                                     sample_id_map=self.taxa_summary_paired_samp_id_map2)
        obs = (obs[0], obs[1], self.remove_nums(obs[2]), obs[3])
        self.assertEqual(obs, exp)

    def test_compare_taxa_summaries_tail_type(self):
        """Functions correctly using various tail types."""
        # High.
        exp = ('# Correlation coefficient: pearson.\n# The parametric '
               'p-value(s) were calculated using a one-sided (positive '
               'association) test of significance using a '
               't-distribution.\n# The nonparametric p-value(s) were calculated '
               'using a one-sided (positive association) permutation test with  '
               'permutations.\n# The '
               'confidence interval(s) were constructed at a confidence level of '
               '.% using Fisher\'s z-transformation (see Sokal and Rohlf rd '
               'edition pg. ). The confidence interval(s) are two-sided.\n# '
               'Number of samples that matched between the taxa summary files: '
               '\nCorrelation coefficient\tParametric p-value\tNonparametric '
               'p-value\tCI (lower)\tCI (upper)\n-.\t.\t.\t-.\t.\n')

        obs = compare_taxa_summaries(self.taxa_summary_paired1,
                                     self.taxa_summary_paired5, 'paired', 'pearson',
                                     tail_type='high', confidence_level=0.80,
                                     perform_detailed_comparisons=False,
                                     sample_id_map=self.taxa_summary_paired_samp_id_map2)
        obs = (self.remove_nums(obs[2]))
        self.assertEqual(obs, exp)

        # Low.
        exp = ('# Correlation coefficient: pearson.\n# The parametric '
               'p-value(s) were calculated using a one-sided (negative '
               'association) test of significance using a '
               't-distribution.\n# The nonparametric p-value(s) were calculated '
               'using a one-sided (negative association) permutation test with  '
               'permutations.\n# The '
               'confidence interval(s) were constructed at a confidence level of '
               '.% using Fisher\'s z-transformation (see Sokal and Rohlf rd '
               'edition pg. ). The confidence interval(s) are two-sided.\n# '
               'Number of samples that matched between the taxa summary files: '
               '\nCorrelation coefficient\tParametric p-value\tNonparametric '
               'p-value\tCI (lower)\tCI (upper)\n-.\t.\t.\t-.\t.\n')

        obs = compare_taxa_summaries(self.taxa_summary_paired1,
                                     self.taxa_summary_paired5, 'paired', 'pearson',
                                     tail_type='low', confidence_level=0.80,
                                     perform_detailed_comparisons=False,
                                     sample_id_map=self.taxa_summary_paired_samp_id_map2)
        obs = (self.remove_nums(obs[2]))
        self.assertEqual(obs, exp)

    def test_compare_taxa_summaries_no_permutations(self):
        """Functions correctly without performing nonparametric test."""
        # Verified with R's cor.test function.
        exp = ('# Correlation coefficient: pearson.\n# The parametric '
               'p-value(s) were calculated using a one-sided (positive '
               'association) test of significance using a '
               't-distribution.\n# The '
               'confidence interval(s) were constructed at a confidence level of '
               '70.0% using Fisher\'s z-transformation (see Sokal and Rohlf 3rd '
               'edition pg. 575). The confidence interval(s) are two-sided.\n# '
               'Number of samples that matched between the taxa summary files: 1'
               '\nCorrelation coefficient\tParametric p-value\tNonparametric '
               'p-value\tCI (lower)\tCI (upper)\n-0.6264\t0.8132\tN/A\t-0.9438\t'
               '0.2922\n')

        obs = compare_taxa_summaries(self.taxa_summary_paired1,
                                     self.taxa_summary_paired5, 'paired', 'pearson',
                                     tail_type='high', num_permutations=0, confidence_level=0.70,
                                     perform_detailed_comparisons=False,
                                     sample_id_map=self.taxa_summary_paired_samp_id_map2)
        obs = (obs[2])
        self.assertEqual(obs, exp)

    def test_compare_taxa_summaries_correct_header(self):
        """Header is constructed properly based on different parameters."""
        exp = ('# Correlation coefficient: pearson.\n# The parametric '
               'p-value(s) were calculated using a one-sided (positive '
               'association) test of significance using a t-distribution.\n# The '
               'confidence interval(s) were constructed at a confidence level of '
               '90.0% using Fisher\'s z-transformation (see Sokal and Rohlf 3rd '
               'edition pg. 575). The confidence interval(s) are two-sided.\n# '
               'Number of samples that matched between the taxa summary files: 1'
               '\nCorrelation coefficient\tParametric p-value\tNonparametric '
               'p-value\tCI (lower)\tCI (upper)')

        obs = compare_taxa_summaries(self.taxa_summary_paired1,
                                     self.taxa_summary_paired5, 'paired', 'pearson',
                                     tail_type='high', num_permutations=0, confidence_level=0.90,
                                     perform_detailed_comparisons=False,
                                     sample_id_map=self.taxa_summary_paired_samp_id_map2)

        # Only look at the header.
        obs = '\n'.join(obs[2].split('\n')[:-2])
        self.assertEqual(obs, exp)

        exp = ('# Correlation coefficient: pearson.\n# The parametric '
               'p-value(s) were calculated using a one-sided (negative '
               'association) test of significance using a t-distribution.\n'
               '# The nonparametric p-value(s) were calculated '
               'using a one-sided (negative association) permutation test with '
               '85 permutations.\n# The '
               'confidence interval(s) were constructed at a confidence level of '
               '5.668% using Fisher\'s z-transformation (see Sokal and Rohlf 3rd '
               'edition pg. 575). The confidence interval(s) are two-sided.\n# '
               'Number of samples that matched between the taxa summary files: 1'
               '\nCorrelation coefficient\tParametric p-value\tNonparametric '
               'p-value\tCI (lower)\tCI (upper)')

        obs = compare_taxa_summaries(self.taxa_summary_paired1,
                                     self.taxa_summary_paired5, 'paired', 'pearson',
                                     tail_type='low', num_permutations=85, confidence_level=0.05668,
                                     perform_detailed_comparisons=False,
                                     sample_id_map=self.taxa_summary_paired_samp_id_map2)

        # Only look at the header.
        obs = '\n'.join(obs[2].split('\n')[:-2])
        self.assertEqual(obs, exp)

    def test_compare_taxa_summaries_invalid_input(self):
        """Throws error on invalid input."""
        # Invalid comparison mode.
        self.assertRaises(ValueError, compare_taxa_summaries,
                          self.taxa_summary_obs1, self.taxa_summary_exp1, 'foo',
                          'pearson')
        # Invalid correlation type.
        self.assertRaises(ValueError, compare_taxa_summaries,
                          self.taxa_summary_obs1, self.taxa_summary_exp1, 'paired',
                          'foo')
        # Invalid tail type.
        self.assertRaises(ValueError, compare_taxa_summaries,
                          self.taxa_summary_obs1, self.taxa_summary_exp1, 'paired',
                          'spearman', 'foo')
        # Invalid number of permutations.
        self.assertRaises(ValueError, compare_taxa_summaries,
                          self.taxa_summary_obs1, self.taxa_summary_exp1, 'paired',
                          'spearman', 'high', -1)
        # Invalid confidence level.
        self.assertRaises(ValueError, compare_taxa_summaries,
                          self.taxa_summary_obs1, self.taxa_summary_exp1, 'paired',
                          'spearman', 'high', 0, 1)

    def test_make_compatible_taxa_summaries(self):
        """Test making compatible taxa summaries works correctly on two ts."""
        exp = [(['Even1'], ['foo'], array([[0.5]])), (['Even1'], ['foo'],
                                                      array([[0.1]]))]
        obs = _make_compatible_taxa_summaries(self.taxa_summary5,
                                              self.taxa_summary6)
        self.compare_multiple_level_array(obs, exp)

    def test_make_compatible_taxa_summaries_already_compatible(self):
        """Test on taxa summaries that are already compatible."""
        # Using two compatible ts.
        obs = _make_compatible_taxa_summaries(self.taxa_summary4,
                                              self.taxa_summary5)
        self.compare_multiple_level_array(obs, (self.taxa_summary4, self.taxa_summary5))

        # Using the same ts twice.
        obs = _make_compatible_taxa_summaries(self.taxa_summary4,
                                              self.taxa_summary4)
        self.compare_multiple_level_array(obs, (self.taxa_summary4, self.taxa_summary4))

        obs = _make_compatible_taxa_summaries(self.taxa_summary5,
                                              self.taxa_summary5)
        self.compare_multiple_level_array(obs, (self.taxa_summary5, self.taxa_summary5))

    def test_make_compatible_taxa_summaries_unordered(self):
        """Test on taxa summaries whose sample IDs are in different orders."""
        # Using two compatible ts.
        exp = [(['Even1', 'Even3'], ['foo'], array([[0.1, 0.6]])),
               (['Even1', 'Even3'], ['foo'], array([[0.77, 0.2]]))]
        obs = _make_compatible_taxa_summaries(self.taxa_summary6,
                                              self.taxa_summary7)
        self.compare_multiple_level_array(obs, exp)

    def test_make_compatible_taxa_summaries_sample_id_map(self):
        """Test making compatible ts using a sample ID map."""
        exp = ((['Even7', 'Even8'], ['Eukarya'], array([[1., 1.]])),
               (['Even1', 'Even2'], ['Eukarya'], array([[0.5, 0.6]])))
        obs = _make_compatible_taxa_summaries(self.taxa_summary3,
                                              self.taxa_summary4, self.sample_id_map1)
        self.compare_multiple_level_array(obs, exp)

    def test_make_compatible_taxa_summaries_sample_id_map_extra_sample(self):
        """Test using a sample ID map with an extra sample ID."""
        exp = ((['Even7', 'Even8'], ['Eukarya'], array([[1., 1.]])),
               (['Even1', 'Even2'], ['Eukarya'], array([[0.5, 0.6]])))
        obs = _make_compatible_taxa_summaries(self.taxa_summary3,
                                              self.taxa_summary4, self.sample_id_map2)
        self.compare_multiple_level_array(obs, exp)

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
            (['Even1', 'Even2', 'Even3'],
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
                    [0., 0., 0.],
                    [0., 0., 0.],
                    [0.192137761955, 0.191095101593, 0.188504131885],
                    [0.0264895739603, 0.0259942669171, 0.0318460745596],
                    [0.491800007824, 0.526186212556, 0.49911159984],
                    [0.0311411916592, 0.0184083913576, 0.0282325481054],
                    [0.166137214246, 0.163087129528, 0.168923372865],
                    [0., 0., 0.],
                    [0.00426952518811, 0.00303205147361, 0.0020105109874]])),
            (['Even4', 'Even5', 'Even6'],
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
             array([[0., 0., 0.],
                    [0.99, 0.11, 0.075],
                    [0.1921, 0.19109, 0.18],
                    [0.192137761955, 0.191095101593, 0.188504131885],
                    [0.0264895739603, 0.0259942669171, 0.0318460745596],
                    [0.491800007824, 0.526186212556, 0.49911159984],
                    [0.0311411916592, 0.0184083913576, 0.0282325481054],
                    [0.166137214246, 0.163087129528, 0.168923372865],
                    [0., 0., 0.],
                    [0.00426952518811, 0.00303205147361, 0.0020105109874]])),
            (['Even7', 'Even8'],
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
             array([[0., 0.],
                    [0., 0.],
                    [0., 0.],
                    [0., 0.],
                    [0., 0.],
                    [0., 0.],
                    [0., 0.],
                    [0., 0.],
                    [1., 1.],
                    [0., 0.]]))
        ]

        obs = _sort_and_fill_taxa_summaries([self.taxa_summary1,
                                             self.taxa_summary2,
                                             self.taxa_summary3])
        self.compare_multiple_level_array(obs, exp)

    def test_sort_and_fill_taxa_summaries_same(self):
        """Test when the taxa summaries are the same."""
        exp = [(['Even7', 'Even8'], ['Eukarya'], array([[1.0, 1.0]])),
               (['Even7', 'Even8'], ['Eukarya'], array([[1.0, 1.0]]))]
        obs = _sort_and_fill_taxa_summaries([self.taxa_summary3,
                                             self.taxa_summary3])
        self.compare_multiple_level_array(obs, exp)

        exp = [(['Even7', 'Even8'], ['Eukarya'], array([[1.0, 1.0]])),
               (['Even1', 'Even2'], ['Eukarya'], array([[0.5, 0.6]]))]
        obs = _sort_and_fill_taxa_summaries([self.taxa_summary3,
                                             self.taxa_summary4])
        self.compare_multiple_level_array(obs, exp)

        # Test the other direction.
        exp = [(['Even1', 'Even2'], ['Eukarya'], array([[0.5, 0.6]])),
               (['Even7', 'Even8'], ['Eukarya'], array([[1.0, 1.0]]))]
        obs = _sort_and_fill_taxa_summaries([self.taxa_summary4,
                                             self.taxa_summary3])
        self.compare_multiple_level_array(obs, exp)

    def test_sort_and_fill_taxa_summaries_no_intersection(self):
        """Test when the taxa summaries have no intersection."""
        exp = [(['Even1', 'Even2'], ['Eukarya', 'foo'], array([[0.5, 0.6],
                                                              [0.0, 0.0]])),
               (['Even1', 'Even2'], ['Eukarya', 'foo'], array([[0.0, 0.0],
                                                              [0.5, 0.6]]))]
        obs = _sort_and_fill_taxa_summaries([self.taxa_summary4,
                                             self.taxa_summary5])
        self.compare_multiple_level_array(obs, exp)

    def test_compute_correlation_expected_pearson(self):
        """Test functions correctly with expected mode and pearson corr."""
        exp = ((0.68824720161169595, 0.31175279838830405, 0.689,
               (-0.80594408245459292, 0.99269848760560575)), None)
        np.random.seed(self.value_for_seed)
        obs = _compute_correlation(self.taxa_summary_obs1,
                                   self.taxa_summary_exp1, 'expected', 'pearson', 'two-sided',
                                   999, 0.95)
        self.compare_multiple_level_array(obs, exp)

    def test_compute_correlation_expected_pearson_detailed(self):
        """Test functions with expected mode, pearson corr, detailed tests."""
        exp = ((0.68824720161169595, 0.31175279838830405, 0.664,
                (-0.66416860387615351, 0.9863313909937903)),
               [('S1', 'Expected', 1.0, 1, 1, 1.0, 1, (None, None)),
                ('S2', 'Expected', 1.0, 1, 1, 1.0, 1, (None, None))])
        np.random.seed(self.value_for_seed)
        obs = _compute_correlation(self.taxa_summary_obs1,
                                   self.taxa_summary_exp1, 'expected', 'pearson', 'two-sided',
                                   999, 0.90, True)
        self.compare_multiple_level_array(obs, exp)

    def test_compute_correlation_expected_spearman_detailed(self):
        """Test functions with expected mode, spearman corr, detailed tests."""
        # Verified with R's cor.test function.

        # No repeats in ranks.
        exp = ((0.70710678118654757, 0.29289321881345232, 0.664,
                (-0.93471237439516763, 0.99801522848859603)),
               [('S1', 'Expected', 1.0, 1, 1, 1.0, 1, (None, None)),
                ('S2', 'Expected', 1.0, 1, 1, 1.0, 1, (None, None))])
        np.random.seed(self.value_for_seed)
        obs = _compute_correlation(self.taxa_summary_obs1,
                                   self.taxa_summary_exp1,
                                   'expected', 'spearman', 'two-sided',
                                   999, 0.99, True)
        self.compare_multiple_level_array(obs, exp)

        # Repeats in ranks.
        exp = ((0.83914639167827365, 0.036729871219315036, 0.137,
                (-0.2625774054977657, 0.99110427297276837)),
               [('S1', 'Expected', 0.86602540378443871, 0.33333333333333326,
                 0.66666666666666652, 0.679, 1, (None, None)),
                ('S2', 'Expected', 1.0, 0, 0, 0.327, 0.654, (None, None))])
        np.random.seed(self.value_for_seed)
        obs = _compute_correlation(self.taxa_summary_obs2,
                                   self.taxa_summary_exp2,
                                   'expected', 'spearman', 'two-sided', 999,
                                   0.99, True)
        self.compare_multiple_level_array(obs, exp)

    def test_compute_correlation_expected_expected_sample_id(self):
        """Test functions with expected mode using an expected sample ID."""
        # Using a single-sample expected ts.
        exp = ((0.83914639167827365, 0.036729,
                0.13786213786213786, (0.032537093928499863,
                                      0.98380431996767537)),
               [('S1', 'Expected', 0.86602540378443871, 0.33333333333333326,
                 0.66666666666666652, 0.6793206793206793, 1, (None, None)),
                ('S2', 'Expected', 1.0, 0, 0, 0.32667332667332666,
                 0.6533466533466533, (None, None))])
        np.random.seed(self.value_for_seed)
        obs = _compute_correlation(self.taxa_summary_obs2,
                                   self.taxa_summary_exp2,
                                   'expected', 'spearman', 'two-sided',
                                   1000, 0.96, True, 'Expected')
        self.compare_multiple_level_array(obs, exp)
        
        # Using a two-sample expected ts.
        exp = ((0.83914639167827365, 0.036729,
                0.13786213786213786, (0.032537093928499863,
                                      0.98380431996767537)),
               [('S1', 'Expected', 0.86602540378443871, 0.33333333333333326,
                 0.66666666666666652, 0.6793206793206793, 1, (None, None)),
                ('S2', 'Expected', 1.0, 0, 0, 0.32667332667332666,
                 0.6533466533466533, (None, None))])

        np.random.seed(self.value_for_seed)
        obs = _compute_correlation(self.taxa_summary_obs2,
                                   self.taxa_summary_exp3,
                                   'expected', 'spearman', 'two-sided',
                                   1000, 0.96, True, 'Expected')
        self.compare_multiple_level_array(obs, exp)

    def test_compute_correlation_paired(self):
        """Test runs correctly on standard input taxa summaries."""
        # Verified using R's cor.test function.
        exp = ((0.90243902439024193, 0.013812916237431808, 0.03,
               (0.33958335414859975, 0.98938788428012969)), None)
        
        np.random.seed(self.value_for_seed)
        obs = _compute_correlation(self.taxa_summary_paired1,
                                   self.taxa_summary_paired2, 'paired', 'pearson', 'two-sided',
                                   999, 0.95)
        self.compare_multiple_level_array(obs, exp)

    def test_compute_correlation_paired_detailed(self):
        """Test runs correctly on standard input ts with detailed tests."""
        # Verified using R's cor.test function.
        exp = ((0.90243902439024193, 0.013812916237431808, 0.038,
                (0.33958335414859975, 0.98938788428012969)),
               [('S1', 'S1', 1.0, 0, 0, 0.373, 0.746, (None, None)),
                ('S2', 'S2', 0.94491118252306538, 0.21229561500966185,
                 0.4245912300193237, 0.342, 0.684, (None, None))])

        np.random.seed(self.value_for_seed)
        obs = _compute_correlation(self.taxa_summary_paired1,
                                   self.taxa_summary_paired2, 'paired', 'pearson', 'two-sided',
                                   999, 0.95, True)

        self.compare_multiple_level_array(obs, exp)

    def test_compute_correlation_paired_mismatched_ids(self):
        """Test runs correctly on taxa summaries with mismatched IDs."""
        # Verified using R's cor.test function.
        exp = ((0.90243902439024193, 0.013812916237431808, 0.038,
                (0.48961251524616184, 0.98476601971494115)),
               [('S1', 'E1', 1.0, 0, 0, 0.373, 0.746, (None, None)),
                ('S2', 'E1', 0.94491118252306538, 0.21229561500966185,
                 0.4245912300193237, 0.342, 0.684, (None, None))])
        np.random.seed(self.value_for_seed)
        obs = _compute_correlation(self.taxa_summary_paired1,
                                   self.taxa_summary_paired3, 'paired', 'pearson', 'two-sided',
                                   999, 0.90, True)

        self.compare_multiple_level_array(obs, exp)

    def test_compute_correlation_paired_no_permutations(self):
        """Test runs correctly on taxa summaries without nonparametric test."""
        # Verified using R's cor.test function.
        exp = ((0.90243902439024193, 0.013812916237431808, None,
                (0.48961251524616184, 0.98476601971494115)),
               [('S1', 'E1', 1.0, 0, 0, None, None, (None, None)),
                ('S2', 'E1', 0.94491118252306538, 0.21229561500966185,
                 0.4245912300193237, None, None, (None, None))])
        obs = _compute_correlation(self.taxa_summary_paired1,
                                   self.taxa_summary_paired3, 'paired', 'pearson', 'two-sided',
                                   0, 0.90, True)
        self.compare_multiple_level_array(obs, exp)

    def test_compute_correlation_paired_tail_type(self):
        """Test runs correctly on taxa summaries with single-tailed test."""
        # Verified using R's cor.test function.
        exp = ((0.50151319303934072, 0.99615167469303467, None,
               (0.21229308395242716, 0.70994854785747563)),
               [('Even1', 'Even4', 0.18206274970805661, 0.68040119701955393,
                1, None, None, (-0.45214511993120049, 0.69399617731638885)),
                ('Even2', 'Even5', 0.8953176947428052, 0.99944917346644746, 1,
                 None, None, (0.65074921580251011, 0.97157247969474569)),
                ('Even3', 'Even6', 0.89969103596277555, 0.99952349856557832,
                 1, None, None, (0.6635260751754386, 0.97280581291498336))])
        ts1, ts2 = _sort_and_fill_taxa_summaries(
            [self.taxa_summary1, self.taxa_summary2])
        obs = _compute_correlation(ts1, ts2, 'paired', 'pearson', 'low', 0,
                                   0.90, True)
        self.compare_multiple_level_array(obs, exp)

    def test_compute_correlation_expected_invalid_sample_count(self):
        """Test running on expected taxa summary without exactly one sample."""
        self.assertRaises(ValueError, _compute_correlation,
                          self.taxa_summary1, self.taxa_summary2, 'expected',
                          'pearson', 'high', 0, 0.95)

    def test_compute_correlation_expected_invalid_expected_sample_id(self):
        """Test running in expected mode with invalid expected sample ID."""
        self.assertRaises(ValueError, _compute_correlation,
                          self.taxa_summary_obs2, self.taxa_summary_exp3, 'expected',
                          'pearson', 'low', 10, 0.1, expected_sample_id='foo')

    def test_compute_correlation_invalid_comparison_mode(self):
        """Test running using an invalid comparison mode."""
        self.assertRaises(ValueError, _compute_correlation,
                          self.taxa_summary1, self.taxa_summary2, 'foo',
                          'pearson', 'two-sided', 999, 0.90)

    def test_compute_correlation_incompatible_taxa(self):
        """Test running on taxa summaries that have mismatched taxa info."""
        self.assertRaises(ValueError, _compute_correlation,
                          self.taxa_summary_obs1_mismatch,
                          self.taxa_summary_exp1,
                          'expected', 'pearson', 'low', 999, 0.5)

    def test_compute_correlation_paired_incompatible_samples(self):
        """Test on incompatible taxa summaries (mismatched sample lengths)."""
        self.assertRaises(ValueError, _compute_correlation,
                          self.taxa_summary1, self.taxa_summary3, 'paired',
                          'spearman', 'high', 9, 0.22222)

    def test_compute_correlation_invalid_num_permutations(self):
        """Test on invalid number of permutations."""
        self.assertRaises(ValueError, _compute_correlation,
                          self.taxa_summary1, self.taxa_summary1, 'paired',
                          'spearman', 'high', -10, 0.22222)

    def test_compute_correlation_invalid_confidence_level(self):
        """Test on invalid confidence level."""
        self.assertRaises(ValueError, _compute_correlation,
                          self.taxa_summary1, self.taxa_summary1, 'paired',
                          'spearman', 'high', 10, 0)

    def test_compute_correlation_invalid_tail_type(self):
        """Test on invalid tail type."""
        self.assertRaises(ValueError, _compute_correlation,
                          self.taxa_summary1, self.taxa_summary1, 'paired',
                          'spearman', 'foo', 10, 0.1)


if __name__ == "__main__":
    main()
