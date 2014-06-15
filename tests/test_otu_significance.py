#!/usr/bin/env python
# File created on 16 Aug 2013
from __future__ import division

__author__ = "Luke Ursell"
__copyright__ = "Copyright 2013, The QIIME project"
__credits__ = ["Luke Ursell", "Will Van Treuren", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Luke Ursell"
__email__ = "lkursell@gmail.com"

from unittest import TestCase, main
from qiime.otu_significance import (get_sample_cats, get_cat_sample_groups,
                                    get_sample_indices, group_significance_row_generator, sort_by_pval,
                                    run_group_significance_test, group_significance_output_formatter,
                                    GROUP_TEST_CHOICES, grouped_correlation_row_generator,
                                    run_grouped_correlation, CORRELATION_TEST_CHOICES,
                                    grouped_correlation_formatter, correlation_row_generator,
                                    run_correlation_test)
from qiime.stats import (assign_correlation_pval, fisher,
                                   fisher_population_correlation)
from numpy import array, hstack, corrcoef
from numpy.random import seed
from numpy.testing import assert_almost_equal
from os import remove
from qiime.parse import parse_mapping_file_to_dict, parse_otu_table
from biom.parse import parse_biom_table

class GroupSignificanceFunctionsTests(TestCase):

    """Tests of group significance functions."""

    def setUp(self):
        """Define values used by all tests."""
        # nothing to do, all tests use different things
        pass

    def test_get_sample_cats(self):
        """Test get_sample_cats."""
        pmf_in = \
            {'Sample1': {'test_cat': 'cat1', 'test_corr': '1', 'test_empty': 'abc'},
             'Sample2':
             {'test_cat': 'cat1', 'test_corr': '1', 'test_empty': 'abc'},
                'Sample3':
                {'test_cat': 'cat2', 'test_corr': '1', 'test_empty': ''},
                'Sample4':
                {'test_cat': 'cat2', 'test_corr': '1', 'test_empty': ''},
                'Sample5':
                {'test_cat': 'cat3', 'test_corr': '1', 'test_empty': 'abc'},
                'Sample6': {'test_cat': 'cat3', 'test_corr': '1', 'test_empty': ''}}
        # test with properly formatted, ordered parsed mapping file
        category = 'test_cat'
        exp = {'Sample1': 'cat1',
               'Sample2': 'cat1',
               'Sample3': 'cat2',
               'Sample4': 'cat2',
               'Sample5': 'cat3',
               'Sample6': 'cat3'}
        obs = get_sample_cats(pmf_in, category)
        self.assertEqual(exp, obs)
        # test with only one category
        category = 'test_corr'
        exp = {'Sample1': '1',
               'Sample2': '1',
               'Sample3': '1',
               'Sample4': '1',
               'Sample5': '1',
               'Sample6': '1'}
        obs = get_sample_cats(pmf_in, category)
        self.assertEqual(exp, obs)
        # test with empty categories
        category = 'test_empty'
        exp = {'Sample1': 'abc',
               'Sample2': 'abc',
               'Sample5': 'abc'}
        obs = get_sample_cats(pmf_in, category)
        self.assertEqual(exp, obs)

    def test_get_cat_sample_groups(self):
        """Test get_cat_sample_groups works."""
        sample_cats = {'Sample1': 'cat1',
                       'Sample2': 'cat1',
                       'Sample3': 'cat2',
                       'Sample4': 'cat2',
                       'Sample5': 'cat3',
                       'Sample6': 'cat3'}
        exp = {'cat1': ['Sample1', 'Sample2'],
               'cat2': ['Sample4', 'Sample3'],
               'cat3': ['Sample5', 'Sample6']}
        obs = get_cat_sample_groups(sample_cats)
        self.assertEqual(exp, obs)

    def test_get_sample_indices(self):
        """Test get_sample_indices works"""
        cat_sample_groups = {'cat1': ['Sample1', 'Sample2'],
                             'cat2': ['Sample4', 'Sample3'],
                             'cat3': ['Sample5', 'Sample6']}
        bt = parse_biom_table(BT_IN_1)
        exp = {'cat1': [0, 1], 'cat2': [3, 2], 'cat3': [4, 5]}
        obs = get_sample_indices(cat_sample_groups, bt)
        self.assertEqual(exp, obs)

    def test_group_significance_row_generator(self):
        """Test group_significance_row_generator works."""
        # run with ordered example
        sample_indices = {'cat1': [0, 1], 'cat2': [3, 2], 'cat3': [4, 5]}
        bt = parse_biom_table(BT_IN_1)
        data = array([bt.data(i, axis='observation') for i in bt.observation_ids])
        obs = list(group_significance_row_generator(bt, sample_indices))
        exp = zip(data.take([0, 1], 1),
                  data.take([3, 2], 1), data.take([4, 5], 1))
        for o, e in zip(obs, exp):
            assert_almost_equal(e, o)
        # run with unequal length example
        sample_indices = {'g0': [0, 1, 4, 5], 'g1': [3], 'g2': [2]}
        obs = list(group_significance_row_generator(bt, sample_indices))
        exp = zip(data.take([0, 1, 4, 5], 1),
                  data.take([3], 1), data.take([2], 1))
        for o, e in zip(obs, exp):
            [assert_almost_equal(i, j) for i, j in
                zip(sorted(hstack(e)), sorted(hstack(o)))]
            # the order of e and o are not the same. this isn't a problem in
            # the actual function because we iterate over the order of the
            # sample keys every time.

    def test_run_group_significance_test(self):
        """Test that all group significance tests can be run."""
        bt = parse_biom_table(BT_IN_1)
        bt_4 = parse_biom_table(BT_4)

        # test with non-paramteric t-test
        sample_indices = {'cat1': [0, 5, 1], 'cat2': [2, 4, 3]}
        row_gen = group_significance_row_generator(bt, sample_indices)
        exp_test_stats = [0.17503798979747345, 0.20029818620053824,
                          -1.5065313062753816, -
                          0.043884559904114794, -1.0631239617935129,
                          -1.2878361428003895]
        # we are expecting 1001 comparisons)
        exp_pvals = map(lambda x: x / 1001., [888, 899, 279, 1001, 489, 299])
        exp_means = [[52.333333333333336, 48.333333333333336],
                     [34.0, 30.333333333333332],
                     [20.0, 49.333333333333336],
                     [55.333333333333336, 56.0],
                     [20.0, 38.0],
                     [30.0, 60.333333333333336]]
        seed(0)  # seed prng for reproducibility
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'nonparametric_t_test',
                                        GROUP_TEST_CHOICES, reps=1000)
        assert_almost_equal(exp_test_stats, obs_test_stats)
        assert_almost_equal(exp_pvals, obs_pvals)
        assert_almost_equal(exp_means, obs_means)
        # test with nonparametric t-test but different ordering
        sample_indices = {'cat1': [0, 1, 5], 'cat2': [4, 3, 2]}
        row_gen = group_significance_row_generator(bt, sample_indices)
        seed(0)  # seed prng for reproducibility
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'nonparametric_t_test',
                                        GROUP_TEST_CHOICES, reps=1000)
        assert_almost_equal(exp_test_stats, obs_test_stats)
        assert_almost_equal(exp_pvals, obs_pvals)
        assert_almost_equal(exp_means, obs_means)
        # test with BT_4 biom table
        sample_indices = {'cat1': [0, 3, 1, 4], 'cat2': [5, 2, 7, 6]}
        row_gen = group_significance_row_generator(bt_4, sample_indices)
        exp_test_stats = [-0.38741397129147953, -0.38334158591463874,
                          0.077468274988510541, -
                          0.2322539745918096, 0.16469600468808282,
                          -0.49589486133213057]
        # we are expecting 1001 comparisons)
        exp_pvals = map(lambda x: x / 1001., [821, 719, 916, 935, 938, 604])
        exp_means = [[43.5, 51.75],
                     [29.75, 34.75],
                     [41.5, 40.0],
                     [50.5, 53.75],
                     [28.0, 25.5],
                     [41.75, 54.0]]
        seed(0)
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'nonparametric_t_test',
                                        GROUP_TEST_CHOICES, reps=1000)
        assert_almost_equal(exp_test_stats, obs_test_stats)
        assert_almost_equal(exp_pvals, obs_pvals)
        assert_almost_equal(exp_means, obs_means)

        # test with parametric t test
        # bt_1 agrees with Prism
        sample_indices = {'cat1': [4, 1, 2], 'cat2': [5, 0, 3]}
        row_gen = group_significance_row_generator(bt, sample_indices)
        exp_test_stats = [-1.0504514628777806, -0.94113003446934629,
                          -
                          0.66264262463016887, 0.17617555832772411, 1.1144416530351877,
                          -1.2483315640812607]
        exp_pvals = [0.3527834167236007, 0.39992473225679626,
                     0.5437923932346147, 0.8687158192049661, 0.32753202812350557,
                     0.27998887149482976]
        exp_means = [[39.666666666666664, 61.0],
                     [24.333333333333332, 40.0],
                     [27.0, 42.333333333333336],
                     [57.0, 54.333333333333336],
                     [38.333333333333336, 19.666666666666668],
                     [30.333333333333332, 60.0]]
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'parametric_t_test',
                                        GROUP_TEST_CHOICES)
        assert_almost_equal(exp_test_stats, obs_test_stats)
        assert_almost_equal(exp_pvals, obs_pvals)
        assert_almost_equal(exp_means, obs_means)
        # test with BT_4
        sample_indices = {'cat1': [0, 1, 2, 3], 'cat2': [4, 5, 6, 7]}
        row_gen = group_significance_row_generator(bt_4, sample_indices)
        exp_test_stats = [0.43577690622483684, -2.5911938781738648,
                          -
                          1.3573515147239095, 1.2101173913086851, 2.137178815882979,
                          0.0099191576638653078]
        exp_pvals = [0.67823972846362579, 0.041145883121579255,
                     0.2235024418313547, 0.27174025956151748, 0.076447615888438444,
                     0.9924073718332862]
        exp_means = [[52.25, 43.0],
                     [20.5, 44.0],
                     [29.25, 52.25],
                     [59.75, 44.5],
                     [39.0, 14.5],
                     [48.0, 47.75]]
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'parametric_t_test',
                                        GROUP_TEST_CHOICES)
        assert_almost_equal(exp_test_stats, obs_test_stats)
        assert_almost_equal(exp_pvals, obs_pvals)
        assert_almost_equal(exp_means, obs_means)

        # test with bootstrapped mann_whitney_u
        sample_indices = {'cat1': [4, 1, 2], 'cat2': [5, 0, 3]}
        row_gen = group_significance_row_generator(bt, sample_indices)
        exp_test_stats = [2.0, 2.0, 2.0, 3.0, 2.0, 2.0]
        exp_pvals = array([0.22977023, 0.2047952, 0.19280719, 0.44255744,
                           0.18781219, 0.22477522])
        exp_means = [[39.666666666666664, 61.0],
                     [24.333333333333332, 40.0],
                     [27.0, 42.333333333333336],
                     [57.0, 54.333333333333336],
                     [38.333333333333336, 19.666666666666668],
                     [30.333333333333332, 60.0]]
        seed(0)  # seed prng for reproducibility
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'bootstrap_mann_whitney_u',
                                        GROUP_TEST_CHOICES, reps=1000)
        assert_almost_equal(exp_test_stats, obs_test_stats)
        assert_almost_equal(exp_pvals, obs_pvals)
        assert_almost_equal(exp_means, obs_means)
        # test with BT_4
        sample_indices = {'cat1': [0, 1, 2, 3], 'cat2': [4, 5, 6, 7]}
        row_gen = group_significance_row_generator(bt_4, sample_indices)
        exp_test_stats = [6.0, 1.0, 5.0, 2.0, 1.0, 7.0]
        exp_pvals = array([0.49050949, 0.01598402, 0.31968032, 0.06193806,
                           0.02997003, 0.6953047])
        exp_means = [[52.25, 43.0],
                     [20.5, 44.0],
                     [29.25, 52.25],
                     [59.75, 44.5],
                     [39.0, 14.5],
                     [48.0, 47.75]]
        seed(0)  # seed prng for reproducibility
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'bootstrap_mann_whitney_u',
                                        GROUP_TEST_CHOICES, reps=1000)
        assert_almost_equal(exp_test_stats, obs_test_stats)
        assert_almost_equal(exp_pvals, obs_pvals)
        assert_almost_equal(exp_means, obs_means)

        # test with parametric mann whitney u
        sample_indices = {'cat1': [0, 3, 1], 'cat2': [4, 2, 5]}
        row_gen = group_significance_row_generator(bt, sample_indices)
        exp_test_stats = [3.0, 3.0, 4.0, 4.0, 3.0, 4.0]
        exp_pvals = array([0.66252058, 0.66252058, 1., 1., 0.66252058, 1.])
        exp_means = [[52.666666666666664, 48.0],
                     [23.666666666666668, 40.666666666666664],
                     [34.0, 35.333333333333336],
                     [56.333333333333336, 55.0],
                     [32.333333333333336, 25.666666666666668],
                     [46.0, 44.333333333333336]]
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'mann_whitney_u',
                                        GROUP_TEST_CHOICES)
        assert_almost_equal(exp_test_stats, obs_test_stats)
        assert_almost_equal(exp_pvals, obs_pvals)
        assert_almost_equal(exp_means, obs_means)
        # test with BT_4
        sample_indices = {'cat1': [0, 1, 2, 3], 'cat2': [4, 5, 6, 7]}
        row_gen = group_significance_row_generator(bt_4, sample_indices)
        exp_test_stats = [6.0, 1.0, 5.0, 2.0, 1.0, 7.0]
        exp_pvals = array([0.66500554, 0.06060197, 0.46782508, 0.1123512,
                           0.06060197, 0.88523391])
        exp_means = [[52.25, 43.0],
                     [20.5, 44.0],
                     [29.25, 52.25],
                     [59.75, 44.5],
                     [39.0, 14.5],
                     [48.0, 47.75]]
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'mann_whitney_u',
                                        GROUP_TEST_CHOICES)
        assert_almost_equal(exp_test_stats, obs_test_stats)
        assert_almost_equal(exp_pvals, obs_pvals)
        assert_almost_equal(exp_means, obs_means)

        # test with ANOVA
        sample_indices = {'cat1': [0, 3], 'cat2': [4, 5], 'cat3': [2, 1]}
        row_gen = group_significance_row_generator(bt, sample_indices)
        exp_test_stats = [0.022340083574413375, 20.028268551236753,
                          2.086854460093897, 0.96500593119810185, 4.8390804597701154,
                          0.54346882684796749]
        exp_pvals = [0.97806870848824634, 0.018391757629969238,
                     0.27043709109167957, 0.47468983920325486, 0.11510587547067222,
                     0.62890473306440042]
        exp_means = [[53.0, 46.5, 51.5],
                     [28.5, 55.5, 12.5],
                     [50.0, 45.5, 8.5],
                     [50.5, 47.5, 69.0],
                     [28.0, 9.0, 50.0],
                     [65.0, 39.5, 31.0]]
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'ANOVA',
                                        GROUP_TEST_CHOICES)
        assert_almost_equal(exp_test_stats, obs_test_stats)
        assert_almost_equal(exp_pvals, obs_pvals)
        assert_almost_equal(exp_means, obs_means)
        # test with uneven group sizes
        sample_indices = {'cat1': [0, 2, 3, 1], 'cat2': [4, 5]}
        row_gen = group_significance_row_generator(bt, sample_indices)
        exp_test_stats = [0.05663963168179019, 16.436058700209646,
                          0.43828937472444823, 0.675244322576109, 4.7713717693836974,
                          0.083541102077687446]
        exp_pvals = [0.8235822412182755, 0.015422975290359022,
                     0.54414414026513325, 0.45738578176242134, 0.094285405564661875,
                     0.78691584834507211]
        exp_means = [[52.25, 46.5],
                     [20.5, 55.5],
                     [29.25, 45.5],
                     [59.75, 47.5],
                     [39.0, 9.0],
                     [48.0, 39.5]]
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'ANOVA',
                                        GROUP_TEST_CHOICES)
        assert_almost_equal(exp_test_stats, obs_test_stats)
        assert_almost_equal(exp_pvals, obs_pvals)
        assert_almost_equal(exp_means, obs_means)
        # test with bt_4
        sample_indices = {'cat1': [0, 1, 2, 3], 'cat2': [4, 5, 6, 7]}
        row_gen = group_significance_row_generator(bt_4, sample_indices)
        exp_test_stats = [0.18990151199889027, 6.7142857142857144,
                          1.8424031345232912, 1.4643841007477372, 4.5675332910589734,
                          9.8389688760617899e-05]
        exp_pvals = [0.6782397284636259, 0.041145883121579234,
                     0.22350244183135481, 0.27174025956151771, 0.076447615888438403,
                     0.9924073718332751]
        exp_means = [[52.25, 43.0],
                     [20.5, 44.0],
                     [29.25, 52.25],
                     [59.75, 44.5],
                     [39.0, 14.5],
                     [48.0, 47.75]]
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'ANOVA',
                                        GROUP_TEST_CHOICES)
        assert_almost_equal(exp_test_stats, obs_test_stats)
        assert_almost_equal(exp_pvals, obs_pvals)
        assert_almost_equal(exp_means, obs_means)

        # test with g goodness of fit
        sample_indices = {'cat1': [0, 3], 'cat2': [4, 5], 'cat3': [2, 1]}
        row_gen = group_significance_row_generator(bt, sample_indices)
        exp_test_stats = [0.46328913071721711,
                          29.810689447160001, 37.234612591840595, 4.7031232724401875,
                          31.207185565457102, 13.332324853339509]
        exp_pvals = [0.79322801392154108,
                     3.3627225458535774e-07, 8.2149818410655555e-09, 0.09522034650579822,
                     1.6728066897036456e-07, 0.00127327567601971]
        exp_means = [[53.0, 46.5, 51.5],
                     [28.5, 55.5, 12.5],
                     [50.0, 45.5, 8.5],
                     [50.5, 47.5, 69.0],
                     [28.0, 9.0, 50.0],
                     [65.0, 39.5, 31.0]]
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'g_test',
                                        GROUP_TEST_CHOICES)
        assert_almost_equal(exp_test_stats, obs_test_stats)
        assert_almost_equal(exp_pvals, obs_pvals)
        assert_almost_equal(exp_means, obs_means)
        # test with individual groups
        sample_indices = {'cat1': [0], 'cat2': [1], 'cat3': [3],
                          'cat4': [2], 'cat5': [5], 'cat6': [4]}
        row_gen = group_significance_row_generator(bt, sample_indices)
        exp_test_stats = [68.7536611489639, 62.908926545455522,
                          115.84654226008865, 26.819713749563704, 84.940231595557307,
                          105.37909384565077]
        exp_pvals = [1.8616725644907271e-13, 3.0403858229558975e-12,
                     2.3772983815049693e-23, 6.1843461955812955e-05,
                     7.7481603433718027e-17, 3.8768150325829967e-21]
        exp_means = [[28.0, 52.0, 78.0, 51.0, 77.0, 16.0],
                     [25.0, 14.0, 32.0, 11.0, 63.0, 48.0],
                     [31.0, 2.0, 69.0, 15.0, 27.0, 64.0],
                     [36.0, 68.0, 65.0, 70.0, 62.0, 33.0],
                     [16.0, 41.0, 40.0, 59.0, 3.0, 15.0],
                     [32.0, 8.0, 98.0, 54.0, 50.0, 29.0]]
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'g_test',
                                        GROUP_TEST_CHOICES)
        assert_almost_equal(exp_test_stats, obs_test_stats)
        assert_almost_equal(exp_pvals, obs_pvals)
        assert_almost_equal(exp_means, obs_means)
        # test with uneven length groups
        sample_indices = {'cat1': [0, 3, 4, 5], 'cat3': [2, 1]}
        row_gen = group_significance_row_generator(bt, sample_indices)
        exp_test_stats = [0.030099778845406742,
                          16.703388149486191, 29.941854048163027, 3.39187772427496,
                          14.935738277477988, 5.4519230964604013]
        exp_pvals = [0.86226402523867973, 4.3702877865113464e-05,
                     4.451983032513133e-08, 0.065518295867083964, 0.00011123571448583719,
                     0.019546798231055287]
        exp_means = [[49.75, 51.5],
                     [42.0, 12.5],
                     [47.75, 8.5],
                     [49.0, 69.0],
                     [18.5, 50.0],
                     [52.25, 31.0]]
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'g_test',
                                        GROUP_TEST_CHOICES)
        assert_almost_equal(exp_test_stats, obs_test_stats)
        assert_almost_equal(exp_pvals, obs_pvals)
        assert_almost_equal(exp_means, obs_means)
        # test with bt_4
        sample_indices = {'cat1': [0, 1, 2, 3], 'cat2': [4, 5, 6, 7]}
        row_gen = group_significance_row_generator(bt_4, sample_indices)
        exp_test_stats = [0.8950130401309585, 8.6948783805472942,
                          6.5397009199496443, 2.2281537448054953, 11.541070115516771,
                          0.00064935138712822981]
        exp_pvals = [0.34412242732851783, 0.0031910540870178925,
                     0.010549308294222293, 0.13551569348660794, 0.00068075444949030543,
                     0.97967020739471489]
        exp_means = [[52.25, 43.0],
                     [20.5, 44.0],
                     [29.25, 52.25],
                     [59.75, 44.5],
                     [39.0, 14.5],
                     [48.0, 47.75]]
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'g_test',
                                        GROUP_TEST_CHOICES)
        assert_almost_equal(exp_test_stats, obs_test_stats)
        assert_almost_equal(exp_pvals, obs_pvals)
        assert_almost_equal(exp_means, obs_means)

        # test with Kruskal Wallis
        sample_indices = {'cat1': [0, 3], 'cat2': [4, 5], 'cat3': [2, 1]}
        row_gen = group_significance_row_generator(bt, sample_indices)
        exp_test_stats = [0.2857142857142847,
                          4.5714285714285694, 3.7142857142857117, 3.7142857142857117,
                          4.5714285714285694, 0.85714285714285765]
        exp_pvals = [0.86687789975018215, 0.10170139230422694,
                     0.15611804531597129, 0.15611804531597129, 0.10170139230422694,
                     0.65143905753105535]
        exp_means = [[53.0, 46.5, 51.5],
                     [28.5, 55.5, 12.5],
                     [50.0, 45.5, 8.5],
                     [50.5, 47.5, 69.0],
                     [28.0, 9.0, 50.0],
                     [65.0, 39.5, 31.0]]
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'kruskal_wallis',
                                        GROUP_TEST_CHOICES)
        assert_almost_equal(exp_test_stats, obs_test_stats)
        assert_almost_equal(exp_pvals, obs_pvals)
        assert_almost_equal(exp_means, obs_means)
        # test with individual groups
        sample_indices = {'cat1': [0], 'cat2': [1], 'cat3': [3],
                          'cat4': [2], 'cat5': [5], 'cat6': [4]}
        row_gen = group_significance_row_generator(bt, sample_indices)
        exp_test_stats = [5.0, 5.0, 5.0, 5.0, 5.0, 5.0]
        exp_pvals = [0.41588018699550794, 0.41588018699550794,
                     0.41588018699550794, 0.41588018699550794, 0.41588018699550794,
                     0.41588018699550794]
        exp_means = [[28.0, 52.0, 78.0, 51.0, 77.0, 16.0],
                     [25.0, 14.0, 32.0, 11.0, 63.0, 48.0],
                     [31.0, 2.0, 69.0, 15.0, 27.0, 64.0],
                     [36.0, 68.0, 65.0, 70.0, 62.0, 33.0],
                     [16.0, 41.0, 40.0, 59.0, 3.0, 15.0],
                     [32.0, 8.0, 98.0, 54.0, 50.0, 29.0]]
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'kruskal_wallis',
                                        GROUP_TEST_CHOICES)
        assert_almost_equal(exp_test_stats, obs_test_stats)
        assert_almost_equal(exp_pvals, obs_pvals)
        assert_almost_equal(exp_means, obs_means)
        # test with uneven length groups
        sample_indices = {'cat1': [0, 3, 4, 5], 'cat3': [2, 1]}
        row_gen = group_significance_row_generator(bt, sample_indices)
        exp_test_stats = [0.0, 3.428571428571427, 3.428571428571427,
                          3.428571428571427, 3.428571428571427, 0.21428571428571175]
        exp_pvals = [1, 0.064077506451059238, 0.064077506451059238,
                     0.064077506451059238, 0.064077506451059238, 0.64342884356362262]
        exp_means = [[49.75, 51.5],
                     [42.0, 12.5],
                     [47.75, 8.5],
                     [49.0, 69.0],
                     [18.5, 50.0],
                     [52.25, 31.0]]
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'kruskal_wallis',
                                        GROUP_TEST_CHOICES)
        assert_almost_equal(exp_test_stats, obs_test_stats)
        assert_almost_equal(exp_pvals, obs_pvals)
        assert_almost_equal(exp_means, obs_means)
        # test with bt_4
        sample_indices = {'cat1': [0, 1, 2, 3], 'cat2': [4, 5, 6, 7]}
        row_gen = group_significance_row_generator(bt_4, sample_indices)
        exp_test_stats = [0.33333333333333215, 4.0833333333333321,
                          0.75903614457831325, 3.0, 4.0833333333333321, 0.083333333333332149]
        exp_pvals = [0.56370286165077377, 0.043308142810792101,
                     0.38363032713198986, 0.08326451666355042, 0.043308142810792101,
                     0.77282999268444919]
        exp_means = [[52.25, 43.0],
                     [20.5, 44.0],
                     [29.25, 52.25],
                     [59.75, 44.5],
                     [39.0, 14.5],
                     [48.0, 47.75]]
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'kruskal_wallis',
                                        GROUP_TEST_CHOICES)
        assert_almost_equal(exp_test_stats, obs_test_stats)
        assert_almost_equal(exp_pvals, obs_pvals)
        assert_almost_equal(exp_means, obs_means)

    def test_group_significance_output_formatter(self):
        """output_formatter works"""
        # Using ANOVA test for example
        otu_table_1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "BIOM-Format 1.1.2","date": "2013-08-16T15:23:02.872397","matrix_type": "sparse","matrix_element_type": "float","shape": [6, 8],"data": [[0,0,28.0],[0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[0,6,73.0],[0,7,6.0],[1,0,25.0],[1,1,14.0],[1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[1,6,27.0],[1,7,38.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],[2,3,69.0],[2,4,64.0],[2,5,27.0],[2,6,64.0],[2,7,54.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],[3,4,33.0],[3,5,62.0],[3,6,60.0],[3,7,23.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],[4,5,3.0],[4,6,35.0],[4,7,5.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0],[5,6,93.0],[5,7,19.0]],"rows": [{"id": "OTU1", "metadata": {"taxonomy": "k__One"}},{"id": "OTU2", "metadata": {"taxonomy": "k__Two"}},{"id": "OTU3", "metadata": {"taxonomy": "k__Three"}},{"id": "OTU4", "metadata": {"taxonomy": "k__Four"}},{"id": "OTU5", "metadata": {"taxonomy": "k__Five"}},{"id": "OTU6", "metadata": {"taxonomy": "k__Six"}}],"columns": [{"id": "Sample1", "metadata": null},{"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},{"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},{"id": "Sample6", "metadata": null},{"id": "Sample7", "metadata": null},{"id": "Sample8", "metadata": null}]}"""
        bt = parse_biom_table(otu_table_1)

        test_stats = [0.18990151199889027,
                      6.7142857142857144,
                      1.8424031345232912,
                      1.4643841007477372,
                      4.5675332910589734,
                      9.8389688760617899e-05]
        pvals = [0.6782397284636259,
                 0.041145883121579234,
                 0.22350244183135481,
                 0.27174025956151771,
                 0.076447615888438403,
                 0.9924073718332751]
        means = [[52.25, 43.0],
                 [20.5, 44.0],
                 [29.25, 52.25],
                 [59.75, 44.5],
                 [39.0, 14.5],
                 [48.0, 47.75]]
        fdr_pvals = [0.8138876741563511,
                     0.24687529872947539,
                     0.44700488366270963,
                     0.40761038934227656,
                     0.22934284766531521,
                     0.9924073718332751]
        bon_pvals = [4.069438370781755,
                     0.2468752987294754,
                     1.3410146509881289,
                     1.6304415573691062,
                     0.4586856953306304,
                     5.95444423099965]
        cat_sample_indices = {'cat1': [0, 1, 2, 3], 'cat2': [4, 5, 6, 7]}

        line_header_out = 'OTU\tTest-Statistic\tP\tFDR_P\tBonferroni_P\tcat1_mean\tcat2_mean\ttaxonomy'
        line_1_out = 'OTU1\t0.189901511999\t0.678239728464\t0.813887674156\t4.06943837078\t52.25\t43.0\tk__One'
        line_6_out = 'OTU6\t9.83896887606e-05\t0.992407371833\t0.992407371833\t5.954444231\t48.0\t47.75\tk__Six'

        lines = group_significance_output_formatter(bt, test_stats, pvals,
                                                    fdr_pvals, bon_pvals, means, cat_sample_indices, md_key='taxonomy')

        self.assertEqual(lines[0], line_header_out)
        self.assertEqual(lines[1], line_1_out)
        self.assertEqual(lines[6], line_6_out)

    def test_sort_by_pval(self):
        """sort_by_pval works"""
        lines = [
            'OTU\tTest-Statistic\tP\tFDR_P\tBonferroni_P\tcat1_mean\tcat2_mean\tTaxonomy',
            'OTU1\t0.189901511999\t0.678239728464\t0.813887674156\t4.06943837078\t52.25\t43.0\tk__One',
            'OTU2\t6.71428571429\t0.0411458831216\t0.246875298729\t0.246875298729\t20.5\t44.0\tk__Two',
            'OTU3\t1.84240313452\t0.223502441831\t0.447004883663\t1.34101465099\t29.25\t52.25\tk__Three',
            'OTU4\t1.46438410075\t0.271740259562\t0.407610389342\t1.63044155737\t59.75\t44.5\tk__Four',
            'OTU5\t4.56753329106\t0.0764476158884\t0.229342847665\t0.458685695331\t39.0\t14.5\tk__Five',
            'OTU6\t9.83896887606e-05\t0.992407371833\t0.992407371833\t5.954444231\t48.0\t47.75\tk__Six']

        lines_sorted_pval_1 = \
            'OTU2\t6.71428571429\t0.0411458831216\t0.246875298729\t0.246875298729\t20.5\t44.0\tk__Two'
        lines_sorted_fdr_1 = \
            'OTU5\t4.56753329106\t0.0764476158884\t0.229342847665\t0.458685695331\t39.0\t14.5\tk__Five'
        lines_sorted_bonf_6 = \
            'OTU6\t9.83896887606e-05\t0.992407371833\t0.992407371833\t5.954444231\t48.0\t47.75\tk__Six'

        lines_pval = sort_by_pval(lines, 2)
        lines_pval_fdr = sort_by_pval(lines, 3)
        lines_pval_bonf = sort_by_pval(lines, 4)

        self.assertEqual(lines_pval[1], lines_sorted_pval_1)
        self.assertEqual(lines_pval_fdr[1], lines_sorted_fdr_1)
        self.assertEqual(lines_pval_bonf[6], lines_sorted_bonf_6)


class GroupedCorrelationTests(TestCase):

    """Tests of grouped correlation significance functions."""

    def setUp(self):
        """Define values used by all tests."""
        self.bt_str = '{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "BIOM-Format 1.2.0-dev","date": "2013-11-28T16:50:27.438635","matrix_type": "sparse","matrix_element_type": "float","shape": [10, 10],"data": [[0,0,22.0],[0,4,15.0],[0,5,74.0],[0,6,34.0],[0,7,76.0],[0,9,48.0],[1,0,70.0],[1,2,30.0],[1,3,37.0],[1,4,24.0],[1,5,77.0],[1,6,71.0],[1,7,58.0],[1,8,43.0],[2,1,2.0],[2,2,90.0],[2,4,48.0],[2,5,54.0],[2,7,22.0],[2,8,91.0],[3,0,80.0],[3,1,86.0],[3,3,78.0],[3,7,12.0],[4,4,68.0],[4,7,76.0],[4,8,57.0],[5,2,23.0],[5,5,66.0],[5,7,51.0],[5,9,77.0],[6,0,31.0],[6,1,47.0],[6,2,16.0],[6,4,96.0],[6,5,9.0],[7,0,17.0],[7,2,52.0],[7,5,11.0],[7,7,22.0],[8,1,74.0],[8,2,7.0],[8,4,80.0],[8,7,59.0],[9,0,6.0],[9,2,34.0],[9,3,63.0],[9,4,77.0],[9,5,8.0],[9,6,38.0],[9,7,73.0],[9,8,98.0],[9,9,45.0]],"rows": [{"id": "o1 ", "metadata": {"taxonomy": ["bug1"]}},{"id": "o2", "metadata": {"taxonomy": ["bug2"]}},{"id": "o3", "metadata": {"taxonomy": ["bug3"]}},{"id": "o4", "metadata": {"taxonomy": ["bug4"]}},{"id": "o5", "metadata": {"taxonomy": ["bug5"]}},{"id": "o6", "metadata": {"taxonomy": ["bug6"]}},{"id": "o7", "metadata": {"taxonomy": ["bug7"]}},{"id": "o8", "metadata": {"taxonomy": ["bug8"]}},{"id": "o9", "metadata": {"taxonomy": ["bug9"]}},{"id": "o10", "metadata": {"taxonomy": ["bug10"]}}],"columns": [{"id": "s1", "metadata": null},{"id": "s2", "metadata": null},{"id": "s6", "metadata": null},{"id": "s4", "metadata": null},{"id": "s5", "metadata": null},{"id": "s10", "metadata": null},{"id": "s7", "metadata": null},{"id": "s8", "metadata": null},{"id": "s9", "metadata": null},{"id": "s3", "metadata": null}]}'
        self.mf_ordered = ['#SampleIDt\thsid\tfield\tval',
                           's1\t1\tf1\t6.1',
                           's3\t1\tf1\t0.0',
                           's7\t1\tf1\t14.2',
                           's9\t1\tf2\t6.5',
                           's2\t1\tf2\t21',
                           's6\t2\tf3\t0.3',
                           's5\t2\tf2\t9.1',
                           's4\t2\tf3\t0.8',
                           's8\t2\tf2\t5.0',
                           's10\t2\tf2\t11.']
        self.mf_non_ordered = ['#SampleIDt\thsid\tval',
                               'Sample9\t1\t6.5',
                               'Sample8\t2\t5.0',
                               'Sample1\t1\t6.1',
                               'Sample2\t1\t21',
                               'Sample6\t2\t0.3',
                               'Sample5\t2\t9.1',
                               'Sample7\t1\t14.2',
                               'Sample4\t2\t0.8',
                               'Sample10\t2\t11.',
                               'Sample3\t1\t0.0']
        self.cvs1 = ['1', '2']
        self.mds1 = [[6.1, 0.0, 14.2, 6.5, 21], [.3, 9.1, .8, 5.0, 11.]]
        self.otus1 = [
            array([[22., 48., 34., 0., 0.],
                   [70., 0., 71., 43., 0.],
                   [0., 0., 0., 91., 2.],
                   [80., 0., 0., 0., 86.],
                   [0., 0., 0., 57., 0.],
                   [0., 77., 0., 0., 0.],
                   [31., 0., 0., 0., 47.],
                   [17., 0., 0., 0., 0.],
                   [0., 0., 0., 0., 74.],
                   [6., 45., 38., 98., 0.]]),
            array([[0., 15., 0., 76., 74.],
                   [30., 24., 37., 58., 77.],
                   [90., 48., 0., 22., 54.],
                   [0., 0., 78., 12., 0.],
                   [0., 68., 0., 76., 0.],
                   [23., 0., 0., 51., 66.],
                   [16., 96., 0., 0., 9.],
                   [52., 0., 0., 22., 11.],
                   [7., 80., 0., 59., 0.],
                   [34., 77., 63., 73., 8.]])]

    def test_grouped_correlation_row_generator(self):
        """Test that group row generator behaves as expected."""
        category = 'val'
        gc_to_samples = \
            {'1': ['s1', 's3', 's7', 's9', 's2'],
             '2': ['s6', 's5', 's4', 's8', 's10']}
        bt = parse_biom_table(self.bt_str)
        pmf = parse_mapping_file_to_dict(self.mf_ordered)[0]
        obs_cvs, obs_mds, obs_otus = grouped_correlation_row_generator(bt, pmf,
                                                                       category, gc_to_samples)
        self.assertEqual(obs_cvs, self.cvs1)
        assert_almost_equal(obs_mds, self.mds1)
        assert_almost_equal(obs_otus, self.otus1)
        # make sure it throws an error on non float md
        self.assertRaises(ValueError, grouped_correlation_row_generator, bt,
                          pmf, 'field', gc_to_samples)

    def test_run_grouped_correlation(self):
        """Test that grouped correlation values are calculated as expected."""
        # hand calculation of spearman and pearson for 01
        # md_g1 = array([6.1, 0.0, 14.2, 6.5, 21])
        # md_g2 = array([.3, 9.1, .8, 5.0, 11])
        # o1_g1 = array([22, 48, 34, 0, 0])
        # o1_g2 = array([0, 15, 0, 76, 74])
        # c1_g1 = -0.6155870112510925 #spearman(md_g1, o1_g1)
        # c2_g2 = 0.66688592885535025 #spearman(md_g2, o1_g2)
        # fisher_population_correlation([-0.6155870112510925,
        # 0.66688592885535025], [5,5])
        # fpc, h = (0.043595171909468329, 0.12776325359984511)
        g1_rhos = [corrcoef(self.otus1[0][i], self.mds1[0])[0][1]
                   for i in range(10)]
        g2_rhos = [corrcoef(self.otus1[1][i], self.mds1[1])[0][1]
                   for i in range(10)]
        exp_rhos = [g1_rhos, g2_rhos]
        g1_pvals = [assign_correlation_pval(g1_rhos[i], 5,
                                            'parametric_t_distribution') for i in range(10)]
        g2_pvals = [assign_correlation_pval(g2_rhos[i], 5,
                                            'parametric_t_distribution') for i in range(10)]
        exp_pvals = [g1_pvals, g2_pvals]
        exp_f_pvals = [fisher([g1_pvals[i], g2_pvals[i]]) for i in range(10)]

        tmp = [fisher_population_correlation([g1_rhos[i], g2_rhos[i]], [5, 5])
               for i in range(10)]
        exp_f_rhos = [x[0] for x in tmp]
        exp_f_hs = [x[1] for x in tmp]

        obs_rhos, obs_pvals, obs_f_pvals, obs_f_rhos, obs_f_hs = \
            run_grouped_correlation(self.mds1, self.otus1, 'pearson',
                                    CORRELATION_TEST_CHOICES, 'parametric_t_distribution')

        assert_almost_equal(obs_rhos, exp_rhos)
        assert_almost_equal(obs_pvals, exp_pvals)
        assert_almost_equal(obs_f_pvals, exp_f_pvals)
        assert_almost_equal(obs_f_rhos, exp_f_rhos)
        assert_almost_equal(obs_f_hs, exp_f_hs)

    def test_grouped_correlation_formatter(self):
        """Test that grouped correlation results are correctly formatted."""
        obs_rhos, obs_pvals, obs_f_pvals, obs_f_rhos, obs_f_hs = \
            run_grouped_correlation(self.mds1, self.otus1, 'pearson',
                                    CORRELATION_TEST_CHOICES, 'parametric_t_distribution')
        obs_lines = grouped_correlation_formatter(
            parse_biom_table(self.bt_str),
            obs_rhos, obs_pvals, obs_f_rhos, obs_f_pvals, obs_f_hs,
            grouping_category='val', category_values=['1', '2'],
            md_key='taxonomy')
        exp_lines = [
            'OTU\tRho_val:1\tRho_val:2\tPval_val:1\tPval_val:2\tFisher population correlation\tFisher combined p\tHomogeneity pval\ttaxonomy',
            'o1 \t-0.549008525652\t0.624559974829\t0.337884839311\t0.26003464759\t0.0576788156554\t0.301540747632\t0.177206041905\tbug1',
            'o2\t-0.0384383440919\t0.498044011456\t0.951070834767\t0.393160408737\t0.248789595012\t0.741753460366\t0.558440916618\tbug2',
            'o3\t-0.193866521776\t0.0709448388675\t0.754716520102\t0.909746057157\t-0.0625618635573\t0.944764075331\t0.789149030542\tbug3',
            'o4\t0.477058782953\t-0.535587095672\t0.416488072837\t0.352242311475\t-0.039368393563\t0.428279994782\t0.263944507542\tbug4',
            'o5\t-0.210109340238\t0.321581922552\t0.734462039847\t0.5977201336\t0.0599902605341\t0.80041158791\t0.584587444442\tbug5',
            'o6\t-0.656420030285\t0.443371522913\t0.228885504776\t0.454564967068\t-0.153808372937\t0.339487962449\t0.206619183728\tbug6',
            'o7\t0.598119568608\t0.439213235895\t0.286674057126\t0.45931604814\t0.523199343913\t0.398633295268\t0.826682300837\tbug7',
            'o8\t-0.237574613471\t-0.479408280331\t0.700380699877\t0.413860864444\t-0.364624217025\t0.648810427511\t0.779467050573\tbug8',
            'o9\t0.785506814483\t0.359680356852\t0.115335393753\t0.552116179836\t0.615702620054\t0.239043351127\t0.494561919352\tbug9',
            'o10\t-0.451199640732\t-0.216326825605\t0.445650543742\t0.72672774463\t-0.339035584667\t0.689001456839\t0.789926403391\tbug10']

        self.assertEqual(obs_lines, exp_lines)
        self.assertEqual(obs_lines, exp_lines)

    def test_correlation_row_generator(self):
        """Test that correlation_row_generator works"""
        otu_table_1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0",
        "format_url": "http://biom-format.org","type": "OTU table","generated_by":
        "BIOM-Format 1.1.2","date": "2013-08-16T10:16:20.131837","matrix_type":
        "sparse","matrix_element_type": "float","shape": [6, 6],"data": [[0,0,28.0],
        [0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[1,0,25.0],[1,1,14.0],
        [1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],
        [2,3,69.0],[2,4,64.0],[2,5,27.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],
        [3,4,33.0],[3,5,62.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],
        [4,5,3.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0]],
        "rows": [{"id": "OTU1", "metadata": {"taxonomy": ["k__One"]}},{"id": "OTU2",
        "metadata": {"taxonomy": ["k__Two"]}},{"id": "OTU3", "metadata": {"taxonomy":
        ["k__Three"]}},{"id": "OTU4", "metadata": {"taxonomy": ["k__Four"]}},{"id":
        "OTU5", "metadata": {"taxonomy": ["k__Five"]}},{"id": "OTU6", "metadata":
        {"taxonomy": ["k__Six"]}}],"columns": [{"id": "Sample1", "metadata": null},
        {"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},
        {"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},
        {"id": "Sample6", "metadata": null}]}"""
        bt = parse_biom_table(otu_table_1)
        pmf = {'Sample1': {'test_cat': 'cat1', 'test_corr': '1'},
               'Sample2': {'test_cat': 'cat1', 'test_corr': '2'},
               'Sample3': {'test_cat': 'cat2', 'test_corr': '3'},
               'Sample4': {'test_cat': 'cat2', 'test_corr': '4'},
               'Sample5': {'test_cat': 'cat3', 'test_corr': '5'},
               'Sample6': {'test_cat': 'cat3', 'test_corr': '6'}}
        data_result = []
        corr_row_gen_data = correlation_row_generator(bt, pmf, 'test_corr')
        for i in corr_row_gen_data:
            data_result.append(i)

        data_output = [(array([28., 52., 51., 78., 16., 77.]),
                        array([1., 2., 3., 4., 5., 6.])),
                       (array([25., 14., 11., 32., 48., 63.]),
                        array([1., 2., 3., 4., 5., 6.])),
                       (array([31., 2., 15., 69., 64., 27.]),
                        array([1., 2., 3., 4., 5., 6.])),
                       (array([36., 68., 70., 65., 33., 62.]),
                        array([1., 2., 3., 4., 5., 6.])),
                       (array([16., 41., 59., 40., 15., 3.]),
                        array([1., 2., 3., 4., 5., 6.])),
                       (array([32., 8., 54., 98., 29., 50.]),
                        array([1., 2., 3., 4., 5., 6.]))]

        assert_almost_equal(data_result, data_output)

    def test_run_correlation_test(self):
        """Test run_correlation_test works."""
        otu_table_1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0",
        "format_url": "http://biom-format.org","type": "OTU table","generated_by":
        "BIOM-Format 1.1.2","date": "2013-08-16T10:16:20.131837","matrix_type":
        "sparse","matrix_element_type": "float","shape": [6, 6],"data": [[0,0,28.0],
        [0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[1,0,25.0],[1,1,14.0],
        [1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],
        [2,3,69.0],[2,4,64.0],[2,5,27.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],
        [3,4,33.0],[3,5,62.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],
        [4,5,3.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0]],
        "rows": [{"id": "OTU1", "metadata": {"taxonomy": ["k__One"]}},{"id": "OTU2",
        "metadata": {"taxonomy": ["k__Two"]}},{"id": "OTU3", "metadata": {"taxonomy":
        ["k__Three"]}},{"id": "OTU4", "metadata": {"taxonomy": ["k__Four"]}},{"id":
        "OTU5", "metadata": {"taxonomy": ["k__Five"]}},{"id": "OTU6", "metadata":
        {"taxonomy": ["k__Six"]}}],"columns": [{"id": "Sample1", "metadata": null},
        {"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},
        {"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},
        {"id": "Sample6", "metadata": null}]}"""
        bt = parse_biom_table(otu_table_1)
        pmf = {'Sample1': {'test_cat': 'cat1', 'test_corr': '1'},
               'Sample2': {'test_cat': 'cat1', 'test_corr': '2'},
               'Sample3': {'test_cat': 'cat2', 'test_corr': '3'},
               'Sample4': {'test_cat': 'cat2', 'test_corr': '4'},
               'Sample5': {'test_cat': 'cat3', 'test_corr': '5'},
               'Sample6': {'test_cat': 'cat3', 'test_corr': '6'}}

        data_gen = correlation_row_generator(bt, pmf, 'test_corr')

        # pearson
        exp_ccs = [0.34884669332532803,
                   0.83015306552396662,
                   0.44037594794853846,
                   0.064224958686639605,
                   -0.41225244540969846,
                   0.34313146667442163]
        exp_pvals = [0.49795623745457107,
                     0.04082210114420709,
                     0.38213734667034305,
                     0.90379502098011888,
                     0.41665291191825937,
                     0.50550281276602604]
        exp_bootstrapped_pvals = [0.48599999999999999,
                                  0.045999999999999999,
                                  0.38200000000000001,
                                  0.88200000000000001,
                                  0.438]

        obs_ccs, obs_pvals = run_correlation_test(data_gen, 'pearson',
                                                  CORRELATION_TEST_CHOICES,
                                                  pval_assignment_method='parametric_t_distribution')
        assert_almost_equal(exp_ccs, obs_ccs)
        assert_almost_equal(exp_pvals, obs_pvals)

        # test with bootstrapped pvalues
        seed(0)
        data_gen = correlation_row_generator(bt, pmf, 'test_corr')
        obs_ccs, obs_pvals = run_correlation_test(data_gen, 'pearson',
                                                  CORRELATION_TEST_CHOICES, pval_assignment_method='bootstrapped',
                                                  permutations=1000)

        # We can use the default epsilon in assertFloatEqual to test all
        # p-values but the last one. This last one needs a different epsilon
        # because the 479th permuted correlation coefficient can differ in
        # equality to the observed correlation coefficient depending on the
        # platform/numpy installation/configuration. Thus, the count of
        # more-extreme correlation coefficients can differ by 1. With 1000
        # permutations used in the test, this difference is 1 / 1000 = 0.001,
        # hence the larger epsilon used here.
        assert_almost_equal(obs_pvals[:-1], exp_bootstrapped_pvals)
        assert_almost_equal(obs_pvals[-1], 0.54600000000000004)

        # spearman
        exp_ccs = [0.25714285714285712,
                   0.77142857142857146,
                   0.31428571428571428,
                   -0.25714285714285712,
                   -0.60000000000000009,
                   0.25714285714285712]
        exp_pvals = [0.62278717201166178,
                     0.07239650145772597,
                     0.54409329446064136,
                     0.62278717201166178,
                     0.20799999999999996,
                     0.62278717201166178]
        exp_bootstrapped_pvals = [0.66200000000000003,
                                  0.105,
                                  0.56299999999999994,
                                  0.65100000000000002,
                                  0.26300000000000001,
                                  0.65500000000000003]

        data_gen = correlation_row_generator(bt, pmf, 'test_corr')
        obs_ccs, obs_pvals = run_correlation_test(data_gen, 'spearman',
                                                  CORRELATION_TEST_CHOICES,
                                                  pval_assignment_method='parametric_t_distribution')

        assert_almost_equal(exp_ccs, obs_ccs)
        assert_almost_equal(exp_pvals, obs_pvals)

        # test with bootstrapped pvalues
        seed(0)
        data_gen = correlation_row_generator(bt, pmf, 'test_corr')
        obs_ccs, obs_pvals = run_correlation_test(data_gen, 'spearman',
                                                  CORRELATION_TEST_CHOICES, pval_assignment_method='bootstrapped',
                                                  permutations=1000)
        assert_almost_equal(exp_bootstrapped_pvals, obs_pvals)

        # kendall
        data_gen = correlation_row_generator(bt, pmf, 'test_corr')

        exp_ccs = [0.2, 0.6, 0.2, -0.2, -0.4666666666666667, 0.2]
        exp_pvals = [0.57302511935539036,
                     0.090873939986249083,
                     0.57302511935539036,
                     0.57302511935539036,
                     0.18848603806737496,
                     0.57302511935539036]
        exp_bootstrapped_pvals = [0.73199999999999998,
                                  0.14000000000000001,
                                  0.72599999999999998,
                                  0.70499999999999996,
                                  0.29899999999999999,
                                  0.71399999999999997]

        data_gen = correlation_row_generator(bt, pmf, 'test_corr')
        obs_ccs, obs_pvals = run_correlation_test(data_gen, 'kendall',
                                                  CORRELATION_TEST_CHOICES,
                                                  pval_assignment_method='kendall')

        assert_almost_equal(exp_ccs, obs_ccs)
        assert_almost_equal(exp_pvals, obs_pvals)

        # test with bootstrapped pvalues
        seed(0)
        data_gen = correlation_row_generator(bt, pmf, 'test_corr')
        obs_ccs, obs_pvals = run_correlation_test(data_gen, 'kendall',
                                                  CORRELATION_TEST_CHOICES, pval_assignment_method='bootstrapped',
                                                  permutations=1000)
        assert_almost_equal(exp_bootstrapped_pvals, obs_pvals)

# globals used by certain tests.
BT_IN_1 = '{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "testCode","date": "2013-08-20T15:48:21.166180","matrix_type": "sparse","matrix_element_type": "float","shape": [6, 6],"data": [[0,0,28.0],[0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[1,0,25.0],[1,1,14.0],[1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],[2,3,69.0],[2,4,64.0],[2,5,27.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],[3,4,33.0],[3,5,62.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],[4,5,3.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0]],"rows": [{"id": "OTU1", "metadata": {"taxonomy": ["k__One"]}},{"id": "OTU2", "metadata": {"taxonomy": ["k__Two"]}},{"id": "OTU3", "metadata": {"taxonomy": ["k__Three"]}},{"id": "OTU4", "metadata": {"taxonomy": ["k__Four"]}},{"id": "OTU5", "metadata": {"taxonomy": ["k__Five"]}},{"id": "OTU6", "metadata": {"taxonomy": ["k__Six"]}}],"columns": [{"id": "Sample1", "metadata": null},{"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},{"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},{"id": "Sample6", "metadata": null}]}'
MF_IN_1 = ['#SampleID\ttest_cat\ttest_corr',
           'Sample1\tcat1\t1',
           'Sample2\tcat1\t2',
           'Sample3\tcat2\t3',
           'Sample4\tcat2\t4',
           'Sample5\tcat3\t5',
           'Sample6\tcat3\t6',
           'NotInOtuTable1\tcat5\t7',
           'NotInOtuTable2\tcat5\t8']
BT_4 = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "BIOM-Format 1.1.2","date": "2013-08-16T15:23:02.872397","matrix_type": "sparse","matrix_element_type": "float","shape": [6, 8],"data": [[0,0,28.0],[0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[0,6,73.0],[0,7,6.0],[1,0,25.0],[1,1,14.0],[1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[1,6,27.0],[1,7,38.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],[2,3,69.0],[2,4,64.0],[2,5,27.0],[2,6,64.0],[2,7,54.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],[3,4,33.0],[3,5,62.0],[3,6,60.0],[3,7,23.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],[4,5,3.0],[4,6,35.0],[4,7,5.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0],[5,6,93.0],[5,7,19.0]],"rows": [{"id": "OTU1", "metadata": {"taxonomy": "k__One"}},{"id": "OTU2", "metadata": {"taxonomy": "k__Two"}},{"id": "OTU3", "metadata": {"taxonomy": "k__Three"}},{"id": "OTU4", "metadata": {"taxonomy": "k__Four"}},{"id": "OTU5", "metadata": {"taxonomy": "k__Five"}},{"id": "OTU6", "metadata": {"taxonomy": "k__Six"}}],"columns": [{"id": "Sample1", "metadata": null},{"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},{"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},{"id": "Sample6", "metadata": null},{"id": "Sample7", "metadata": null},{"id": "Sample8", "metadata": null}]}"""

# run unit tests if run from command-line
if __name__ == '__main__':
    main()
