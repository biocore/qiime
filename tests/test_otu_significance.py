#!/usr/bin/env python
# File created on 16 Aug 2013
from __future__ import division

__author__ = "Luke Ursell"
__copyright__ = "Copyright 2013, The QIIME project"
__credits__ = ["Luke Ursell, Will Van Treuren"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Luke Ursell"
__email__ = "lkursell@gmail.com"
__status__ = "Development"

from cogent.util.unit_test import TestCase, main 
from qiime.otu_significance import (get_sample_cats, get_cat_sample_groups, 
    get_sample_indices, group_significance_row_generator, sort_by_pval,
    run_group_significance_test, group_significance_output_formatter,
    GROUP_TEST_CHOICES, longitudinal_row_generator, 
    run_longitudinal_correlation_test, CORRELATION_TEST_CHOICES)
from numpy import array, hstack
from numpy.random import seed
from cogent.util.dict2d import Dict2D
from qiime.util import get_tmp_filename
from os import remove
from qiime.parse import parse_mapping_file_to_dict, parse_otu_table
from biom.parse import parse_biom_table
from qiime.format import format_biom_table

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

class TestGroupSignificanceFunctions(TestCase):
    """Tests of group significance functions."""

    def setUp(self):
        """Define values used by all tests."""
        # nothing to do, all tests use different things
        pass

    def test_get_sample_cats(self):
        """Test get_sample_cats."""
        pmf_in = {'Sample1': {'test_cat': 'cat1', 'test_corr': '1', 'test_empty': 'abc'},
            'Sample2': {'test_cat': 'cat1', 'test_corr': '1', 'test_empty': 'abc'},
            'Sample3': {'test_cat': 'cat2', 'test_corr': '1', 'test_empty': ''},
            'Sample4': {'test_cat': 'cat2', 'test_corr': '1', 'test_empty': ''},
            'Sample5': {'test_cat': 'cat3', 'test_corr': '1', 'test_empty': 'abc'},
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
        # note that the order within a group doesn't 
        # matter for the group significance row generator because each set of 
        # values is whats important, not their order

    def test_get_sample_indices(self):
        """Test get_sample_indices works"""
        cat_sample_groups = {'cat1': ['Sample1', 'Sample2'],
            'cat2': ['Sample4', 'Sample3'],
            'cat3': ['Sample5', 'Sample6']}
        bt = parse_biom_table(BT_IN_1)
        # note that the order within a group doesn't 
        # matter for the group significance row generator because each set of 
        # values is whats important, not their order
        exp = {'cat1': [0, 1], 'cat2': [3, 2], 'cat3': [4, 5]}
        obs = get_sample_indices(cat_sample_groups, bt)
        self.assertEqual(exp, obs)

    def test_group_significance_row_generator(self):
        """Test group_significance_row_generator works."""
        # run with ordered example
        sample_indices = {'cat1': [0, 1], 'cat2': [3, 2], 'cat3': [4, 5]}
        bt = parse_biom_table(BT_IN_1)
        data = array([bt.observationData(i) for i in bt.ObservationIds])
        obs = list(group_significance_row_generator(bt, sample_indices))
        exp = zip(data.take([0,1],1), data.take([3,2],1), data.take([4,5],1))
        for o,e in zip(obs, exp):
            self.assertFloatEqual(e,o)
        # run with unequal length example
        sample_indices = {'g0': [0,1,4,5], 'g1': [3], 'g2': [2]}
        obs = list(group_significance_row_generator(bt, sample_indices))
        exp = zip(data.take([0,1,4,5],1), data.take([3],1), data.take([2],1))
        for o,e in zip(obs, exp):
            [self.assertFloatEqual(i,j) for i, j in zip(sorted(hstack(e)),sorted(hstack(o)))]
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
            -1.5065313062753816, -0.043884559904114794, -1.0631239617935129,
            -1.2878361428003895]
        exp_pvals = [0.887, 0.898, .278, 1.0, .488, .298]
        exp_means = [[52.333333333333336, 48.333333333333336],
            [34.0, 30.333333333333332],
            [20.0, 49.333333333333336],
            [55.333333333333336, 56.0],
            [20.0, 38.0],
            [30.0, 60.333333333333336]]
        seed(0) # seed prng for reproducibility
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'nonparametric_t_test', 
                GROUP_TEST_CHOICES, reps=1000)
        self.assertFloatEqual(exp_test_stats, obs_test_stats)
        self.assertFloatEqual(exp_pvals, obs_pvals)
        self.assertFloatEqual(exp_means, obs_means)
        # test with nonparametric t-test but different ordering 
        sample_indices = {'cat1': [0, 1, 5], 'cat2': [4, 3, 2]}
        row_gen = group_significance_row_generator(bt, sample_indices)
        seed(0) # seed prng for reproducibility
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'nonparametric_t_test', 
                GROUP_TEST_CHOICES, reps=1000)
        self.assertFloatEqual(exp_test_stats, obs_test_stats)
        self.assertFloatEqual(exp_pvals, obs_pvals)
        self.assertFloatEqual(exp_means, obs_means)
        # test with BT_4 biom table
        sample_indices = {'cat1': [0,3,1,4], 'cat2': [5,2,7,6]}
        row_gen = group_significance_row_generator(bt_4, sample_indices)
        exp_test_stats = [-0.38741397129147953, -0.38334158591463874,
            0.077468274988510541, -0.2322539745918096, 0.16469600468808282,
            -0.49589486133213057]
        exp_pvals = [.82,.718,.915,.934,.937,.603]
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
        self.assertFloatEqual(exp_test_stats, obs_test_stats)
        self.assertFloatEqual(exp_pvals, obs_pvals)
        self.assertFloatEqual(exp_means, obs_means)

        # test with parametric t test
        # bt_1 agrees with Prism
        sample_indices = {'cat1': [4, 1, 2], 'cat2': [5, 0, 3]}
        row_gen = group_significance_row_generator(bt, sample_indices)
        exp_test_stats = [-1.0504514628777806, -0.94113003446934629,
            -0.66264262463016887, 0.17617555832772411,  1.1144416530351877,
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
        self.assertFloatEqual(exp_test_stats, obs_test_stats)
        self.assertFloatEqual(exp_pvals, obs_pvals)
        self.assertFloatEqual(exp_means, obs_means)
        # test with BT_4
        sample_indices = {'cat1': [0, 1, 2, 3], 'cat2': [4, 5, 6, 7]}
        row_gen = group_significance_row_generator(bt_4, sample_indices)
        exp_test_stats = [0.43577690622483684, -2.5911938781738648,
            -1.3573515147239095, 1.2101173913086851, 2.137178815882979,
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
        self.assertFloatEqual(exp_test_stats, obs_test_stats)
        self.assertFloatEqual(exp_pvals, obs_pvals)
        self.assertFloatEqual(exp_means, obs_means)

        # test with bootstrapped mann_whitney_u
        sample_indices = {'cat1': [4, 1, 2], 'cat2': [5, 0, 3]}
        row_gen = group_significance_row_generator(bt, sample_indices)
        exp_test_stats = [7.0, 7.0, 7.0, 6.0, 7.0, 7.0]
        exp_pvals = [0.333, 0.305, 0.3, 0.623, 0.295, 0.334]
        exp_means = [[39.666666666666664, 61.0],
            [24.333333333333332, 40.0],
            [27.0, 42.333333333333336],
            [57.0, 54.333333333333336],
            [38.333333333333336, 19.666666666666668],
            [30.333333333333332, 60.0]]
        seed(0) # seed prng for reproducibility
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'bootstrap_mann_whitney_u', 
                GROUP_TEST_CHOICES, reps=1000)
        self.assertFloatEqual(exp_test_stats, obs_test_stats)
        self.assertFloatEqual(exp_pvals, obs_pvals)
        self.assertFloatEqual(exp_means, obs_means)
        # test with BT_4
        sample_indices = {'cat1': [0, 1, 2, 3], 'cat2': [4, 5, 6, 7]}
        row_gen = group_significance_row_generator(bt_4, sample_indices)  
        exp_test_stats = [10.0, 15.0, 11.0, 14.0, 15.0, 9.0]
        exp_pvals = [0.605, 0.033, 0.414, 0.097, 0.041, 0.814]
        exp_means = [[52.25, 43.0],
             [20.5, 44.0],
             [29.25, 52.25],
             [59.75, 44.5],
             [39.0, 14.5],
             [48.0, 47.75]]
        seed(0) # seed prng for reproducibility
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'bootstrap_mann_whitney_u', 
                GROUP_TEST_CHOICES, reps=1000)
        self.assertFloatEqual(exp_test_stats, obs_test_stats)
        self.assertFloatEqual(exp_pvals, obs_pvals)
        self.assertFloatEqual(exp_means, obs_means)

        # test with parametric mann whitney u
        sample_indices = {'cat1': [0, 3, 1], 'cat2': [4, 2, 5]}
        row_gen = group_significance_row_generator(bt, sample_indices)
        exp_test_stats = [6.0, 6.0, 5.0, 5.0, 6.0, 5.0]
        exp_pvals = [0.51269076026192328, 0.51269076026192328,
            0.82725934656271127, 0.82725934656271127, 0.51269076026192328,
            0.82725934656271127]
        exp_means = [[52.666666666666664, 48.0],
            [23.666666666666668, 40.666666666666664],
            [34.0, 35.333333333333336],
            [56.333333333333336, 55.0],
            [32.333333333333336, 25.666666666666668],
            [46.0, 44.333333333333336]]
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'mann_whitney_u',
                GROUP_TEST_CHOICES)
        self.assertFloatEqual(exp_test_stats, obs_test_stats)
        self.assertFloatEqual(exp_pvals, obs_pvals)
        self.assertFloatEqual(exp_means, obs_means)
        # test with BT_4
        sample_indices = {'cat1': [0, 1, 2, 3], 'cat2': [4, 5, 6, 7]}
        row_gen = group_significance_row_generator(bt_4, sample_indices)
        exp_test_stats = [10.0, 15.0, 11.0, 14.0, 15.0, 9.0]
        exp_pvals = [0.5637028616507731, 0.043308142810791955,
            0.38363032713198975, 0.083264516663550406, 0.043308142810791955,
            0.77282999268444752]
        exp_means = [[52.25, 43.0],
             [20.5, 44.0],
             [29.25, 52.25],
             [59.75, 44.5],
             [39.0, 14.5],
             [48.0, 47.75]]
        obs_test_stats, obs_pvals, obs_means = \
            run_group_significance_test(row_gen, 'mann_whitney_u',
                GROUP_TEST_CHOICES)
        self.assertFloatEqual(exp_test_stats, obs_test_stats)
        self.assertFloatEqual(exp_pvals, obs_pvals)
        self.assertFloatEqual(exp_means, obs_means)

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
        self.assertFloatEqual(exp_test_stats, obs_test_stats)
        self.assertFloatEqual(exp_pvals, obs_pvals)
        self.assertFloatEqual(exp_means, obs_means)
        # test with uneven group sizes
        sample_indices = {'cat1': [0, 2 ,3, 1], 'cat2': [4, 5]}
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
        self.assertFloatEqual(exp_test_stats, obs_test_stats)
        self.assertFloatEqual(exp_pvals, obs_pvals)
        self.assertFloatEqual(exp_means, obs_means)
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
        self.assertFloatEqual(exp_test_stats, obs_test_stats)
        self.assertFloatEqual(exp_pvals, obs_pvals)
        self.assertFloatEqual(exp_means, obs_means)

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
        self.assertFloatEqual(exp_test_stats, obs_test_stats)
        self.assertFloatEqual(exp_pvals, obs_pvals)
        self.assertFloatEqual(exp_means, obs_means)
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
        self.assertFloatEqual(exp_test_stats, obs_test_stats)
        self.assertFloatEqual(exp_pvals, obs_pvals)
        self.assertFloatEqual(exp_means, obs_means)
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
        self.assertFloatEqual(exp_test_stats, obs_test_stats)
        self.assertFloatEqual(exp_pvals, obs_pvals)
        self.assertFloatEqual(exp_means, obs_means)
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
        self.assertFloatEqual(exp_test_stats, obs_test_stats)
        self.assertFloatEqual(exp_pvals, obs_pvals)
        self.assertFloatEqual(exp_means, obs_means)

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
        self.assertFloatEqual(exp_test_stats, obs_test_stats)
        self.assertFloatEqual(exp_pvals, obs_pvals)
        self.assertFloatEqual(exp_means, obs_means)
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
        self.assertFloatEqual(exp_test_stats, obs_test_stats)
        self.assertFloatEqual(exp_pvals, obs_pvals)
        self.assertFloatEqual(exp_means, obs_means)
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
        self.assertFloatEqual(exp_test_stats, obs_test_stats)
        self.assertFloatEqual(exp_pvals, obs_pvals)
        self.assertFloatEqual(exp_means, obs_means)
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
        self.assertFloatEqual(exp_test_stats, obs_test_stats)
        self.assertFloatEqual(exp_pvals, obs_pvals)
        self.assertFloatEqual(exp_means, obs_means)

    def test_group_significance_output_formatter(self):
        """output_formatter works"""
        #Using ANOVA test for example
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

        line_header_out = 'OTU\tTest-Statistic\tP\tFDR_P\tBonferroni_P\tcat1_mean\tcat2_mean\tTaxonomy'
        line_1_out = 'OTU1\t0.189901511999\t0.678239728464\t0.813887674156\t4.06943837078\t52.25\t43.0\tk__One'
        line_6_out = 'OTU6\t9.83896887606e-05\t0.992407371833\t0.992407371833\t5.954444231\t48.0\t47.75\tk__Six'

        lines = group_significance_output_formatter(bt, test_stats, pvals, 
            fdr_pvals, bon_pvals, means, cat_sample_indices)

        self.assertEqual(lines[0], line_header_out)
        self.assertEqual(lines[1], line_1_out)
        self.assertEqual(lines[6], line_6_out)

    def test_sort_by_pval(self):
        """sort_by_pval works"""
        lines = ['OTU\tTest-Statistic\tP\tFDR_P\tBonferroni_P\tcat1_mean\tcat2_mean\tTaxonomy',
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

    # def test_longitudinal_row_generator(self):
    #     """Test that longitudinal row generator behaves as expected."""
    #     MF_2 = ['#SampleIDt\thsid\tval',
    #         'Sample1\t1\t1',
    #         'Sample2\t1\t2',
    #         'Sample3\t2\t3',
    #         'Sample4\t2\t4',
    #         'Sample5\t3\t5',
    #         'Sample6\t3\t6']
    #     bt = parse_biom_table(BT_IN_1)
    #     pmf, _ = parse_mapping_file_to_dict(MF_2)
    #     category = 'val'
    #     hsid_to_samples = {'1': ['Sample1', 'Sample2'],
    #         '2': ['Sample4', 'Sample3'],
    #         '3': ['Sample5', 'Sample6']}
    #     hsid_to_sample_indices = {'1': [0, 1], '2': [3, 2], '3': [4, 5]}
    #     # for convinience 
    #     data = array([bt.observationData(i) for i in bt.ObservationIds])
    #     exp_otus = [[i.take([0,1]), i.take([4,5]), i.take([3,2])] for i in data]
    #     exp_grad_vals = [array([1,2]), array([5,6]), array([4,3])]
    #     obs_data = longitudinal_row_generator(bt, pmf, category, 
    #         hsid_to_samples, hsid_to_sample_indices)
    #     for i,j in zip(obs_data, exp_otus):
    #         self.assertFloatEqual(i[0], j)
    #         self.assertFloatEqual(i[1], exp_grad_vals)
    #     # test with hand calculated example
    #     bt_str = '{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "testCode","date": "2013-08-22T18:23:28.199054","matrix_type": "sparse","matrix_element_type": "float","shape": [7, 6],"data": [[0,0,1.0],[0,1,2.0],[0,2,3.0],[0,3,4.0],[0,4,5.0],[0,5,6.0],[1,0,12.0],[1,1,11.0],[1,2,10.0],[1,3,9.0],[1,4,8.0],[1,5,7.0],[2,0,13.0],[2,1,14.0],[2,2,15.0],[2,3,16.0],[2,4,17.0],[2,5,18.0],[3,0,24.0],[3,1,23.0],[3,2,22.0],[3,3,21.0],[3,4,20.0],[3,5,19.0],[4,0,25.0],[4,1,26.0],[4,2,27.0],[4,3,28.0],[4,4,29.0],[4,5,30.0],[5,0,36.0],[5,1,35.0],[5,2,34.0],[5,3,33.0],[5,4,32.0],[5,5,31.0],[6,0,37.0],[6,1,38.0],[6,2,39.0],[6,3,40.0],[6,4,41.0],[6,5,42.0]],"rows": [{"id": "a", "metadata": null},{"id": "b", "metadata": null},{"id": "c", "metadata": null},{"id": "d", "metadata": null},{"id": "e", "metadata": null},{"id": "f", "metadata": null},{"id": "g", "metadata": null}],"columns": [{"id": "A", "metadata": null},{"id": "C", "metadata": null},{"id": "D", "metadata": null},{"id": "E", "metadata": null},{"id": "F", "metadata": null},{"id": "B", "metadata": null}]}'
    #     bt = parse_biom_table(bt_str)
    #     mf = ['#SampleIDt\thsid\tval\tph',
    #         'A\ta\tdummy\t-17.8',
    #         'B\tb\tdummy\t6.9',
    #         'C\ta\tdummy\t3.1',
    #         'D\tb\tdummy\t4.44',
    #         'E\ta\tdummy\t52.0',
    #         'F\ta\tdummy\t13.4']
    #     pmf, _ = parse_mapping_file_to_dict(mf)
    #     category = 'ph'
    #     hsid_to_samples = {'a':['A', 'C', 'E', 'F'], 'b':['B', 'D']}
    #     hsid_to_sample_indices = {'a':[0,1,3,4], 'b':[5, 2]}
    #     data = array([bt.observationData(i) for i in bt.ObservationIds])
    #     exp_otus = [[i.take([0,1,3,4]), i.take([5,2])] for i in data]
    #     exp_grad_vals = [array([-17.8, 3.1, 52.0, 13.4]), array([6.9, 4.44])]
    #     obs_data = longitudinal_row_generator(bt, pmf, category, 
    #         hsid_to_samples, hsid_to_sample_indices)
    #     for i,j in zip(obs_data, exp_otus):
    #         self.assertFloatEqual(i[0], j)
    #         self.assertFloatEqual(i[1], exp_grad_vals)

    # def test_run_longitudinal_correlation_test(self):
    #     """Test the longitudinal correlations are calculated correctly."""
    #     bt_str = '{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "testCode","date": "2013-08-22T18:23:28.199054","matrix_type": "sparse","matrix_element_type": "float","shape": [7, 6],"data": [[0,0,1.0],[0,1,2.0],[0,2,3.0],[0,3,4.0],[0,4,5.0],[0,5,6.0],[1,0,12.0],[1,1,11.0],[1,2,10.0],[1,3,9.0],[1,4,8.0],[1,5,7.0],[2,0,13.0],[2,1,14.0],[2,2,15.0],[2,3,16.0],[2,4,17.0],[2,5,18.0],[3,0,24.0],[3,1,23.0],[3,2,22.0],[3,3,21.0],[3,4,20.0],[3,5,19.0],[4,0,25.0],[4,1,26.0],[4,2,27.0],[4,3,28.0],[4,4,29.0],[4,5,30.0],[5,0,36.0],[5,1,35.0],[5,2,34.0],[5,3,33.0],[5,4,32.0],[5,5,31.0],[6,0,37.0],[6,1,38.0],[6,2,39.0],[6,3,40.0],[6,4,41.0],[6,5,42.0]],"rows": [{"id": "a", "metadata": null},{"id": "b", "metadata": null},{"id": "c", "metadata": null},{"id": "d", "metadata": null},{"id": "e", "metadata": null},{"id": "f", "metadata": null},{"id": "g", "metadata": null}],"columns": [{"id": "A", "metadata": null},{"id": "C", "metadata": null},{"id": "D", "metadata": null},{"id": "E", "metadata": null},{"id": "F", "metadata": null},{"id": "B", "metadata": null}]}'
    #     bt = parse_biom_table(bt_str)
    #     mf = ['#SampleIDt\thsid\tval\tph',
    #         'A\ta\tdummy\t-17.8',
    #         'B\tb\tdummy\t6.9',
    #         'C\ta\tdummy\t3.1',
    #         'D\tb\tdummy\t4.44',
    #         'E\ta\tdummy\t52.0',
    #         'F\ta\tdummy\t13.4']
    #     pmf, _ = parse_mapping_file_to_dict(mf)
    #     hsid_to_samples = {'a':['A', 'C', 'E', 'F'], 'b':['B', 'D']}
    #     hsid_to_sample_indices = {'a':[0,1,3,4], 'b':[5, 2]}
    #     data_feed = longitudinal_row_generator(bt, pmf, category, 
    #         hsid_to_samples, hsid_to_sample_indices)
    #     # test pearson
    #     exp = [[0.69462346989315615, 1.0],
    #         [-0.69462346989315593, -0.99999999999999445],
    #         [0.69462346989315671, 1.0],
    #         [-0.69462346989315593, -1.0],
    #         [0.69462346989315593, 0.99999999999998868],
    #         [-0.69462346989315737, -0.99999999999998868]]
    #     obs_test_stats = run_longitudinal_correlation_test(data_feed, 'pearson',
    #         CORRELATION_TEST_CHOICES)
    #     # test spearman
    #     data_feed = longitudinal_row_generator(bt, pmf, category, 
    #         hsid_to_samples, hsid_to_sample_indices)
    #     obs_test_stats = run_longitudinal_correlation_test(data_feed, 'spearman',
    #         CORRELATION_TEST_CHOICES)
    #     exp = [[0.80000000000000004, 1.0],
    #         [-0.80000000000000004, -1.0],
    #         [0.80000000000000004, 1.0],
    #         [-0.80000000000000004, -1.0],
    #         [0.80000000000000004, 1.0],
    #         [-0.80000000000000004, -1.0]]
    #     # test kendalls_tau
    #     data_feed = longitudinal_row_generator(bt, pmf, category, 
    #         hsid_to_samples, hsid_to_sample_indices)
    #     self.assertRaises(AssertionError, run_longitudinal_correlation_test, 
    #         data_feed, 'kendall', CORRELATION_TEST_CHOICES)
    #     hsid_to_samples = {'a':['A', 'C', 'F'], 'b':['B','E','D']}
    #     hsid_to_sample_indices = {'a':[0,1,4], 'b':[5, 3, 2]}
    #     data_feed = longitudinal_row_generator(bt, pmf, category, 
    #         hsid_to_samples, hsid_to_sample_indices)

    #     exp = [[0.6666666666666666, 1.0],
    #         [-0.6666666666666666, -1.0],
    #         [0.6666666666666666, 1.0],
    #         [-0.6666666666666666, -1.0],
    #         [0.6666666666666666, 1.0],
    #         [-0.6666666666666666, -1.0]]


    # def test_corerlation_row_generator(self):
    #     """correlation_row_generator works"""
    #     # test once Will updates with bt.iterObservationData()
    #     otu_table_1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0",
    #     "format_url": "http://biom-format.org","type": "OTU table","generated_by": 
    #     "BIOM-Format 1.1.2","date": "2013-08-16T10:16:20.131837","matrix_type": 
    #     "sparse","matrix_element_type": "float","shape": [6, 6],"data": [[0,0,28.0],
    #     [0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[1,0,25.0],[1,1,14.0],
    #     [1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],
    #     [2,3,69.0],[2,4,64.0],[2,5,27.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],
    #     [3,4,33.0],[3,5,62.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],
    #     [4,5,3.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0]],
    #     "rows": [{"id": "OTU1", "metadata": {"taxonomy": ["k__One"]}},{"id": "OTU2", 
    #     "metadata": {"taxonomy": ["k__Two"]}},{"id": "OTU3", "metadata": {"taxonomy": 
    #     ["k__Three"]}},{"id": "OTU4", "metadata": {"taxonomy": ["k__Four"]}},{"id": 
    #     "OTU5", "metadata": {"taxonomy": ["k__Five"]}},{"id": "OTU6", "metadata": 
    #     {"taxonomy": ["k__Six"]}}],"columns": [{"id": "Sample1", "metadata": null},
    #     {"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},
    #     {"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},
    #     {"id": "Sample6", "metadata": null}]}"""
    #     bt = parse_biom_table(otu_table_1) 
    #     pmf = {'Sample1': {'test_cat': 'cat1', 'test_corr': '1'},
    #      'Sample2': {'test_cat': 'cat1', 'test_corr': '2'},
    #      'Sample3': {'test_cat': 'cat2', 'test_corr': '3'},
    #      'Sample4': {'test_cat': 'cat2', 'test_corr': '4'},
    #      'Sample5': {'test_cat': 'cat3', 'test_corr': '5'},
    #      'Sample6': {'test_cat': 'cat3', 'test_corr': '6'}}
    #     data_result = []
    #     corr_row_gen_data = correlation_row_generator(bt, pmf, 'test_corr')
    #     for i in corr_row_gen_data:
    #         data_result.append(i)

    #     data_output = [(array([ 28.,  52.,  51.,  78.,  16.,  77.]),
    #       array([ 1.,  2.,  3.,  4.,  5.,  6.])),
    #      (array([ 25.,  14.,  11.,  32.,  48.,  63.]),
    #       array([ 1.,  2.,  3.,  4.,  5.,  6.])),
    #      (array([ 31.,   2.,  15.,  69.,  64.,  27.]),
    #       array([ 1.,  2.,  3.,  4.,  5.,  6.])),
    #      (array([ 36.,  68.,  70.,  65.,  33.,  62.]),
    #       array([ 1.,  2.,  3.,  4.,  5.,  6.])),
    #      (array([ 16.,  41.,  59.,  40.,  15.,   3.]),
    #       array([ 1.,  2.,  3.,  4.,  5.,  6.])),
    #      (array([ 32.,   8.,  54.,  98.,  29.,  50.]),
    #       array([ 1.,  2.,  3.,  4.,  5.,  6.]))]

    #     self.assertEqual(data_result, data_output)

    # def test_run_correlation_test_pearson(self):
    #     """run_correlation_test_pearson works"""
    #     test_choices = {'pearson': pearson}
    #     otu_table_1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0",
    #     "format_url": "http://biom-format.org","type": "OTU table","generated_by": 
    #     "BIOM-Format 1.1.2","date": "2013-08-16T10:16:20.131837","matrix_type": 
    #     "sparse","matrix_element_type": "float","shape": [6, 6],"data": [[0,0,28.0],
    #     [0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[1,0,25.0],[1,1,14.0],
    #     [1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],
    #     [2,3,69.0],[2,4,64.0],[2,5,27.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],
    #     [3,4,33.0],[3,5,62.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],
    #     [4,5,3.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0]],
    #     "rows": [{"id": "OTU1", "metadata": {"taxonomy": ["k__One"]}},{"id": "OTU2", 
    #     "metadata": {"taxonomy": ["k__Two"]}},{"id": "OTU3", "metadata": {"taxonomy": 
    #     ["k__Three"]}},{"id": "OTU4", "metadata": {"taxonomy": ["k__Four"]}},{"id": 
    #     "OTU5", "metadata": {"taxonomy": ["k__Five"]}},{"id": "OTU6", "metadata": 
    #     {"taxonomy": ["k__Six"]}}],"columns": [{"id": "Sample1", "metadata": null},
    #     {"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},
    #     {"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},
    #     {"id": "Sample6", "metadata": null}]}"""
    #     bt = parse_biom_table(otu_table_1) 
    #     pmf = {'Sample1': {'test_cat': 'cat1', 'test_corr': '1'},
    #      'Sample2': {'test_cat': 'cat1', 'test_corr': '2'},
    #      'Sample3': {'test_cat': 'cat2', 'test_corr': '3'},
    #      'Sample4': {'test_cat': 'cat2', 'test_corr': '4'},
    #      'Sample5': {'test_cat': 'cat3', 'test_corr': '5'},
    #      'Sample6': {'test_cat': 'cat3', 'test_corr': '6'}}

    #     data_gen = correlation_row_generator(bt, pmf, 'test_corr')

    #     corr_coeffs_output = [0.34884669332532803,
    #      0.83015306552396662,
    #      0.44037594794853846,
    #      0.064224958686639605,
    #      -0.41225244540969846,
    #      0.34313146667442163]
    #     p_pvals_output = [0.49795623745457107,
    #      0.04082210114420709,
    #      0.38213734667034305,
    #      0.90379502098011888,
    #      0.41665291191825937,
    #      0.50550281276602604]
    #     np_pvals_output = [0.48599999999999999,
    #      0.045999999999999999,
    #      0.38200000000000001,
    #      0.88200000000000001,
    #      0.438,
    #      0.54700000000000004]
    #     ci_highs_output = [0.90437105797013595,
    #      0.98087748393368879,
    #      0.92231069900787022,
    #      0.83239950278244235,
    #      0.60007469750146891,
    #      0.90318173092399967]
    #     ci_lows_output = [-0.64544754081743805,
    #      0.056981095244566363,
    #      -0.57762332473129863,
    #      -0.78843132325670628,
    #      -0.91701105702411645,
    #      -0.64921936409076253]

    #     seed(0)
    #     corr_coeffs_result, p_pvals_result, np_pvals_result, ci_highs_result, \
    #     ci_lows_result = run_correlation_test(data_gen, 'pearson', test_choices)

    #     self.assertFloatEqual(corr_coeffs_result, corr_coeffs_output)
    #     self.assertFloatEqual(p_pvals_result, p_pvals_output)
    #     self.assertFloatEqual(np_pvals_result, np_pvals_output)
    #     self.assertFloatEqual(ci_highs_result, ci_highs_output)
    #     self.assertFloatEqual(ci_lows_result, ci_lows_output)

    # def test_run_correlation_test_spearman(self):
    #     """run_correlation_test_spearman works"""
    #     test_choices = {'spearman': spearman}
    #     otu_table_1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0",
    #     "format_url": "http://biom-format.org","type": "OTU table","generated_by": 
    #     "BIOM-Format 1.1.2","date": "2013-08-16T10:16:20.131837","matrix_type": 
    #     "sparse","matrix_element_type": "float","shape": [6, 6],"data": [[0,0,28.0],
    #     [0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[1,0,25.0],[1,1,14.0],
    #     [1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],
    #     [2,3,69.0],[2,4,64.0],[2,5,27.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],
    #     [3,4,33.0],[3,5,62.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],
    #     [4,5,3.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0]],
    #     "rows": [{"id": "OTU1", "metadata": {"taxonomy": ["k__One"]}},{"id": "OTU2", 
    #     "metadata": {"taxonomy": ["k__Two"]}},{"id": "OTU3", "metadata": {"taxonomy": 
    #     ["k__Three"]}},{"id": "OTU4", "metadata": {"taxonomy": ["k__Four"]}},{"id": 
    #     "OTU5", "metadata": {"taxonomy": ["k__Five"]}},{"id": "OTU6", "metadata": 
    #     {"taxonomy": ["k__Six"]}}],"columns": [{"id": "Sample1", "metadata": null},
    #     {"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},
    #     {"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},
    #     {"id": "Sample6", "metadata": null}]}"""
    #     bt = parse_biom_table(otu_table_1) 
    #     pmf = {'Sample1': {'test_cat': 'cat1', 'test_corr': '1'},
    #      'Sample2': {'test_cat': 'cat1', 'test_corr': '2'},
    #      'Sample3': {'test_cat': 'cat2', 'test_corr': '3'},
    #      'Sample4': {'test_cat': 'cat2', 'test_corr': '4'},
    #      'Sample5': {'test_cat': 'cat3', 'test_corr': '5'},
    #      'Sample6': {'test_cat': 'cat3', 'test_corr': '6'}}

    #     data_gen = correlation_row_generator(bt, pmf, 'test_corr')

    #     corr_coeffs_output = [0.25714285714285712,
    #      0.77142857142857146,
    #      0.31428571428571428,
    #      -0.25714285714285712,
    #      -0.60000000000000009,
    #      0.25714285714285712]
    #     p_pvals_output = [0.62278717201166178,
    #      0.07239650145772597,
    #      0.54409329446064136,
    #      0.62278717201166178,
    #      0.20799999999999996,
    #      0.62278717201166178]
    #     np_pvals_output = [0.66200000000000003,
    #      0.085999999999999993,
    #      0.56299999999999994,
    #      0.65100000000000002,
    #      0.23100000000000001,
    #      0.65500000000000003]
    #     ci_highs_output = [0.88418587382825442,
    #      0.97351163763440296,
    #      0.89704483813983049,
    #      0.70063116821712201,
    #      0.41234932644798467,
    #      0.88418587382825442]
    #     ci_lows_output = [-0.70063116821712201,
    #      -0.10732436824253763,
    #      -0.66753963648898318,
    #      -0.88418587382825442,
    #      -0.94930820848352271,
    #      -0.70063116821712201]

    #     seed(0)
    #     corr_coeffs_result, p_pvals_result, np_pvals_result, ci_highs_result, \
    #     ci_lows_result = run_correlation_test(data_gen, 'spearman', test_choices)

    #     self.assertFloatEqual(corr_coeffs_result, corr_coeffs_output)
    #     self.assertFloatEqual(p_pvals_result, p_pvals_output)
    #     self.assertFloatEqual(np_pvals_result, np_pvals_output)
    #     self.assertFloatEqual(ci_highs_result, ci_highs_output)
    #     self.assertFloatEqual(ci_lows_result, ci_lows_output)

    # def test_run_correlation_test_kendall(self):
    #     """run_correlation_test_kendall works"""
    #     test_choices = {'kendall': kendall_correlation}
    #     otu_table_1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0",
    #     "format_url": "http://biom-format.org","type": "OTU table","generated_by": 
    #     "BIOM-Format 1.1.2","date": "2013-08-16T10:16:20.131837","matrix_type": 
    #     "sparse","matrix_element_type": "float","shape": [6, 6],"data": [[0,0,28.0],
    #     [0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[1,0,25.0],[1,1,14.0],
    #     [1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],
    #     [2,3,69.0],[2,4,64.0],[2,5,27.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],
    #     [3,4,33.0],[3,5,62.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],
    #     [4,5,3.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0]],
    #     "rows": [{"id": "OTU1", "metadata": {"taxonomy": ["k__One"]}},{"id": "OTU2", 
    #     "metadata": {"taxonomy": ["k__Two"]}},{"id": "OTU3", "metadata": {"taxonomy": 
    #     ["k__Three"]}},{"id": "OTU4", "metadata": {"taxonomy": ["k__Four"]}},{"id": 
    #     "OTU5", "metadata": {"taxonomy": ["k__Five"]}},{"id": "OTU6", "metadata": 
    #     {"taxonomy": ["k__Six"]}}],"columns": [{"id": "Sample1", "metadata": null},
    #     {"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},
    #     {"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},
    #     {"id": "Sample6", "metadata": null}]}"""
    #     bt = parse_biom_table(otu_table_1) 
    #     pmf = {'Sample1': {'test_cat': 'cat1', 'test_corr': '1'},
    #      'Sample2': {'test_cat': 'cat1', 'test_corr': '2'},
    #      'Sample3': {'test_cat': 'cat2', 'test_corr': '3'},
    #      'Sample4': {'test_cat': 'cat2', 'test_corr': '4'},
    #      'Sample5': {'test_cat': 'cat3', 'test_corr': '5'},
    #      'Sample6': {'test_cat': 'cat3', 'test_corr': '6'}}

    #     data_gen = correlation_row_generator(bt, pmf, 'test_corr')    

    #     corr_coeffs_output = [0.2, 0.6, 0.2, -0.2, -0.4666666666666667, 0.2]
    #     p_pvals_output = [0.7194444444444446,
    #      0.13611111111111107,
    #      0.7194444444444446,
    #      0.7194444444444444,
    #      0.2722222222222222,
    #      0.7194444444444446]
    #     np_pvals_output = [0.73199999999999998,
    #      0.14000000000000001,
    #      0.72599999999999998,
    #      0.70499999999999996,
    #      0.29899999999999999,
    #      0.71399999999999997]
    #     ci_highs_output = [0.87030079356638179,
    #      0.94930820848352271,
    #      0.87030079356638179,
    #      0.73005876315123075,
    #      0.55514322722082421,
    #      0.87030079356638179]
    #     ci_lows_output = [-0.73005876315123075,
    #      -0.41234932644798483,
    #      -0.73005876315123075,
    #      -0.87030079356638179,
    #      -0.92710628349420865,
    #      -0.73005876315123075]

    #     seed(0)
    #     corr_coeffs_result, p_pvals_result, np_pvals_result, ci_highs_result, \
    #     ci_lows_result = run_correlation_test(data_gen, 'kendall', test_choices)

    #     self.assertFloatEqual(corr_coeffs_result, corr_coeffs_output)
    #     self.assertFloatEqual(p_pvals_result, p_pvals_output)
    #     self.assertFloatEqual(np_pvals_result, np_pvals_output)
    #     self.assertFloatEqual(ci_highs_result, ci_highs_output)
    #     self.assertFloatEqual(ci_lows_result, ci_lows_output)

    # def test_correlation_output_formatter(self):
    #     """correlation_output_formatter works"""
    #     otu_table_1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0",
    #     "format_url": "http://biom-format.org","type": "OTU table","generated_by": 
    #     "BIOM-Format 1.1.2","date": "2013-08-16T10:16:20.131837","matrix_type": 
    #     "sparse","matrix_element_type": "float","shape": [6, 6],"data": [[0,0,28.0],
    #     [0,1,52.0],[0,2,51.0],[0,3,78.0],[0,4,16.0],[0,5,77.0],[1,0,25.0],[1,1,14.0],
    #     [1,2,11.0],[1,3,32.0],[1,4,48.0],[1,5,63.0],[2,0,31.0],[2,1,2.0],[2,2,15.0],
    #     [2,3,69.0],[2,4,64.0],[2,5,27.0],[3,0,36.0],[3,1,68.0],[3,2,70.0],[3,3,65.0],
    #     [3,4,33.0],[3,5,62.0],[4,0,16.0],[4,1,41.0],[4,2,59.0],[4,3,40.0],[4,4,15.0],
    #     [4,5,3.0],[5,0,32.0],[5,1,8.0],[5,2,54.0],[5,3,98.0],[5,4,29.0],[5,5,50.0]],
    #     "rows": [{"id": "OTU1", "metadata": {"taxonomy": ["k__One"]}},{"id": "OTU2", 
    #     "metadata": {"taxonomy": ["k__Two"]}},{"id": "OTU3", "metadata": {"taxonomy": 
    #     ["k__Three"]}},{"id": "OTU4", "metadata": {"taxonomy": ["k__Four"]}},{"id": 
    #     "OTU5", "metadata": {"taxonomy": ["k__Five"]}},{"id": "OTU6", "metadata": 
    #     {"taxonomy": ["k__Six"]}}],"columns": [{"id": "Sample1", "metadata": null},
    #     {"id": "Sample2", "metadata": null},{"id": "Sample3", "metadata": null},
    #     {"id": "Sample4", "metadata": null},{"id": "Sample5", "metadata": null},
    #     {"id": "Sample6", "metadata": null}]}"""
    #     bt = parse_biom_table(otu_table_1)

    #     corr_coeffs = [0.2, 0.6, 0.2, -0.2, -0.4666666666666667, 0.2]
    #     p_pvals = [0.7194444444444446,
    #      0.13611111111111107,
    #      0.7194444444444446,
    #      0.7194444444444444,
    #      0.2722222222222222,
    #      0.7194444444444446]
    #     p_pvals_fdr = [1.0791666666666668,
    #      0.8166666666666664,
    #      0.8633333333333334,
    #      1.4388888888888889,
    #      0.8166666666666667,
    #      0.7194444444444446]
    #     p_pvals_bon = [4.316666666666667,
    #      0.8166666666666664,
    #      4.316666666666667,
    #      4.316666666666666,
    #      1.6333333333333333,
    #      4.316666666666667]
    #     np_pvals = [0.73199999999999998,
    #      0.14000000000000001,
    #      0.72599999999999998,
    #      0.70499999999999996,
    #      0.29899999999999999,
    #      0.71399999999999997]
    #     np_pvals_fdr = [0.73199999999999998,
    #      0.84000000000000008,
    #      0.87119999999999997,
    #      1.4099999999999999,
    #      0.89700000000000002,
    #      1.071]
    #     np_pvals_bon = [4.3919999999999995,
    #      0.8400000000000001,
    #      4.356,
    #      4.2299999999999995,
    #      1.794,
    #      4.284]
    #     ci_highs = [0.87030079356638179,
    #      0.94930820848352271,
    #      0.87030079356638179,
    #      0.73005876315123075,
    #      0.55514322722082421,
    #      0.87030079356638179]
    #     ci_lows = [-0.73005876315123075,
    #      -0.41234932644798483,
    #      -0.73005876315123075,
    #      -0.87030079356638179,
    #      -0.92710628349420865,
    #      -0.73005876315123075]

        

    #     lines_result = correlation_output_formatter(bt, corr_coefs, p_pvals, p_pvals_fdr, \
    #         p_vals_bon, np_pvals, np_pvals_fdr, np_pvals_bon, ci_highs, \
    #         ci_lows)



#run unit tests if run from command-line
if __name__ == '__main__':
    main()













    
