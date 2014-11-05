#!/usr/bin/env python
from __future__ import division

__author__ = "Michael Dwan"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Michael Dwan", "Logan Knecht",
               "Damien Coy", "Levi McCracken", "Andrew Cochran",
               "Will Van Treuren"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Test suite for classes, methods and functions of the stats module."""

from shutil import rmtree
from os.path import exists, join
from string import digits
from tempfile import mkdtemp
from StringIO import StringIO
from unittest import TestCase, main
from warnings import filterwarnings
from itertools import izip
from types import StringType, ListType, FloatType, TupleType


from skbio.util import remove_files
from numpy.testing import assert_almost_equal, assert_allclose
from numpy import (array, asarray, roll, median, nan, arange, matrix,
                   concatenate, nan, ndarray, number, ones,
                   reshape, testing, tril, var, log, fill_diagonal)
from numpy.random import permutation, shuffle, seed
from biom import Table
from biom.util import biom_open

from qiime.stats import (all_pairs_t_test, _perform_pairwise_tests,
                         CorrelationStats,
                         DistanceMatrixStats, MantelCorrelogram, Mantel,
                         PartialMantel, quantile, _quantile,
                         paired_difference_analyses,
                         G_2_by_2, g_fit, t_paired, t_one_sample,
                         t_two_sample, mc_t_two_sample,
                         _permute_observations,
                         correlation_t, ZeroExpectedError, fisher,
                         safe_sum_p_log_p, permute_2d, mantel,
                         mantel_t, _flatten_lower_triangle, pearson,
                         spearman, ANOVA_one_way, mw_t,
                         mw_boot, is_symmetric_and_hollow,
                         tail, fdr_correction,
                         benjamini_hochberg_step_down,
                         bonferroni_correction, fisher_z_transform,
                         fisher_population_correlation,
                         inverse_fisher_z_transform,
                         z_transform_pval, kruskal_wallis, kendall,
                         kendall_pval, assign_correlation_pval,
                         cscore, williams_correction, t_one_observation,
                         normprob, tprob, fprob, chi2prob)

from skbio.stats.distance import (DissimilarityMatrix, DistanceMatrix)

from qiime.util import MetadataMap, get_qiime_temp_dir


class TestHelper(TestCase):

    """Helper class that instantiates some commonly-used objects.

    This class should be subclassed by any test classes that want to use its
    members.
    """

    def compare_multiple_level_array(self, observed, expected):
        """ Compare multiple level arrays.

        It expecte observed and expected arrays, where each element is an
        array of elements.
        """
        if isinstance(observed, (TupleType, ListType)):
            for obs, exp in izip(observed, expected):
                self.compare_multiple_level_array(obs, exp)
        elif observed is not None and isinstance(observed, (number, ndarray, FloatType)):
            assert_almost_equal(observed, expected, decimal=5)
        else:
            self.assertEqual(observed, expected)

    def setUp(self):
        """Define some useful test objects."""
        # The unweighted unifrac distance matrix from the overview tutorial.
        self.overview_dm_str = ["\tPC.354\tPC.355\tPC.356\tPC.481\tPC.593\
                                \tPC.607\tPC.634\tPC.635\tPC.636",
                                "PC.354\t0.0\t0.595483768391\t0.618074717633\
                                \t0.582763100909\t0.566949022108\
                                \t0.714717232268\t0.772001731764\
                                \t0.690237118413\t0.740681707488",
                                "PC.355\t0.595483768391\t0.0\t0.581427669668\
                                \t0.613726772383\t0.65945132763\
                                \t0.745176523638\t0.733836123821\
                                \t0.720305073505\t0.680785600439",
                                "PC.356\t0.618074717633\t0.581427669668\t0.0\
                                \t0.672149021573\t0.699416863323\
                                \t0.71405573754\t0.759178215168\
                                \t0.689701276341\t0.725100672826",
                                "PC.481\t0.582763100909\t0.613726772383\
                                \t0.672149021573\t0.0\t0.64756120797\
                                \t0.666018240373\t0.66532968784\
                                \t0.650464714994\t0.632524644216",
                                "PC.593\t0.566949022108\t0.65945132763\
                                \t0.699416863323\t0.64756120797\t0.0\
                                \t0.703720200713\t0.748240937349\
                                \t0.73416971958\t0.727154987937",
                                "PC.607\t0.714717232268\t0.745176523638\
                                \t0.71405573754\t0.666018240373\
                                \t0.703720200713\t0.0\t0.707316869557\
                                \t0.636288883818\t0.699880573956",
                                "PC.634\t0.772001731764\t0.733836123821\
                                \t0.759178215168\t0.66532968784\
                                \t0.748240937349\t0.707316869557\t0.0\
                                \t0.565875193399\t0.560605525642",
                                "PC.635\t0.690237118413\t0.720305073505\
                                \t0.689701276341\t0.650464714994\
                                \t0.73416971958\t0.636288883818\
                                \t0.565875193399\t0.0\t0.575788039321",
                                "PC.636\t0.740681707488\t0.680785600439\
                                \t0.725100672826\t0.632524644216\
                                \t0.727154987937\t0.699880573956\
                                \t0.560605525642\t0.575788039321\t0.0"]
        self.overview_dm = DistanceMatrix.read(\
            StringIO('\n'.join(self.overview_dm_str)))

        # The overview tutorial's metadata mapping file.
        self.overview_map_str = ["#SampleID\tBarcodeSequence\tTreatment\tDOB",
                                 "PC.354\tAGCACGAGCCTA\tControl\t20061218",
                                 "PC.355\tAACTCGTCGATG\tControl\t20061218",
                                 "PC.356\tACAGACCACTCA\tControl\t20061126",
                                 "PC.481\tACCAGCGACTAG\tControl\t20070314",
                                 "PC.593\tAGCAGCACTTGT\tControl\t20071210",
                                 "PC.607\tAACTGTGCGTAC\tFast\t20071112",
                                 "PC.634\tACAGAGTCGGCT\tFast\t20080116",
                                 "PC.635\tACCGCAGAGTCA\tFast\t20080116",
                                 "PC.636\tACGGTGAGTGTC\tFast\t20080116"]
        self.overview_map = MetadataMap.parseMetadataMap(self.overview_map_str)

        self.test_map_str = [
            "#SampleID\tBarcodeSequence\tFoo\tBar\tDescription",
            "PC.354\tAGCACGAGCCTA\tfoo\ta\t354",
            "PC.355\tAACTCGTCGATG\tfoo\ta\t355",
            "PC.356\tACAGACCACTCA\tbar\ta\t356",
            "PC.481\tACCAGCGACTAG\tfoo\ta\t481",
            "PC.593\tAGCAGCACTTGT\tbar\ta\t593",
            "PC.607\tAACTGTGCGTAC\tbar\ta\t607",
            "PC.634\tACAGAGTCGGCT\tbar\ta\t634",
            "PC.635\tACCGCAGAGTCA\tfoo\ta\t635",
            "PC.636\tACGGTGAGTGTC\tbar\ta\t636"]
        self.test_map = MetadataMap.parseMetadataMap(self.test_map_str)

        # A 1x1 dm.
        self.single_ele_dm = DistanceMatrix([[0]], ['s1'])

        # How many times to test a p-value.
        self.p_val_tests = 10

    def assertCorrectPValue(self, exp_min, exp_max, fn, num_perms=None,
                            p_val_key='p_value'):
        """Tests that the stochastic p-value falls in the specified range.

        Performs the test self.p_val_tests times and fails if the observed
        p-value does not fall into the specified range at least once. Each
        p-value is also tested that it falls in the range 0.0 to 1.0.

        This method assumes that fn is callable, and will pass num_perms to fn
        if num_perms is provided. p_val_key specifies the key that will be used
        to retrieve the p-value from the results dict that is returned by fn.
        """
        found_match = False
        for i in range(self.p_val_tests):
            if num_perms is not None:
                obs = fn(num_perms)
            else:
                obs = fn()
            p_val = obs[p_val_key]
            self.assertTrue(0.0 <= p_val < 1.0)
            if p_val >= exp_min and p_val <= exp_max:
                found_match = True
                break
        self.assertTrue(found_match)


class NonRandomShuffler(object):

    """Helper class for testing p-values that are calculated by permutations.

    Since p-values rely on randomness, it may be useful to use a non-random
    function (such as that provided by this class) to generate permutations
    so that p-values can be accurately tested.

    This code is heavily based on Andrew Cochran's original version.
    """

    def __init__(self):
        """Default constructor initializes the number of calls to zero."""
        self.num_calls = 0

    def permutation(self, x):
        """Non-random permutation function to test p-test code.

        Returns the 'permuted' version of x.

        Arguments:
            x - the array to be 'permuted'
        """
        x = array(x)
        x = roll(x, self.num_calls)
        self.num_calls += 1
        return x


class StatsTests(TestHelper):

    """Tests for top-level functions in the stats module."""

    def setUp(self):
        """Set up data that will be used by the tests."""
        self.value_for_seed = 20

        # Single comp.
        self.labels1 = ['foo', 'bar']
        self.dists1 = [[1, 2, 3], [7, 8]]

        # Multiple comps.
        self.labels2 = ['foo', 'bar', 'baz']
        self.dists2 = [[1, 2, 3], [7, 8], [9, 10, 11]]

        # Too few obs.
        self.labels3 = ['foo', 'bar', 'baz']
        self.dists3 = [[1], [7], [9, 10, 11]]

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

    def test_all_pairs_t_test(self):
        """Test performing Monte Carlo tests on valid dataset."""
        # We aren't testing the numeric values here, as they've already been
        # tested in the functions that compute them. We are interested in the
        # format of the returned string.
        exp = """# The tests of significance were performed using a two-sided Student's two-sample t-test.
# Alternative hypothesis: Group 1 mean != Group 2 mean
# The nonparametric p-values were calculated using 999 Monte Carlo permutations.
# The nonparametric p-values contain the correct number of significant digits.
# Entries marked with "N/A" could not be calculated because at least one of the groups
# of distances was empty, both groups each contained only a single distance, or
# the test could not be performed (e.g. no variance in groups with the same mean).
Group 1	Group 2	t statistic	Parametric p-value	Parametric p-value (Bonferroni-corrected)	Nonparametric p-value	Nonparametric p-value (Bonferroni-corrected)
foo	bar	-6.6	0.00708047956412	0.0212414386924	0.095	0.285
foo	baz	-9.79795897113	0.000608184944463	0.00182455483339	0.101	0.303
bar	baz	-3.0	0.0576688856224	0.173006656867	0.217	0.651
"""
        obs = all_pairs_t_test(self.labels2, self.dists2)
        self.assertEqual(self.remove_nums(obs), self.remove_nums(exp))

    def test_all_pairs_t_test_no_perms(self):
        """Test performing Monte Carlo tests on valid dataset with no perms."""
        exp = """# The tests of significance were performed using a two-sided Student's two-sample t-test.
# Alternative hypothesis: Group 1 mean != Group 2 mean
# Entries marked with "N/A" could not be calculated because at least one of the groups
# of distances was empty, both groups each contained only a single distance, or
# the test could not be performed (e.g. no variance in groups with the same mean).
Group 1	Group 2	t statistic	Parametric p-value	Parametric p-value (Bonferroni-corrected)	Nonparametric p-value	Nonparametric p-value (Bonferroni-corrected)
foo	bar	-6.6	0.00708047956412	0.0212414386924	N/A	N/A
foo	baz	-9.79795897113	0.000608184944463	0.00182455483339	N/A	N/A
bar	baz	-3.0	0.0576688856224	0.173006656867	N/A	N/A
"""
        obs = all_pairs_t_test(self.labels2, self.dists2,
                               num_permutations=0)
        self.assertEqual(self.remove_nums(obs), self.remove_nums(exp))

    def test_all_pairs_t_test_few_perms(self):
        """Test performing Monte Carlo tests on dataset with a few perms."""
        exp = """# The tests of significance were performed using a one-sided (low) Student's two-sample t-test.
# Alternative hypothesis: Group 1 mean < Group 2 mean
# The nonparametric p-values were calculated using 5 Monte Carlo permutations.
# The nonparametric p-values contain the correct number of significant digits.
# Entries marked with "N/A" could not be calculated because at least one of the groups
# of distances was empty, both groups each contained only a single distance, or
# the test could not be performed (e.g. no variance in groups with the same mean).
Group 1	Group 2	t statistic	Parametric p-value	Parametric p-value (Bonferroni-corrected)	Nonparametric p-value	Nonparametric p-value (Bonferroni-corrected)
foo	bar	-6.6	0.00354023978206	0.0106207193462	Too few iters to compute p-value (num_iters=5)	Too few iters to compute p-value (num_iters=5)
foo	baz	-9.79795897113	0.000304092472232	0.000912277416695	Too few iters to compute p-value (num_iters=5)	Too few iters to compute p-value (num_iters=5)
bar	baz	-3.0	0.0288344428112	0.0865033284337	Too few iters to compute p-value (num_iters=5)	Too few iters to compute p-value (num_iters=5)
"""
        obs = all_pairs_t_test(self.labels2, self.dists2,
                               num_permutations=5, tail_type='low')
        self.assertEqual(self.remove_nums(obs), self.remove_nums(exp))

    def test_all_pairs_t_test_invalid_tests(self):
        """Test performing Monte Carlo tests with some invalid tests."""
        exp = """# The tests of significance were performed using a one-sided (high) Student's two-sample t-test.
# Alternative hypothesis: Group 1 mean > Group 2 mean
# The nonparametric p-values were calculated using 20 Monte Carlo permutations.
# The nonparametric p-values contain the correct number of significant digits.
# Entries marked with "N/A" could not be calculated because at least one of the groups
# of distances was empty, both groups each contained only a single distance, or
# the test could not be performed (e.g. no variance in groups with the same mean).
Group 1	Group 2	t statistic	Parametric p-value	Parametric p-value (Bonferroni-corrected)	Nonparametric p-value	Nonparametric p-value (Bonferroni-corrected)
foo	bar	N/A	N/A	N/A	N/A	N/A
"""
        obs = all_pairs_t_test(['foo', 'bar'], [[], [1, 2, 4]],
                               'high', 20)
        self.assertEqual(self.remove_nums(obs), self.remove_nums(exp))

    def test_all_pairs_t_test_invalid_input(self):
        """Test performing Monte Carlo tests on invalid input."""
        # Number of labels and distance groups do not match.
        self.assertRaises(ValueError, all_pairs_t_test,
                          ['foo', 'bar'], [[1, 2, 3], [4, 5, 6], [7, 8]])

        # Invalid tail type.
        self.assertRaises(ValueError, all_pairs_t_test,
                          ['foo', 'bar'], [[1, 2, 3], [4, 5, 6]], 'foo')

        # Invalid number of permutations.
        self.assertRaises(ValueError, all_pairs_t_test,
                          ['foo', 'bar'], [[1, 2, 3], [4, 5, 6]], num_permutations=-1)

    def test_perform_pairwise_tests_single_comp(self):
        """Test on valid dataset w/ 1 comp."""
        # Verified with R's t.test function.
        exp = [['foo', 'bar', -6.5999999999999996, 0.0070804795641244006,
                0.0070804795641244006, 0.100000000001, 0.10000000000001]]
        seed(self.value_for_seed)
        obs = _perform_pairwise_tests(self.labels1, self.dists1, 'two-sided',
                                      999)
        self.compare_multiple_level_array(obs, exp)

    def test_perform_pairwise_tests_multi_comp(self):
        """Test on valid dataset w/ multiple comps."""
        # Verified with R's t.test function.
        exp = [['foo', 'bar', -6.5999999999999996, 0.0070804795641244006,
                0.021241438692373202, nan, nan], ['foo', 'baz',
                                                  -
                                                  9.7979589711327115, 0.00060818494446333643, 0.0018245548333900093,
                                                  nan, nan], ['bar', 'baz', -3.0, 0.05766888562243732,
                                                              0.17300665686731195, nan, nan]]
        obs = _perform_pairwise_tests(self.labels2, self.dists2, 'two-sided',
                                      0)
        self.compare_multiple_level_array(obs, exp)

    def test_perform_pairwise_tests_too_few_obs(self):
        """Test on dataset w/ too few observations."""
        exp = [['foo', 'bar', nan, nan, nan, nan, nan],
               ['foo', 'baz', -7.794228634059948, 0.008032650971672552,
                0.016065301943345104, nan, nan],
               ['bar', 'baz', -2.598076211353316, 0.060844967173160069,
                0.12168993434632014, nan, nan]]

        obs = _perform_pairwise_tests(self.labels3, self.dists3, 'low', 0)
        self.compare_multiple_level_array(obs, exp)

        exp = [['foo', 'bar', nan, nan, nan, nan, nan]]
        obs = _perform_pairwise_tests(['foo', 'bar'], [[], [1, 2, 4]], 'high',
                                      20)
        self.compare_multiple_level_array(obs, exp)


class DistanceMatrixStatsTests(TestHelper):

    """Tests for the DistanceMatrixStats class."""

    def setUp(self):
        """Define some dm stats instances that will be used by the tests."""
        super(DistanceMatrixStatsTests, self).setUp()

        self.empty_dms = DistanceMatrixStats([])
        self.single_dms = DistanceMatrixStats([self.overview_dm])
        self.double_dms = DistanceMatrixStats(
            [self.overview_dm, self.single_ele_dm])
        # For testing the requirement that two distance matrices are set.
        self.two_dms = DistanceMatrixStats(
            [self.overview_dm, self.single_ele_dm], 2)
        # For testing the requirement that the distance matrices meet the
        # minimum size requirements.
        self.size_dms = DistanceMatrixStats(
            [self.overview_dm, self.overview_dm], 2, 4)

    def test_DistanceMatrices_getter(self):
        """Test getter for distmats."""
        self.assertEqual(self.empty_dms.DistanceMatrices, [])
        self.assertEqual(self.single_dms.DistanceMatrices, [self.overview_dm])
        self.assertEqual(self.double_dms.DistanceMatrices,
                         [self.overview_dm, self.single_ele_dm])

    def test_DistanceMatrices_setter(self):
        """Test setter for dms on valid input data."""
        self.empty_dms.DistanceMatrices = []
        self.assertEqual(self.empty_dms.DistanceMatrices, [])

        self.empty_dms.DistanceMatrices = [self.overview_dm]
        self.assertEqual(self.empty_dms.DistanceMatrices, [self.overview_dm])

        self.empty_dms.DistanceMatrices = [self.overview_dm, self.overview_dm]
        self.assertEqual(self.empty_dms.DistanceMatrices,
                         [self.overview_dm, self.overview_dm])

    def test_DistanceMatrices_setter_invalid(self):
        """Test setter for dms on invalid input data."""
        # Allows testing of non-callable property setter that raises errors.
        # Idea was obtained from http://stackoverflow.com/a/3073049
        self.assertRaises(TypeError, setattr, self.empty_dms,
                          'DistanceMatrices', None)
        self.assertRaises(TypeError, setattr, self.empty_dms,
                          'DistanceMatrices', 10)
        self.assertRaises(TypeError, setattr, self.empty_dms,
                          'DistanceMatrices', 20.0)
        self.assertRaises(TypeError, setattr, self.empty_dms,
                          'DistanceMatrices', "foo")
        self.assertRaises(TypeError, setattr, self.empty_dms,
                          'DistanceMatrices', {})
        self.assertRaises(TypeError, setattr, self.empty_dms,
                          'DistanceMatrices', self.overview_dm)
        self.assertRaises(TypeError, setattr, self.empty_dms,
                          'DistanceMatrices', [1])
        self.assertRaises(TypeError, setattr, self.empty_dms,
                          'DistanceMatrices',
                          [DissimilarityMatrix(
                              array([[0, 2], [3, 0]]), ['foo', 'bar']),
                           DissimilarityMatrix(
                               array([[0, 2], [3.5, 0]]), ['foo', 'bar'])])

        # Test constructor as well.
        self.assertRaises(TypeError, DistanceMatrixStats, None)
        self.assertRaises(TypeError, DistanceMatrixStats, 10)
        self.assertRaises(TypeError, DistanceMatrixStats, 20.0)
        self.assertRaises(TypeError, DistanceMatrixStats, "foo")
        self.assertRaises(TypeError, DistanceMatrixStats, {})
        self.assertRaises(TypeError, DistanceMatrixStats, self.overview_dm)
        self.assertRaises(TypeError, DistanceMatrixStats, [1])
        self.assertRaises(TypeError, DistanceMatrixStats,
                          [DissimilarityMatrix(
                              array([[0, 2], [3, 0]]), ['foo', 'bar']),
                           DissimilarityMatrix(
                               array([[0, 2], [3.5, 0]]), ['foo', 'bar'])])

    def test_DistanceMatrices_setter_wrong_number(self):
        """Test setting an invalid number of distance matrices."""
        self.assertRaises(ValueError, setattr, self.two_dms,
                          'DistanceMatrices', [self.overview_dm])
        self.assertRaises(ValueError, setattr, self.two_dms,
                          'DistanceMatrices', [self.overview_dm, self.overview_dm,
                                               self.overview_dm])

    def test_DistanceMatrices_setter_too_small(self):
        """Test setting distance matrices that are too small."""
        self.assertRaises(ValueError, setattr, self.size_dms,
                          'DistanceMatrices', [self.single_ele_dm, self.single_ele_dm])

    def test_call(self):
        """Test __call__() returns an empty result set."""
        self.assertEqual(self.single_dms(), {})
        self.assertEqual(self.single_dms(10), {})
        self.assertEqual(self.single_dms(0), {})

    def test_call_bad_perms(self):
        """Test __call__() fails upon receiving invalid number of perms."""
        self.assertRaises(ValueError, self.single_dms, -1)


class CorrelationStatsTests(TestHelper):

    """Tests for the CorrelationStats class."""

    def setUp(self):
        """Set up correlation stats instances for use in tests."""
        super(CorrelationStatsTests, self).setUp()
        self.cs = CorrelationStats([self.overview_dm, self.overview_dm])

    def test_DistanceMatrices_setter(self):
        """Test setting valid distance matrices."""
        dms = [self.overview_dm, self.overview_dm]
        self.cs.DistanceMatrices = dms
        self.assertEqual(self.cs.DistanceMatrices, dms)

        dms = [self.overview_dm, self.overview_dm, self.overview_dm]
        self.cs.DistanceMatrices = dms
        self.assertEqual(self.cs.DistanceMatrices, dms)

    def test_DistanceMatrices_setter_mismatched_labels(self):
        """Test setting dms with mismatching sample ID labels."""
        mismatch = DistanceMatrix(array([[0]]), ['s2'])

        self.assertRaises(ValueError, setattr, self.cs, 'DistanceMatrices',
                          [self.single_ele_dm, mismatch])
        # Also test that constructor raises this error.
        self.assertRaises(ValueError, CorrelationStats, [self.single_ele_dm,
                                                         mismatch])

    def test_DistanceMatrices_setter_wrong_dims(self):
        """Test setting dms with mismatching dimensions."""
        self.assertRaises(ValueError, setattr, self.cs, 'DistanceMatrices',
                          [self.overview_dm, self.single_ele_dm])
        # Also test that constructor raises this error.
        self.assertRaises(ValueError, CorrelationStats, [self.overview_dm,
                                                         self.single_ele_dm])

    def test_DistanceMatrices_setter_too_few(self):
        """Test setting dms with not enough of them."""
        self.assertRaises(ValueError, setattr, self.cs, 'DistanceMatrices', [])
        # Also test that constructor raises this error.
        self.assertRaises(ValueError, CorrelationStats, [])

    def test_call(self):
        """Test __call__() returns an empty result set."""
        self.assertEqual(self.cs(), {})


class MantelCorrelogramTests(TestHelper):
    """Tests for the MantelCorrelogram class."""

    def setUp(self):
        """Set up mantel correlogram instances for use in tests."""
        super(MantelCorrelogramTests, self).setUp()

        # Mantel correlogram test using the overview tutorial's unifrac dm as
        # both inputs.
        self.mc = MantelCorrelogram(self.overview_dm, self.overview_dm)

        # Smallest test case: 3x3 matrices.
        ids = ['s1', 's2', 's3']
        dm1 = DistanceMatrix(array([[0, 1, 2], [1, 0, 3], [2, 3, 0]]), ids)
        dm2 = DistanceMatrix(array([[0, 2, 5], [2, 0, 8], [5, 8, 0]]), ids)

        self.small_mc = MantelCorrelogram(dm1, dm2)

        # For testing variable-sized bins.
        self.small_mc_var_bins = MantelCorrelogram(dm1, dm2,
                                                   variable_size_distance_classes=True)

    def test_Alpha_getter(self):
        """Test retrieving the value of alpha."""
        self.assertEqual(self.mc.Alpha, 0.05)

    def test_Alpha_setter(self):
        """Test setting the value of alpha."""
        self.mc.Alpha = 0.01
        self.assertEqual(self.mc.Alpha, 0.01)

    def test_Alpha_setter_invalid(self):
        """Test setting the value of alpha with an invalid value."""
        self.assertRaises(ValueError, setattr, self.mc, 'Alpha', -5)
        self.assertRaises(ValueError, setattr, self.mc, 'Alpha', 2)

    def test_DistanceMatrices_setter(self):
        """Test setting a valid number of distance matrices."""
        dms = [self.overview_dm, self.overview_dm]
        self.mc.DistanceMatrices = dms
        self.assertEqual(self.mc.DistanceMatrices, dms)

    def test_DistanceMatrices_setter_wrong_number(self):
        """Test setting an invalid number of distance matrices."""
        self.assertRaises(ValueError, setattr, self.mc, 'DistanceMatrices',
                          [self.overview_dm])
        self.assertRaises(ValueError, setattr, self.mc, 'DistanceMatrices',
                          [self.overview_dm, self.overview_dm, self.overview_dm])

    def test_DistanceMatrices_setter_too_small(self):
        """Test setting distance matrices that are too small."""
        self.assertRaises(ValueError, setattr, self.mc, 'DistanceMatrices',
                          [self.single_ele_dm, self.single_ele_dm])

    def test_call(self):
        """Test running a Mantel correlogram analysis on valid input."""
        # A lot of the returned numbers are based on random permutations and
        # thus cannot be tested for exact values. We'll test what we can
        # exactly, and then test for "sane" values for the "random" values. The
        # matplotlib Figure object cannot be easily tested either, so we'll try
        # our best to make sure it appears sane.
        obs = self.mc()

        exp_method_name = 'Mantel Correlogram'
        self.assertEqual(obs['method_name'], exp_method_name)

        exp_class_index = [0.5757052546507142, 0.60590471266814283,
                           0.63610417068557146, 0.66630362870299997, 0.69650308672042849,
                           0.72670254473785723, 0.75690200275528574]
        assert_almost_equal(obs['class_index'], exp_class_index)

        exp_num_dist = [12, 6, 8, 10, 12, 16, 8]
        self.assertEqual(obs['num_dist'], exp_num_dist)

        exp_mantel_r = [0.73244729118260765, 0.31157641757444593,
                        0.17627427296718071, None, None, None, None]
        self.compare_multiple_level_array(obs['mantel_r'], exp_mantel_r)

        # Test matplotlib Figure for a sane state.
        obs_fig = obs['correlogram_plot']
        obs_ax = obs_fig.get_axes()[0]
        self.assertEqual(obs_ax.get_title(), "Mantel Correlogram")
        self.assertEqual(obs_ax.get_xlabel(), "Distance class index")
        self.assertEqual(obs_ax.get_ylabel(), "Mantel correlation statistic")
        assert_almost_equal(obs_ax.get_xticks(), [0.57, 0.58, 0.59, 0.6,
                                                  0.61, 0.62, 0.63, 0.64, 0.65])
        assert_almost_equal(obs_ax.get_yticks(), [0.1, 0.2, 0.3, 0.4, 0.5,
                                                  0.6, 0.7, 0.8, 0.9])

        # Test p-values and corrected p-values.
        found_match = False
        for i in range(self.p_val_tests):
            obs = self.mc()
            p_vals = obs['mantel_p']
            corr_p_vals = obs['mantel_p_corr']
            self.assertEqual(len(p_vals), 7)
            self.assertEqual(p_vals[3:], [None, None, None, None])
            self.assertTrue(0.0 <= p_vals[0] <= 1.0)
            self.assertTrue(0.0 <= p_vals[1] <= 1.0)
            self.assertTrue(0.0 <= p_vals[2] <= 1.0)
            self.compare_multiple_level_array(corr_p_vals,
                                              [p_val * 3 if p_val is not None else None for p_val in p_vals])

            if (p_vals[0] >= 0 and p_vals[0] <= 0.01 and p_vals[1] > 0.01 and
                    p_vals[1] <= 0.1 and p_vals[2] > 0.1 and p_vals[2] <= 0.5):
                found_match = True
                break
        self.assertTrue(found_match)

    def test_call_small(self):
        """Test running a Mantel correlogram analysis on the smallest input."""
        # The expected output was verified with vegan's mantel correlogram
        # function.
        obs = self.small_mc()

        exp_method_name = 'Mantel Correlogram'
        self.assertEqual(obs['method_name'], exp_method_name)

        exp_class_index = [3.0, 5.0, 7.0]
        assert_almost_equal(obs['class_index'], exp_class_index)

        exp_num_dist = [2, 2, 2]
        self.assertEqual(obs['num_dist'], exp_num_dist)

        exp_mantel_r = [0.86602540378443871, None, None]
        self.compare_multiple_level_array(obs['mantel_r'], exp_mantel_r)

        # Test matplotlib Figure for a sane state.
        obs_fig = obs['correlogram_plot']
        obs_ax = obs_fig.get_axes()[0]
        self.assertEqual(obs_ax.get_title(), "Mantel Correlogram")
        self.assertEqual(obs_ax.get_xlabel(), "Distance class index")
        self.assertEqual(obs_ax.get_ylabel(), "Mantel correlation statistic")
        assert_almost_equal(obs_ax.get_xticks(), [2.85, 2.9, 2.95, 3., 3.05,
                                                  3.1, 3.15, 3.2])
        assert_almost_equal(obs_ax.get_yticks(), [0.82, 0.83, 0.84, 0.85,
                                                  0.86, 0.87, 0.88, 0.89, 0.9, 0.91])

        # Test p-values and corrected p-values.
        found_match = False
        for i in range(self.p_val_tests):
            obs = self.small_mc()
            p_vals = obs['mantel_p']
            corr_p_vals = obs['mantel_p_corr']
            self.assertEqual(len(p_vals), 3)
            self.assertEqual(p_vals[1:], [None, None])
            self.assertTrue(0.0 <= p_vals[0] <= 1.0)
            self.compare_multiple_level_array(corr_p_vals, p_vals)

            if p_vals[0] >= 0 and p_vals[0] <= 0.5:
                found_match = True
                break
        self.assertTrue(found_match)

    def test_find_distance_classes(self):
        """Test finding the distance classes a matrix's elements are in."""
        exp = (array([[-1, 0, 1], [0, -1, 2], [1, 2, -1]]),
               [3.0, 5.0, 7.0])
        obs = self.small_mc._find_distance_classes(
            self.small_mc.DistanceMatrices[1], 3)
        self.compare_multiple_level_array(obs, exp)

        exp = (array([[-1, 1, 2, 0, 0, 5, 7, 4, 6],
                      [1, -1, 0, 2, 3, 6, 6, 6, 4],
                      [2, 0, -1, 4, 5, 5, 7, 4, 6],
                      [0, 2, 4, -1, 3, 3, 3, 3, 2],
                      [0, 3, 5, 3, -1, 5, 7, 6, 6],
                      [5, 6, 5, 3, 5, -1, 5, 2, 5],
                      [7, 6, 7, 3, 7, 5, -1, 0, 0],
                      [4, 6, 4, 3, 6, 2, 0, -1, 0],
                      [6, 4, 6, 2, 6, 5, 0, 0, -1]]),
               [0.57381779, 0.60024231, 0.62666684, 0.65309137, 0.67951589,
                0.70594042, 0.73236494, 0.75878947])
        obs = self.mc._find_distance_classes(
            self.mc.DistanceMatrices[1], 8)
        self.compare_multiple_level_array(obs, exp)

    def test_find_distance_classes_variable_size_bins(self):
        """Test finding distance classes with variable-size bins."""
        # Single distance class.
        exp = (array([[-1, 0, 0], [0, -1, 0], [0, 0, -1]]), [5.0])
        obs = self.small_mc_var_bins._find_distance_classes(
            self.small_mc_var_bins.DistanceMatrices[1], 1)
        self.compare_multiple_level_array(obs, exp)

        # Multiple distance classes (even #).
        exp = (array([[-1, 0, 0], [0, -1, 1], [0, 1, -1]]), [3.5, 6.5])
        obs = self.small_mc_var_bins._find_distance_classes(
            self.small_mc_var_bins.DistanceMatrices[1], 2)
        self.compare_multiple_level_array(obs, exp)

        # Multiple distance classes (odd #).
        exp = (array([[-1, 0, 1], [0, -1, 2], [1, 2, -1]]),
               [2.0, 3.5, 6.5])
        obs = self.small_mc_var_bins._find_distance_classes(
            self.small_mc_var_bins.DistanceMatrices[1], 3)
        self.compare_multiple_level_array(obs, exp)

        # More classes than distances.
        exp = (array([[-1, 0, 1], [0, -1, 2], [1, 2, -1]]),
               [2.0, 3.5, 6.5, 8])
        obs = self.small_mc_var_bins._find_distance_classes(
            self.small_mc_var_bins.DistanceMatrices[1], 4)
        self.compare_multiple_level_array(obs, exp)

    def test_find_distance_classes_invalid_num_classes(self):
        """Test finding the distance classes for a bad number of classes."""
        self.assertRaises(ValueError, self.mc._find_distance_classes,
                          self.mc.DistanceMatrices[1], 0)
        self.assertRaises(ValueError, self.mc._find_distance_classes,
                          self.mc.DistanceMatrices[1], -1)

    def test_find_row_col_indices(self):
        """Test finds the row and col based on a flattened-list index."""
        obs = self.mc._find_row_col_indices(0)
        self.assertEqual(obs, (1, 0))

        obs = self.mc._find_row_col_indices(1)
        self.assertEqual(obs, (2, 0))

        obs = self.mc._find_row_col_indices(2)
        self.assertEqual(obs, (2, 1))

        obs = self.mc._find_row_col_indices(3)
        self.assertEqual(obs, (3, 0))

        obs = self.mc._find_row_col_indices(4)
        self.assertEqual(obs, (3, 1))

        obs = self.mc._find_row_col_indices(5)
        self.assertEqual(obs, (3, 2))

        obs = self.mc._find_row_col_indices(6)
        self.assertEqual(obs, (4, 0))

        self.assertRaises(IndexError, self.mc._find_row_col_indices, -1)

    def test_find_break_points(self):
        """Test finding equal-spaced breakpoints in a range."""
        exp = [-2.2204460492503131e-16, 1.0, 2.0, 3.0, 4.0, 5.0]
        obs = self.mc._find_break_points(0, 5, 5)
        assert_almost_equal(obs, exp)

        exp = [-2.0, -1.66666666667, -1.33333333333, -1.0]
        obs = self.mc._find_break_points(-2, -1, 3)
        assert_almost_equal(obs, exp)

        exp = [-1.0, -0.5, 0.0, 0.5, 1.0]
        obs = self.mc._find_break_points(-1, 1, 4)
        assert_almost_equal(obs, exp)

        exp = [-1.0, 1.0]
        obs = self.mc._find_break_points(-1, 1, 1)
        assert_almost_equal(obs, exp)

    def test_find_break_points_invalid_range(self):
        """Test finding breakpoints on an invalid range."""
        self.assertRaises(ValueError, self.mc._find_break_points, 1, 0, 5)
        self.assertRaises(ValueError, self.mc._find_break_points, 1, 1, 5)

    def test_find_break_points_invalid_num_classes(self):
        """Test finding breakpoints with an invalid number of classes."""
        self.assertRaises(ValueError, self.mc._find_break_points, 0, 1, 0)
        self.assertRaises(ValueError, self.mc._find_break_points, 0, 1, -1)

    def test_correct_p_values(self):
        """Test p-value correction for a small list of p-values."""
        exp = [0.003, 0.006, 0.003]
        obs = self.mc._correct_p_values([0.001, 0.002, 0.001])
        assert_almost_equal(obs, exp)

    def test_correct_p_values_all_None(self):
        """Test p-value correction for all None p-values."""
        exp = [None, None]
        obs = self.mc._correct_p_values([None, None])
        self.assertEqual(obs, exp)

    def test_correct_p_values_mixed(self):
        """Test p-value correction for mixture of None and valid p-values."""
        exp = [None, 0.008, 0.01, None]
        obs = self.mc._correct_p_values([None, 0.004, 0.005, None])
        self.assertEqual(obs, exp)

    def test_correct_p_values_no_change(self):
        """Test p-value correction where none is needed."""
        exp = [None, 0.008]
        obs = self.mc._correct_p_values([None, 0.008])
        self.assertEqual(obs, exp)
        exp = [0.007]
        obs = self.mc._correct_p_values([0.007])
        assert_almost_equal(obs, exp)

    def test_correct_p_values_large_correction(self):
        """Test p-value correction that exceeds 1.0."""
        exp = [1, None, 0.03, 0.03]
        obs = self.mc._correct_p_values([0.5, None, 0.01, 0.01])
        self.compare_multiple_level_array(obs, exp)

    def test_correct_p_values_empty(self):
        """Test p-value correction on empty list."""
        exp = []
        obs = self.mc._correct_p_values([])
        assert_almost_equal(obs, exp)

    def test_generate_correlogram(self):
        """Test creating a correlogram plot."""
        obs_fig = self.mc._generate_correlogram([0, 1, 2], [-0.9, 0, 0.9],
                                                [0.001, 0.1, 0.9])
        obs_ax = obs_fig.get_axes()[0]
        self.assertEqual(obs_ax.get_title(), "Mantel Correlogram")
        self.assertEqual(obs_ax.get_xlabel(), "Distance class index")
        self.assertEqual(obs_ax.get_ylabel(), "Mantel correlation statistic")
        assert_almost_equal(obs_ax.get_xticks(), [0., 0.5, 1., 1.5, 2.])
        assert_almost_equal(obs_ax.get_yticks(), [-1., -0.5, 0., 0.5, 1.])

    def test_generate_correlogram_empty(self):
        """Test creating a correlogram plot with no data."""
        obs_fig = self.mc._generate_correlogram([], [], [])
        obs_ax = obs_fig.get_axes()[0]
        self.assertEqual(obs_ax.get_title(), "Mantel Correlogram")
        self.assertEqual(obs_ax.get_xlabel(), "Distance class index")
        self.assertEqual(obs_ax.get_ylabel(), "Mantel correlation statistic")


class MantelTests(TestHelper):

    """Tests for the Mantel class."""

    def setUp(self):
        """Set up Mantel instances for use in tests."""
        super(MantelTests, self).setUp()

        # Create two small test distance matrices.
        sample_ids = ["S1", "S2", "S3"]
        m1 = array([[0, 1, 2], [1, 0, 3], [2, 3, 0]])
        m2 = array([[0, 2, 7], [2, 0, 6], [7, 6, 0]])
        m1_dm = DistanceMatrix(m1, sample_ids)
        m2_dm = DistanceMatrix(m2, sample_ids)

        self.small_mantel = Mantel(m1_dm, m2_dm, 'less')
        self.overview_mantel = Mantel(self.overview_dm, self.overview_dm,
                                      'greater')

    def test_DistanceMatrices_setter(self):
        """Test setting matrices using a valid number of distance matrices."""
        dms = [self.overview_dm, self.overview_dm]
        self.overview_mantel.DistanceMatrices = dms
        self.assertEqual(self.overview_mantel.DistanceMatrices, dms)

    def test_DistanceMatrices_setter_wrong_number(self):
        """Test setting an invalid number of distance matrices."""
        self.assertRaises(ValueError, setattr, self.overview_mantel,
                          'DistanceMatrices', [self.overview_dm])
        self.assertRaises(ValueError, setattr, self.overview_mantel,
                          'DistanceMatrices', [self.overview_dm, self.overview_dm,
                                               self.overview_dm])

    def test_DistanceMatrices_setter_too_small(self):
        """Test setting distance matrices that are too small."""
        self.assertRaises(ValueError, setattr, self.overview_mantel,
                          'DistanceMatrices', [self.single_ele_dm, self.single_ele_dm])

    def test_call_overview(self):
        """Runs mantel test on the overview dm when compared to itself.

        Expected R output:
            Mantel statistic r: 1
            Significance: 0.001

        Based on 999 permutations
        """
        expected_method_name = "Mantel"
        expected_p_value = 0.001
        expected_r_value = 1.0
        expected_perm_stats_len = 999
        expected_number_of_permutations = 999
        expected_tail_type = "greater"

        overview_mantel_output = self.overview_mantel(999)

        obs_method_name = overview_mantel_output['method_name']
        obs_num_permutations = overview_mantel_output['num_perms']
        obs_r_value = overview_mantel_output['r_value']
        obs_perm_stats_len = len(overview_mantel_output['perm_stats'])
        obs_tail_type = overview_mantel_output['tail_type']

        self.assertEqual(expected_method_name, obs_method_name)
        assert_almost_equal(expected_r_value, obs_r_value)
        assert_almost_equal(expected_perm_stats_len, obs_perm_stats_len)
        self.assertEqual(expected_number_of_permutations, obs_num_permutations)
        self.assertEqual(expected_tail_type, obs_tail_type)
        self.assertCorrectPValue(0, 0.006, self.overview_mantel, 999)

    def test_call_small(self):
        """Test one-sided mantel test (less) on small example dataset."""
        # This test output was verified by R (their mantel function does a
        # one-sided greater test, but I modified their output to do a one-sided
        # less test).
        results = self.small_mantel(999)

        self.assertEqual(results['method_name'], 'Mantel')
        self.assertEqual(results['num_perms'], 999)
        self.assertEqual(results['tail_type'], 'less')
        assert_almost_equal(results['r_value'], 0.755928946018)
        self.assertEqual(len(results['perm_stats']), 999)
        self.assertCorrectPValue(0.6, 1.0, self.small_mantel, 999)


class PartialMantelTests(TestHelper):

    """Tests for the PartialMantel class."""

    def setUp(self):
        """Set up PartialMantel instances for use in tests."""
        super(PartialMantelTests, self).setUp()

        # Test partial Mantel using the unifrac dm from the overview tutorial
        # as all three inputs (should be a small value).
        self.pm = PartialMantel(self.overview_dm, self.overview_dm,
                                self.overview_dm)

        # Just a small matrix that is easy to edit and observe.
        smpl_ids = ['s1', 's2', 's3']
        self.small_pm = PartialMantel(
            DistanceMatrix(array([[0, 1, 4], [1, 0, 3], [4, 3, 0]]), smpl_ids),
            DistanceMatrix(array([[0, 2, 5], [2, 0, 8], [5, 8, 0]]), smpl_ids),
            DistanceMatrix(array([[0, 9, 10], [9, 0, 2], [10, 2, 0]]),
                           smpl_ids))

        self.small_pm_diff = PartialMantel(
            DistanceMatrix(array([[0, 1, 4], [1, 0, 3], [4, 3, 0]]), smpl_ids),
            DistanceMatrix(array([[0, 20, 51], [20, 0, 888], [51, 888, 0]]),
                           smpl_ids),
            DistanceMatrix(array([[0, 9, 10], [9, 0, 2], [10, 2, 0]]),
                           smpl_ids))

        smpl_ids = ['s1', 's2', 's3', 's4', 's5']
        self.small_pm_diff2 = PartialMantel(
            DistanceMatrix(array([[0, 1, 2, 3, 1.4],
                                  [1, 0, 1.5, 1.6, 1.7],
                                  [2, 1.5, 0, 0.8, 1.9],
                                  [3, 1.6, 0.8, 0, 1.0],
                                  [1.4, 1.7, 1.9, 1.0, 0]]), smpl_ids),
            DistanceMatrix(array([[0, 1, 2, 3, 4.1],
                                  [1, 0, 5, 6, 7],
                                  [2, 5, 0, 8, 9],
                                  [3, 6, 8, 0, 10],
                                  [4.1, 7, 9, 10, 0]]), smpl_ids),
            DistanceMatrix(array([[0, 1, 2, 3, 4],
                                  [1, 0, 5, 6, 7],
                                  [2, 5, 0, 8, 9.1],
                                  [3, 6, 8, 0, 10],
                                  [4, 7, 9.1, 10, 0]]), smpl_ids))

    def test_DistanceMatrices_setter(self):
        """Test setting matrices using a valid number of distance matrices."""
        dms = [self.overview_dm, self.overview_dm, self.overview_dm]
        self.pm.DistanceMatrices = dms
        self.assertEqual(self.pm.DistanceMatrices, dms)

    def test_DistanceMatrices_setter_wrong_number(self):
        """Test setting an invalid number of distance matrices."""
        self.assertRaises(ValueError, setattr, self.pm,
                          'DistanceMatrices', [self.overview_dm])
        self.assertRaises(ValueError, setattr, self.pm,
                          'DistanceMatrices', [self.overview_dm, self.overview_dm])

    def test_DistanceMatrices_setter_too_small(self):
        """Test setting distance matrices that are too small."""
        self.assertRaises(ValueError, setattr, self.pm, 'DistanceMatrices',
                          [self.single_ele_dm, self.single_ele_dm, self.single_ele_dm])

    def test_call_small(self):
        """Test the running of partial Mantel analysis on small input."""
        obs = self.small_pm()
        exp_method_name = 'Partial Mantel'
        self.assertEqual(obs['method_name'], exp_method_name)

        exp_mantel_r = 0.99999999999999944
        assert_almost_equal(obs['mantel_r'], exp_mantel_r)
        # We're not testing that this p-value falls between a certain range
        # because this test has poor stability across platforms/numpy
        # configurations. Just make sure the p-value is between 0 and 1.
        self.assertTrue(0.0 <= obs['mantel_p'] <= 1.0)

        obs = self.small_pm_diff()
        exp_method_name = 'Partial Mantel'
        self.assertEqual(obs['method_name'], exp_method_name)

        exp_mantel_r = 0.99999999999999734
        assert_almost_equal(obs['mantel_r'], exp_mantel_r)
        self.assertCorrectPValue(0.25, 0.4, self.small_pm_diff,
                                 p_val_key='mantel_p')

        obs = self.small_pm_diff2()
        exp_method_name = 'Partial Mantel'
        self.assertEqual(obs['method_name'], exp_method_name)

        exp_mantel_r = -0.350624881409
        assert_almost_equal(obs['mantel_r'], exp_mantel_r)
        self.assertCorrectPValue(0.8, 1.0, self.small_pm_diff2,
                                 p_val_key='mantel_p')


class TopLevelTests(TestHelper):

    def setUp(self):
        pass

    def test_quantile(self):
        """checks for correct quantile statistic values"""

        # suffle the data to be sure, it is getting sorted
        sample_data = array(range(1, 11))
        shuffle(sample_data)

        # regular cases
        expected_output = [1.9, 2.8, 3.25, 5.5, 7.75, 7.93]
        list_of_quantiles = [0.1, 0.2, 0.25, 0.5, 0.75, 0.77]
        output = quantile(sample_data, list_of_quantiles)
        assert_almost_equal(expected_output, output)

        sample_data = array([42, 32, 24, 57, 15, 34, 83, 24, 60, 67, 55, 17,
                             83, 17, 80, 65, 14, 34, 39, 53])
        list_of_quantiles = [0.5]
        output = quantile(sample_data, list_of_quantiles)
        assert_almost_equal(output, median(sample_data))

        # quantiles must be between [0, 1]
        with self.assertRaises(AssertionError):
            output = quantile(sample_data, [0.1, 0.2, -0.1, 2, 0.3, 0.5])

        # quantiles must be a list or a numpy array
        with self.assertRaises(AssertionError):
            output = quantile(sample_data, 1)

        # the data must be a list or a numpy array
        with self.assertRaises(AssertionError):
            output = quantile(1, [0])

    def test__quantile(self):
        """checks for correct quantiles according to R. type 7 algorithm"""
        # regular cases
        sample_data = array(range(25, 42))
        assert_almost_equal(_quantile(sample_data, 0.5), median(sample_data))

        # sorted data is assumed for this function
        sample_data = sorted(
            array([0.17483293, 0.99891939, 0.81377467, 0.8137437,
                   0.51990174, 0.35521497, 0.98751461]))
        assert_almost_equal(_quantile(sample_data, 0.10), 0.283062154)


class PairedDifferenceTests(TestHelper):

    def setUp(self):
        self.personal_ids_to_state_values1 = \
            {'firmicutes-abundance':
             {'subject1': [0.45, 0.55],
              'subject2': [0.11, 0.52]},
             'bacteroidetes-abundance':
             {'subject1': [0.28, 0.21],
              'subject2': [0.11, 0.01]}
             }
        self.personal_ids_to_state_values2 = \
            {'firmicutes-abundance':
             {'subject1': [0.45, 0.55],
              'subject2': [0.11, None]},
             'bacteroidetes-abundance':
             {'subject1': [0.28, 0.21],
              'subject2': [0.11, 0.01]}
             }
        self.files_to_remove = []
        self.dirs_to_remove = []
        tmp_dir = get_qiime_temp_dir()
        self.test_out = mkdtemp(dir=tmp_dir,
                                prefix='qiime_paired_diff_tests_',
                                suffix='')
        self.dirs_to_remove.append(self.test_out)

    def tearDown(self):

        remove_files(self.files_to_remove)
        # remove directories last, so we don't get errors
        # trying to remove files which may be in the directories
        for d in self.dirs_to_remove:
            if exists(d):
                rmtree(d)

    def test_paired_difference_analyses(self):
        """paired_difference_analyses functions as expected
        """
        actual = paired_difference_analyses(
            self.personal_ids_to_state_values1,
            ['firmicutes-abundance',
             'bacteroidetes-abundance'],
            ['Pre', 'Post'],
            output_dir=self.test_out,
            ymin=0.0,
            ymax=1.0)
        self.assertTrue(exists(join(self.test_out,
                                    'paired_difference_comparisons.txt')))
        self.assertTrue(
            exists(join(self.test_out, 'firmicutes-abundance.pdf')))
        self.assertTrue(
            exists(join(self.test_out, 'bacteroidetes-abundance.pdf')))
        # three output paths returned
        self.assertEqual(len(actual[0]), 5)
        # expected t values returned, they should be less than (firmicutes) or
        # greater (bacteroidetes) than 2
        self.assertLess(abs(actual[1]['firmicutes-abundance'][4]), 2)
        self.assertLess(2, abs(actual[1]['bacteroidetes-abundance'][4]))

    def test_paired_difference_analyses_biom_output(self):
        """paired_difference_analyses generates correct biom tables
        """
        actual = paired_difference_analyses(
            self.personal_ids_to_state_values1,
            ['firmicutes-abundance',
             'bacteroidetes-abundance'],
            ['Pre', 'Post'],
            output_dir=self.test_out,
            ymin=0.0,
            ymax=1.0)
        biom_table_fp = join(self.test_out, 'differences.biom')
        self.assertTrue(exists(biom_table_fp))
        self.assertTrue(exists(join(self.test_out, 'differences_sids.txt')))
        with biom_open(biom_table_fp) as biom_file:
            table = Table.from_hdf5(biom_file)
        self.assertItemsEqual(table.ids(), ['subject1', 'subject2'])
        self.assertItemsEqual(table.ids(axis='observation'),
                              ['firmicutes-abundance', 'bacteroidetes-abundance'])
        assert_almost_equal(table
                              [(table.index('firmicutes-abundance',
                                            axis='observation'),
                               table.index('subject1', axis='sample'))],
                              0.1, 2)
        assert_almost_equal(table
                              [(table.index('bacteroidetes-abundance',
                                            axis='observation'),
                              table.index('subject1', axis='sample'))],
                              -0.07, 2)
        assert_almost_equal(table
                              [(table.index('firmicutes-abundance',
                                            axis='observation'),
                              table.index('subject2', axis='sample'))],
                              0.41, 2)
        assert_almost_equal(table
                              [(table.index('bacteroidetes-abundance',
                                            axis='observation'),
                              table.index('subject2', axis='sample'))],
                              -0.10, 2)

        # missing data should raise ValueError
        self.assertRaises(ValueError, paired_difference_analyses,
                          self.personal_ids_to_state_values2,
                          ['firmicutes-abundance',
                           'bacteroidetes-abundance'],
                          ['Pre', 'Post'],
                          output_dir=self.test_out,
                          ymin=0.0,
                          ymax=1.0)

    def test_paired_difference_analyses_wo_ymin_ymax(self):
        """paired_difference_analyses functions as expected w/o ymin/ymax
        """
        # runs successfully with ymin/ymax
        actual = paired_difference_analyses(
            self.personal_ids_to_state_values1,
            ['firmicutes-abundance',
             'bacteroidetes-abundance'],
            ['Pre', 'Post'],
            output_dir=self.test_out,
            ymin=None,
            ymax=None)
        self.assertTrue(exists(join(self.test_out,
                                    'paired_difference_comparisons.txt')))
        self.assertTrue(
            exists(join(self.test_out, 'firmicutes-abundance.pdf')))
        self.assertTrue(
            exists(join(self.test_out, 'bacteroidetes-abundance.pdf')))
        # three output paths returned
        self.assertEqual(len(actual[0]), 5)
        # expected t values returned, they should be less than (firmicutes) or
        # greater (bacteroidetes) than 2
        self.assertLess(0, actual[1]['firmicutes-abundance'][4])
        self.assertLess(actual[1]['bacteroidetes-abundance'][4], 0)

    def test_paired_difference_analyses_analysis_cat_subset(self):
        """paired_difference_analyses fns w a subset of analysis categories
        """
        actual = paired_difference_analyses(
            self.personal_ids_to_state_values1,
            ['firmicutes-abundance'],
            ['Pre', 'Post'],
            output_dir=self.test_out,
            ymin=0.0,
            ymax=1.0)
        self.assertTrue(exists(join(self.test_out,
                                    'paired_difference_comparisons.txt')))
        self.assertTrue(
            exists(join(self.test_out, 'firmicutes-abundance.pdf')))
        self.assertFalse(
            exists(join(self.test_out, 'bacteroidetes-abundance.pdf')))
        # three output paths returned
        self.assertEqual(len(actual[0]), 4)
        # expected t values returned
        assert_almost_equal(actual[1]['firmicutes-abundance'][4], 1.645, 3)


class TestsHelper(TestCase):

    """Class with utility methods useful for other tests."""

    # How many times a p-value should be tested to fall in a given range
    # before failing the test.
    p_val_tests = 20

    def assertCorrectPValue(self, exp_min, exp_max, fn, args=None,
                            kwargs=None, p_val_idx=0):
        """Tests that the stochastic p-value falls in the specified range.

        Performs the test self.p_val_tests times and fails if the observed
        p-value does not fall into the specified range at least once. Each
        p-value is also tested that it falls in the range 0.0 to 1.0.

        This method assumes that fn is callable, and will unpack and pass args
        and kwargs to fn if they are provided. It also assumes that fn returns
        a single value (the p-value to be tested) or a tuple of results (any
        length greater than or equal to 1), with the p-value at position
        p_val_idx.

        This is primarily used for testing the Mantel and correlation_test
        functions.
        """

        found_match = False
        for i in range(self.p_val_tests):
            if args is not None and kwargs is not None:
                obs = fn(*args, **kwargs)
            elif args is not None:
                obs = fn(*args)
            elif kwargs is not None:
                obs = fn(**kwargs)
            else:
                obs = fn()

            try:
                p_val = float(obs)
            except TypeError:
                p_val = obs[p_val_idx]
            self.assertTrue(0.0 <= p_val <= 1.0)
            if p_val >= exp_min and p_val <= exp_max:
                found_match = True
                break
        self.assertTrue(found_match)


class TestsTests(TestCase):

    """Tests miscellaneous functions."""

    def test_tail(self):
        """tail should return x/2 if test is true; 1-(x/2) otherwise"""
        assert_allclose(tail(0.25, 'a' == 'a'), 0.25 / 2)
        assert_allclose(tail(0.25, 'a' != 'a'), 1 - (0.25 / 2))

    def test_fisher(self):
        """fisher results should match p 795 Sokal and Rohlf"""
        assert_allclose(fisher([0.073, 0.086, 0.10, 0.080, 0.060]),
                        0.0045957946540917905, atol=10e-7)

    def test_permute_2d(self):
        """permute_2d permutes rows and cols of a matrix."""
        a = reshape(arange(9), (3, 3))
        assert_allclose(permute_2d(a, [0, 1, 2]), a)
        assert_allclose(permute_2d(a, [2, 1, 0]),
                        array([[8, 7, 6], [5, 4, 3],
                               [2, 1, 0]]))
        assert_allclose(permute_2d(a, [1, 2, 0]),
                        array([[4, 5, 3], [7, 8, 6],
                               [1, 2, 0]]))


class GTests(TestCase):

    """Tests implementation of the G tests for fit and independence."""

    def test_G_2_by_2_2tailed_equal(self):
        """G_2_by_2 should return 0 if all cell counts are equal"""
        assert_allclose(0, G_2_by_2(1, 1, 1, 1, False, False)[0])
        assert_allclose(0, G_2_by_2(100, 100, 100, 100, False,
                                    False)[0])
        assert_allclose(0, G_2_by_2(100, 100, 100, 100, True,
                                    False)[0])

    def test_G_2_by_2_bad_data(self):
        """G_2_by_2 should raise ValueError if any counts are negative"""
        self.assertRaises(ValueError, G_2_by_2, 1, -1, 1, 1)

    def test_G_2_by_2_2tailed_examples(self):
        """G_2_by_2 values should match examples in Sokal & Rohlf"""
        # example from p 731, Sokal and Rohlf (1995)
        # without correction
        assert_allclose(G_2_by_2(12, 22, 16, 50, False, False)[0],
                        1.33249, 0.0001)
        assert_allclose(G_2_by_2(12, 22, 16, 50, False, False)[1],
                        0.24836, 0.0001)
        # with correction
        assert_allclose(G_2_by_2(12, 22, 16, 50, True, False)[0],
                        1.30277, 0.0001)
        assert_allclose(G_2_by_2(12, 22, 16, 50, True, False)[1],
                        0.25371, 0.0001)

    def test_G_2_by_2_1tailed_examples(self):
        """G_2_by_2 values should match values from codon_binding program"""
        # first up...the famous arginine case
        assert_allclose(G_2_by_2(36, 16, 38, 106), (29.111609, 0),
                        atol=10e-7)
        # then some other miscellaneous positive and negative values
        assert_allclose(
            G_2_by_2(0, 52, 12, 132), (-7.259930, 0.996474), atol=10e-7)
        assert_allclose(
            G_2_by_2(5, 47, 14, 130), (-0.000481, 0.508751), atol=10e-7)
        assert_allclose(
            G_2_by_2(5, 47, 36, 108), (-6.065167, 0.993106), atol=10e-7)

    def test_g_fit(self):
        """Test G fit is correct with and without Williams correction."""
        # test with williams correction
        data = [array(i) for i in [63, 31, 28, 12, 39, 16, 40, 12]]
        exp_G = 69.030858949133162 / 1.00622406639
        exp_p = 2.8277381487281706e-12
        obs_G, obs_p = g_fit(data, williams=True)
        assert_allclose(obs_G, exp_G)
        assert_allclose(obs_p, exp_p, atol=1e-7)
        # test with hand computed example and williams correction
        data = [array([75, 65, 48]), array([200]), array([10, 250, 13,
                                                          85])]
        exp_G = 85.90859811005285 / 1.0018930430667
        exp_p = 2.4012235241479195e-19
        obs_G, obs_p = g_fit(data, williams=True)
        assert_allclose(obs_G, exp_G)
        assert_allclose(obs_p, exp_p, atol=1e-7)
        # test without williams correction on another hand computed example
        data = [array([10, 12, 15, 7]), array([15, 12, 17, 18]),
                array([6, 9, 13])]
        exp_G = 1.6610421781232
        exp_p = 0.43582212499949591
        obs_G, obs_p = g_fit(data, williams=False)
        assert_allclose(obs_G, exp_G)
        assert_allclose(obs_p, exp_p, atol=1e-7)

    def test_williams_correction(self):
        """Test that the Williams correction is correctly computed."""
        n = 100
        a = 10
        G = 10.5783
        exp = 10.387855973813421
        assert_allclose(williams_correction(n, a, G), exp,
                        rtol=1e-5)
        # test with an example from Sokal and Rohlf pg 699
        n = 241
        a = 8
        G = 8.82396
        exp = 8.76938
        assert_allclose(williams_correction(n, a, G), exp,
                        rtol=1e-5)

    def test_safe_sum_p_log_p(self):
        """safe_sum_p_log_p should ignore zero elements, not raise error"""
        m = array([2, 4, 0, 8])
        self.assertEqual(safe_sum_p_log_p(m, 2), 2 * 1 + 4 * 2 + 8 * 3)


class StatTests(TestsHelper):

    """Tests that the t and z tests are implemented correctly"""

    def setUp(self):
        super(StatTests, self).setUp()

        self.x = [7.33, 7.49, 7.27, 7.93, 7.56, 7.81, 7.46, 6.94, 7.49, 7.44,
                  7.95, 7.47, 7.04, 7.10, 7.64]
        self.y = [7.53, 7.70, 7.46, 8.21, 7.81, 8.01, 7.72, 7.13, 7.68, 7.66,
                  8.11, 7.66, 7.20, 7.25, 7.79]

    def test_t_paired_2tailed(self):
        """t_paired should match values from Sokal & Rohlf p 353"""
        x, y = self.x, self.y
        # check value of t and the probability for 2-tailed
        assert_allclose(t_paired(y, x)[0], 19.7203, 1e-4)
        assert_allclose(t_paired(y, x)[1], 1.301439e-11, 1e-4)

    def test_t_paired_no_variance(self):
        """t_paired should return None if lists are invariant"""
        x = [1, 1, 1]
        y = [0, 0, 0]
        assert_allclose(t_paired(x, x), (nan, nan))
        assert_allclose(t_paired(x, y), (nan, nan))

    def test_t_paired_1tailed(self):
        """t_paired should match pre-calculated 1-tailed values"""
        x, y = self.x, self.y
        # check probability for 1-tailed low and high
        assert_allclose(
            t_paired(y, x, "low")[1], 1 - (1.301439e-11 / 2), 1e-4)
        assert_allclose(
            t_paired(x, y, "high")[1], 1 - (1.301439e-11 / 2), 1e-4)
        assert_allclose(
            t_paired(y, x, "high")[1], 1.301439e-11 / 2, 1e-4)
        assert_allclose(
            t_paired(x, y, "low")[1], 1.301439e-11 / 2, 1e-4)

    def test_t_paired_specific_difference(self):
        """t_paired should allow a specific difference to be passed"""
        x, y = self.x, self.y
        # difference is 0.2, so test should be non-significant if 0.2 passed
        self.assertFalse(t_paired(y, x, exp_diff=0.2)[0] > 1e-10)
        # same, except that reversing list order reverses sign of difference
        self.assertFalse(t_paired(x, y, exp_diff=-0.2)[0] > 1e-10)
        # check that there's no significant difference from the true mean
        assert_allclose(
            t_paired(y, x, exp_diff=0.2)[1], 1, 1e-4)

    def test_t_paired_bad_data(self):
        """t_paired should raise ValueError on lists of different lengths"""
        self.assertRaises(ValueError, t_paired, self.y, [1, 2, 3])

    def test_t_two_sample(self):
        """t_two_sample should match example on p.225 of Sokal and Rohlf"""
        I = array([7.2, 7.1, 9.1, 7.2, 7.3, 7.2, 7.5])
        II = array([8.8, 7.5, 7.7, 7.6, 7.4, 6.7, 7.2])
        assert_allclose(t_two_sample(I, II, 'two-sided'),
                        (-0.1184, 0.45385 * 2),
                        atol=10e-3)

    def test_t_two_sample_no_variance(self):
        """t_two_sample should properly handle lists that are invariant"""
        # By default should return (None, None) to mimic R's t.test.
        x = array([1, 1., 1])
        y = array([0, 0, 0.0])
        self.assertEqual(t_two_sample(x, x), (nan, nan))
        self.assertEqual(t_two_sample(x, y), (nan, nan))

        # Should still receive (nan, nan) if the lists have no variance and
        # have the same single value.
        self.assertEqual(t_two_sample(x, x), (nan, nan))
        self.assertEqual(t_two_sample(x, [1, 1]), (nan, nan))

    def test_t_one_sample(self):
        """t_one_sample results should match those from R"""
        x = array(range(-5, 5))
        y = array(range(-1, 10))
        assert_allclose(t_one_sample(x), (-0.5222, 0.6141), atol=10e-3)
        assert_allclose(t_one_sample(y), (4, 0.002518), atol=10e-3)
        # do some one-tailed tests as well
        assert_allclose(t_one_sample(y, tails='low'), (4, 0.9987), atol=10e-3)
        assert_allclose(
            t_one_sample(y, tails='high'), (4, 0.001259), atol=10e-3)

    def test_t_two_sample_switch(self):
        """t_two_sample should call t_one_observation if 1 item in sample."""
        sample = array([4.02, 3.88, 3.34, 3.87, 3.18])
        x = array([3.02])
        assert_allclose(t_two_sample(x, sample), (-1.5637254, 0.1929248))
        assert_allclose(t_two_sample(sample, x), (-1.5637254, 0.1929248))

        # can't do the test if both samples have single item
        assert_allclose(t_two_sample(x, x), (nan, nan))

        # Test special case if t=0.
        assert_allclose(t_two_sample([2], [1, 2, 3]), (0.0, 1.0))
        assert_allclose(t_two_sample([1, 2, 3], [2]), (0.0, 1.0))

    def test_t_one_observation(self):
        """t_one_observation should match p. 228 of Sokal and Rohlf"""
        sample = array([4.02, 3.88, 3.34, 3.87, 3.18])
        x = 3.02
        # note that this differs after the 3rd decimal place from what's in
        # the book, because Sokal and Rohlf round their intermediate steps...
        assert_allclose(t_one_observation(x, sample), (-1.5637254, 0.1929248))

    def test_t_one_observation_no_variance(self):
        """t_one_observation should correctly handle an invariant list."""
        sample = array([1.0, 1.0, 1.0])

        assert_allclose(t_one_observation(1, sample), (nan, nan))
        assert_allclose(t_one_observation(2, sample, exp_diff=3), (nan, nan))
        assert_allclose(t_one_observation(2, sample, tails='low'), (nan, nan))

    def test_mc_t_two_sample(self):
        """Test gives correct results with valid input data."""
        # Verified against R's t.test() and Deducer::perm.t.test().

        # With numpy array as input.
        exp = (-0.11858541225631833, 0.90756579317867436)
        I = array([7.2, 7.1, 9.1, 7.2, 7.3, 7.2, 7.5])
        II = array([8.8, 7.5, 7.7, 7.6, 7.4, 6.7, 7.2])
        obs = mc_t_two_sample(I, II)
        assert_allclose(obs[:2], exp)
        self.assertEqual(len(obs[2]), 999)
        self.assertCorrectPValue(0.8, 0.9, mc_t_two_sample, [I, II],
                                 p_val_idx=3)

        # With python list as input.
        exp = (-0.11858541225631833, 0.90756579317867436)
        I = [7.2, 7.1, 9.1, 7.2, 7.3, 7.2, 7.5]
        II = [8.8, 7.5, 7.7, 7.6, 7.4, 6.7, 7.2]
        obs = mc_t_two_sample(I, II)
        assert_allclose(obs[:2], exp)
        self.assertEqual(len(obs[2]), 999)
        self.assertCorrectPValue(0.8, 0.9, mc_t_two_sample, [I, II],
                                 p_val_idx=3)

        exp = (-0.11858541225631833, 0.45378289658933718)
        obs = mc_t_two_sample(I, II, tails='low')
        assert_allclose(obs[:2], exp)
        self.assertEqual(len(obs[2]), 999)
        self.assertCorrectPValue(0.4, 0.47, mc_t_two_sample, [I, II],
                                 {'tails': 'low'}, p_val_idx=3)

        exp = (-0.11858541225631833, 0.54621710341066287)
        obs = mc_t_two_sample(I, II, tails='high', permutations=99)
        assert_allclose(obs[:2], exp)
        self.assertEqual(len(obs[2]), 99)
        self.assertCorrectPValue(0.4, 0.62, mc_t_two_sample, [I, II],
                                 {'tails': 'high', 'permutations': 99},
                                 p_val_idx=3)

        exp = (-2.8855783649036986, 0.99315596652421401)
        obs = mc_t_two_sample(I, II, tails='high',
                              permutations=99, exp_diff=1)
        assert_allclose(obs[:2], exp)
        self.assertEqual(len(obs[2]), 99)
        self.assertCorrectPValue(0.55, 0.99, mc_t_two_sample, [I, II],
                                 {'tails': 'high', 'permutations': 99,
                                  'exp_diff': 1}, p_val_idx=3)

    def test_mc_t_two_sample_unbalanced_obs(self):
        """Test gives correct results with unequal number of obs per sample."""
        # Verified against R's t.test() and Deducer::perm.t.test().
        exp = (-0.10302479888889175, 0.91979753020527177)
        I = array([7.2, 7.1, 9.1, 7.2, 7.3, 7.2])
        II = array([8.8, 7.5, 7.7, 7.6, 7.4, 6.7, 7.2])
        obs = mc_t_two_sample(I, II)
        assert_allclose(obs[:2], exp)
        self.assertEqual(len(obs[2]), 999)
        self.assertCorrectPValue(0.8, 0.9, mc_t_two_sample, [I, II],
                                 p_val_idx=3)

    def test_mc_t_two_sample_single_obs_sample(self):
        """Test works correctly with one sample having a single observation."""
        sample = array([4.02, 3.88, 3.34, 3.87, 3.18])
        x = array([3.02])
        exp = (-1.5637254, 0.1929248)
        obs = mc_t_two_sample(x, sample)
        assert_allclose(obs[:2], exp, atol=1e-6)
        assert_allclose(len(obs[2]), 999)
        self.assertTrue(0.0 <= obs[3] <= 1.0)

        # Test the case where we can have no variance in the permuted lists.
        x = array([1, 1, 2])
        y = array([1])
        exp = (-0.5, 0.666666666667)
        obs = mc_t_two_sample(x, y)
        assert_allclose(obs[:2], exp)
        assert_allclose(len(obs[2]), 999)
        self.assertTrue(0.0 <= obs[3] <= 1.0)

    def test_mc_t_two_sample_no_perms(self):
        """Test gives empty permutation results if no perms are given."""
        exp = (-0.11858541225631833, 0.90756579317867436, [], nan)
        I = array([7.2, 7.1, 9.1, 7.2, 7.3, 7.2, 7.5])
        II = array([8.8, 7.5, 7.7, 7.6, 7.4, 6.7, 7.2])
        obs = mc_t_two_sample(I, II, permutations=0)
        assert_allclose(obs[0], exp[0])
        assert_allclose(obs[1], exp[1])
        self.assertEqual(obs[2], exp[2])
        assert_allclose(obs[3], exp[3])

    def test_mc_t_two_sample_no_mc(self):
        """Test no MC stats if initial t-test is bad."""
        x = array([1, 1, 1])
        y = array([0, 0, 0])
        self.assertEqual(mc_t_two_sample(x, y), (nan, nan, [], nan))

    def test_mc_t_two_sample_no_variance(self):
        """Test input with no variance. Should match Deducer::perm.t.test."""
        x = array([1, 1, 1])
        y = array([2, 2, 2])

        exp = (nan, nan)
        obs = mc_t_two_sample(x, y, permutations=1000)
        self.assertEqual(obs[:2], exp)

    def test_mc_t_two_sample_no_permuted_variance(self):
        """Test with chance of getting no variance with some perms."""
        # Verified against R's t.test() and Deducer::perm.t.test().
        x = array([1, 1, 2])
        y = array([2, 2, 1])

        exp = (-0.70710678118654791, 0.51851851851851838)
        obs = mc_t_two_sample(x, y, permutations=1000)

        assert_allclose(obs[:2], exp)
        self.assertEqual(len(obs[2]), 1000)
        self.assertCorrectPValue(0.90, 1.0, mc_t_two_sample, [x, y],
                                 {'permutations': 1000}, p_val_idx=3)

    def test_mc_t_two_sample_invalid_input(self):
        """Test fails on various invalid input."""
        # self.assertRaises(ValueError, mc_t_two_sample, [1, 2, 3],
        #                   [4., 5., 4.], tails='foo')
        # self.assertRaises(ValueError, mc_t_two_sample, [1, 2, 3],
        #                   [4., 5., 4.], permutations=-1)
        self.assertRaises(ValueError, mc_t_two_sample, [1], [4.])
        self.assertRaises(ValueError, mc_t_two_sample, [1, 2], [])

    def test_permute_observations(self):
        """Test works correctly on small input dataset."""
        I = [10, 20., 1]
        II = [2, 4, 5, 7]
        obs = _permute_observations(I, II, 1)
        self.assertEqual(len(obs[0]), 1)
        self.assertEqual(len(obs[1]), 1)
        self.assertEqual(len(obs[0][0]), len(I))
        self.assertEqual(len(obs[1][0]), len(II))
        assert_allclose(sorted(concatenate((obs[0][0],
                                            obs[1][0]))),
                        sorted(I + II))

    def test_tail(self):
        """tail should return prob/2 if test is true, or 1-(prob/2) if false
        """
        assert_allclose(tail(0.25, True), 0.125)
        assert_allclose(tail(0.25, False), 0.875)
        assert_allclose(tail(1, True), 0.5)
        assert_allclose(tail(1, False), 0.5)
        assert_allclose(tail(0, True), 0)
        assert_allclose(tail(0, False), 1)


class CorrelationTests(TestsHelper):

    """Tests of correlation coefficients and Mantel test."""

    def setUp(self):
        """Sets up variables used in the tests."""
        super(CorrelationTests, self).setUp()

        # For testing spearman and correlation_test using method='spearman'.
        # Taken from the Spearman wikipedia article. Also used for testing
        # Pearson (verified with R).
        self.data1 = [106, 86, 100, 101, 99, 103, 97, 113, 112, 110]
        self.data2 = [7, 0, 27, 50, 28, 29, 20, 12, 6, 17]

        # For testing spearman.
        self.a = [1, 2, 4, 3, 1, 6, 7, 8, 10, 4]
        self.b = [2, 10, 20, 1, 3, 7, 5, 11, 6, 13]
        self.c = [7, 1, 20, 13, 3, 57, 5, 121, 2, 9]
        self.r = (1.7, 10, 20, 1.7, 3, 7, 5, 11, 6.5, 13)
        self.x = (1, 2, 4, 3, 1, 6, 7, 8, 10, 4, 100, 2, 3, 77)

        # Ranked copies for testing spearman.
        self.b_ranked = [2, 7, 10, 1, 3, 6, 4, 8, 5, 9]
        self.c_ranked = [5, 1, 8, 7, 3, 9, 4, 10, 2, 6]

        # silence the warnings that will tests for correlation_test
        filterwarnings('ignore', category=RuntimeWarning)

    def test_mantel(self):
        """mantel should be significant for same matrix, not for random"""
        a = reshape(arange(25), (5, 5))
        a = tril(a) + tril(a).T
        fill_diagonal(a, 0)
        b = a.copy()
        # closely related -- should be significant
        self.assertCorrectPValue(0.0, 0.049, mantel, (a, b, 1000))

        c = reshape(ones(25), (5, 5))
        c[0, 1] = 3.0
        c[1, 0] = 3.0
        fill_diagonal(c, 0)
        # not related -- should not be significant
        self.assertCorrectPValue(0.06, 1.0, mantel, (a, c, 1000))

    def test_mantel_test_one_sided_greater(self):
        """Test one-sided mantel test (greater)."""
        # This test output was verified by R (their mantel function does a
        # one-sided greater test).
        m1 = array([[0, 1, 2], [1, 0, 3], [2, 3, 0]])
        m2 = array([[0, 2, 7], [2, 0, 6], [7, 6, 0]])
        p, stat, perms = mantel_t(m1, m1, 999, alt='greater')
        assert_allclose(stat, 1.0)
        self.assertEqual(len(perms), 999)

        self.assertCorrectPValue(0.09, 0.25, mantel_t, (m1, m1, 999),
                                 {'alt': 'greater'})

        p, stat, perms = mantel_t(m1, m2, 999, alt='greater')
        assert_allclose(stat, 0.755928946018)
        self.assertEqual(len(perms), 999)
        self.assertCorrectPValue(0.2, 0.5, mantel_t, (m1, m2, 999),
                                 {'alt': 'greater'})

    def test_mantel_test_one_sided_less(self):
        """Test one-sided mantel test (less)."""
        # This test output was verified by R (their mantel function does a
        # one-sided greater test, but I modified their output to do a one-sided
        # less test).
        m1 = array([[0, 1, 2], [1, 0, 3], [2, 3, 0]])
        m2 = array([[0, 2, 7], [2, 0, 6], [7, 6, 0]])
        m3 = array([[0, 0.5, 0.25], [0.5, 0, 0.1], [0.25, 0.1, 0]])
        p, stat, perms = mantel_t(m1, m1, 999, alt='less')
        assert_allclose(p, 1.0)
        assert_allclose(stat, 1.0)
        self.assertEqual(len(perms), 999)

        p, stat, perms = mantel_t(m1, m2, 999, alt='less')
        assert_allclose(stat, 0.755928946018)
        self.assertEqual(len(perms), 999)
        self.assertCorrectPValue(0.6, 1.0, mantel_t, (m1, m2, 999),
                                 {'alt': 'less'})

        p, stat, perms = mantel_t(m1, m3, 999, alt='less')
        assert_allclose(stat, -0.989743318611)
        self.assertEqual(len(perms), 999)
        self.assertCorrectPValue(0.1, 0.25, mantel_t, (m1, m3, 999),
                                 {'alt': 'less'})

    def test_mantel_test_two_sided(self):
        """Test two-sided mantel test."""
        # This test output was verified by R (their mantel function does a
        # one-sided greater test, but I modified their output to do a two-sided
        # test).
        m1 = array([[0, 1, 2], [1, 0, 3], [2, 3, 0]])
        m2 = array([[0, 2, 7], [2, 0, 6], [7, 6, 0]])
        m3 = array([[0, 0.5, 0.25], [0.5, 0, 0.1], [0.25, 0.1, 0]])
        p, stat, perms = mantel_t(m1, m1, 999, alt='two sided')
        assert_allclose(stat, 1.0)
        self.assertEqual(len(perms), 999)
        self.assertCorrectPValue(0.20, 0.45, mantel_t, (m1, m1, 999),
                                 {'alt': 'two sided'})

        p, stat, perms = mantel_t(m1, m2, 999, alt='two sided')
        assert_allclose(stat, 0.755928946018)
        self.assertEqual(len(perms), 999)
        self.assertCorrectPValue(0.6, 0.75, mantel_t, (m1, m2, 999),
                                 {'alt': 'two sided'})

        p, stat, perms = mantel_t(m1, m3, 999, alt='two sided')
        assert_allclose(stat, -0.989743318611)
        self.assertEqual(len(perms), 999)
        self.assertCorrectPValue(0.2, 0.45, mantel_t, (m1, m3, 999),
                                 {'alt': 'two sided'})

    def test_mantel_test_invalid_distance_matrix(self):
        """Test mantel test with invalid distance matrix."""
        # Single asymmetric, non-hollow distance matrix.
        self.assertRaises(ValueError, mantel_t, array([[1, 2], [3, 4]]),
                          array([[0, 0], [0, 0]]), 999)

        # Two asymmetric distance matrices.
        self.assertRaises(ValueError, mantel_t, array([[0, 2], [3, 0]]),
                          array([[0, 1], [0, 0]]), 999)

    def test_mantel_test_invalid_input(self):
        """Test mantel test with invalid input."""
        self.assertRaises(ValueError, mantel_t, array([[1]]),
                          array([[1]]), 999, alt='foo')
        self.assertRaises(ValueError, mantel_t, array([[1]]),
                          array([[1, 2], [3, 4]]), 999)
        self.assertRaises(ValueError, mantel_t, array([[1]]),
                          array([[1]]), 0)
        self.assertRaises(ValueError, mantel_t, array([[1]]),
                          array([[1]]), -1)

    def test_is_symmetric_and_hollow(self):
        """Should correctly test for symmetry and hollowness of dist mats."""
        self.assertTrue(is_symmetric_and_hollow(array([[0, 1], [1, 0]])))
        self.assertTrue(is_symmetric_and_hollow(matrix([[0, 1], [1, 0]])))
        self.assertTrue(is_symmetric_and_hollow(matrix([[0.0, 0],
                                                        [0.0, 0]])))
        self.assertTrue(not is_symmetric_and_hollow(
            array([[0.001, 1], [1, 0]])))
        self.assertTrue(not is_symmetric_and_hollow(
            array([[0, 1.1], [1, 0]])))
        self.assertTrue(not is_symmetric_and_hollow(
            array([[0.5, 1.1], [1, 0]])))

    def test_flatten_lower_triangle(self):
        """Test flattening various dms' lower triangulars."""
        self.assertEqual(_flatten_lower_triangle(
            array([[8]])), [])
        self.assertEqual(_flatten_lower_triangle(
            array([[1, 2], [3, 4]])), [3])
        self.assertEqual(_flatten_lower_triangle(
            array([[1, 2, 3], [4, 5, 6], [7, 8, 9]])), [4, 7, 8])

    def test_pearson(self):
        """Test pearson correlation method on valid data."""
        # This test output was verified by R.
        assert_allclose(pearson([1, 2], [1, 2]), 1.0)
        assert_allclose(pearson([1, 2, 3], [1, 2, 3]), 1.0)
        assert_allclose(pearson([1, 2, 3], [1, 2, 4]), 0.9819805)

    def test_pearson_invalid_input(self):
        """Test running pearson on bad input."""
        self.assertRaises(ValueError, pearson, [1.4, 2.5], [5.6, 8.8, 9.0])
        self.assertRaises(ValueError, pearson, [1.4], [5.6])

    def test_spearman(self):
        """Test the spearman function with valid input."""
        # One vector has no ties.
        exp = 0.3719581
        obs = spearman(self.a, self.b)
        assert_allclose(obs, exp)

        # Both vectors have no ties.
        exp = 0.2969697
        obs = spearman(self.b, self.c)
        assert_allclose(obs, exp)

        # Both vectors have ties.
        exp = 0.388381
        obs = spearman(self.a, self.r)
        assert_allclose(obs, exp)

        exp = -0.17575757575757578
        obs = spearman(self.data1, self.data2)
        assert_allclose(obs, exp)

    def test_spearman_no_variation(self):
        """Test the spearman function with a vector having no variation."""
        exp = nan
        obs = spearman([1, 1, 1], [1, 2, 3])
        assert_allclose(obs, exp)

    def test_spearman_ranked(self):
        """Test the spearman function with a vector that is already ranked."""
        exp = 0.2969697
        obs = spearman(self.b_ranked, self.c_ranked)
        assert_allclose(obs, exp)

    def test_spearman_one_obs(self):
        """Test running spearman on a single observation."""
        self.assertRaises(ValueError, spearman, [1.0], [5.0])

    def test_spearman_invalid_input(self):
        """Test the spearman function with invalid input."""
        self.assertRaises(ValueError, spearman, [], [])
        self.assertRaises(ValueError, spearman, self.a, [])

    def test_correlation_test_pearson(self):
        """Test correlation_t using pearson on valid input."""
        # These results were verified with R.
        # Test with non-default confidence level and permutations.
        obs = correlation_t(self.data1, self.data2, method='pearson',
                            confidence_level=0.90, permutations=990)
        assert_allclose(obs[:2], (-0.03760147,
                                  0.91786297277172868), atol=10e-7)
        self.assertEqual(len(obs[2]), 990)
        for r in obs[2]:
            self.assertTrue(r >= -1.0 and r <= 1.0)
        self.assertCorrectPValue(0.9, 0.93, correlation_t,
                                 (self.data1, self.data2),
                                 {'method': 'pearson',
                                  'confidence_level': 0.90,
                                  'permutations': 990},
                                 p_val_idx=3)
        assert_allclose(obs[4], (-0.5779077, 0.5256224))

        # Test with non-default tail type.
        obs = correlation_t(self.data1, self.data2, method='pearson',
                            confidence_level=0.90, permutations=990,
                            tails='low')
        assert_allclose(obs[:2], (-0.03760147,
                                  0.45893148638586434), atol=10e-7)
        self.assertEqual(len(obs[2]), 990)
        for r in obs[2]:
            self.assertTrue(r >= -1.0 and r <= 1.0)
        self.assertCorrectPValue(0.41, 0.46, correlation_t,
                                 (self.data1, self.data2),
                                 {'method': 'pearson',
                                  'confidence_level': 0.90,
                                  'permutations': 990,
                                  'tails': 'low'},
                                 p_val_idx=3)
        assert_allclose(obs[4], (-0.5779077, 0.5256224))

    def test_correlation_test_spearman(self):
        """Test correlation_t using spearman on valid input."""
        # This example taken from Wikipedia page:
        # http://en.wikipedia.org/wiki/Spearman's_rank_correlation_coefficient
        obs = correlation_t(self.data1, self.data2, method='spearman',
                            tails='high')
        assert_allclose(obs[:2], (-0.17575757575757578,
                                  0.686405827612))
        self.assertEqual(len(obs[2]), 999)
        for rho in obs[2]:
            self.assertTrue(rho >= -1.0 and rho <= 1.0)
        self.assertCorrectPValue(0.67, 0.7, correlation_t,
                                 (self.data1, self.data2),
                                 {'method': 'spearman',
                                  'tails': 'high'},
                                 p_val_idx=3)
        assert_allclose(obs[4], (-0.7251388558041697,
                                 0.51034422964834503))

        # The p-value is off because the example uses a one-tailed test, while
        # we use a two-tailed test. Someone confirms the answer that we get
        # here for a two-tailed test:
        # http://stats.stackexchange.com/questions/22816/calculating-p-value-
        #     for-spearmans-rank-correlation-coefficient-example-on-wikip
        obs = correlation_t(self.data1, self.data2, method='spearman',
                            tails='two-sided')
        assert_allclose(obs[:2], (-0.17575757575757578,
                                  0.62718834477648433))
        self.assertEqual(len(obs[2]), 999)
        for rho in obs[2]:
            self.assertTrue(rho >= -1.0 and rho <= 1.0)
        self.assertCorrectPValue(0.60, 0.64, correlation_t,
                                 (self.data1, self.data2),
                                 {'method': 'spearman', 'tails': 'two-sided'},
                                 p_val_idx=3)
        assert_allclose(obs[4], (-0.7251388558041697,
                                 0.51034422964834503))

    def test_correlation_test_invalid_input(self):
        """Test correlation_t using invalid input."""
        self.assertRaises(ValueError, correlation_t, self.data1, self.data2,
                          method='foo')
        # self.assertRaises(ValueError, correlation_t, self.data1, self.data2,
        #                   tails='foo')
        self.assertRaises(ValueError, correlation_t, self.data1, self.data2,
                          permutations=-1)
        self.assertRaises(ValueError, correlation_t, self.data1, self.data2,
                          confidence_level=-1)
        self.assertRaises(ValueError, correlation_t, self.data1, self.data2,
                          confidence_level=1.1)
        self.assertRaises(ValueError, correlation_t, self.data1, self.data2,
                          confidence_level=0)
        self.assertRaises(ValueError, correlation_t, self.data1, self.data2,
                          confidence_level=0.0)
        self.assertRaises(ValueError, correlation_t, self.data1, self.data2,
                          confidence_level=1)
        self.assertRaises(ValueError, correlation_t, self.data1, self.data2,
                          confidence_level=1.0)

    def test_correlation_test_no_permutations(self):
        """Test correlation_t with no permutations."""
        # These results were verified with R.
        exp = (-0.2581988897471611, 0.7418011102528389, [], None,
               (-0.97687328610475876, 0.93488023560400879))
        obs = correlation_t([1, 2, 3, 4], [1, 2, 1, 1], permutations=0)
        assert_allclose(obs[0], exp[0])
        assert_allclose(obs[1], exp[1])
        assert_allclose(obs[2], exp[2])
        self.assertEqual(obs[3], exp[3])
        assert_allclose(obs[4], exp[4])

    def test_correlation_test_perfect_correlation(self):
        """Test correlation_t with perfectly-correlated input vectors."""
        # These results were verified with R.
        obs = correlation_t([1, 2, 3, 4], [1, 2, 3, 4])
        assert_allclose(obs[:2], (1.0, 0.0))
        self.assertEqual(len(obs[2]), 999)
        for r in obs[2]:
            self.assertTrue(r >= -1.0 and r <= 1.0)
        self.assertCorrectPValue(0.06, 0.09, correlation_t,
                                 ([1, 2, 3, 4], [1, 2, 3, 4]),
                                 p_val_idx=3)
        assert_allclose(obs[4], (0.99999999999998879, 1.0))

    def test_correlation_test_small_obs(self):
        """Test correlation_t with a small number of observations."""
        # These results were verified with R.
        obs = correlation_t([1, 2, 3], [1, 2, 3])
        assert_allclose(obs[:2], (1.0, 0))
        self.assertEqual(len(obs[2]), 999)
        for r in obs[2]:
            self.assertTrue(r >= -1.0 and r <= 1.0)
        self.assertCorrectPValue(0.3, 0.4, correlation_t,
                                 ([1, 2, 3], [1, 2, 3]),
                                 p_val_idx=3)
        self.assertEqual(obs[4], (None, None))

        obs = correlation_t([1, 2, 3], [1, 2, 3], method='spearman')
        assert_allclose(obs[:2], (1.0, 0))
        self.assertEqual(len(obs[2]), 999)
        for r in obs[2]:
            self.assertTrue(r >= -1.0 and r <= 1.0)
        self.assertCorrectPValue(0.3, 0.4, correlation_t,
                                 ([1, 2, 3], [1, 2, 3]),
                                 {'method': 'spearman'}, p_val_idx=3)
        self.assertEqual(obs[4], (None, None))

    def test_mw_test(self):
        """mann-whitney test results should match Sokal & Rohlf"""
        # using Sokal and Rolhf and R wilcox.test
        # x <- c(104, 109, 112, 114, 116, 118, 118, 119, 121, 123, 125, 126,
        #        126, 128, 128, 128)
        # y <- c(100, 105, 107, 107, 108, 111, 116, 120, 121, 123)
        # wilcox.test(x,y)
        # W = 123.5, p-value = 0.0232
        x = [104, 109, 112, 114, 116, 118, 118, 119, 121, 123, 125, 126, 126,
             128, 128, 128]
        y = [100, 105, 107, 107, 108, 111, 116, 120, 121, 123]
        u, p = mw_t(x, y, continuity=True, two_sided=True)
        # a return of 123.5 would also be okay, there is a consensus to use the
        # smaller U statistic, but the probability calculated from each is the
        # same
        self.assertTrue(u == 36.5 or u == 123.5)
        assert_allclose(p, .0232, rtol=1e-3)

    def test_mw_boot(self):
        """excercising the Monte-carlo variant of mann-whitney"""
        x = [104, 109, 112, 114, 116, 118, 118, 119, 121, 123, 125, 126, 126,
             128, 128, 128]
        y = [100, 105, 107, 107, 108, 111, 116, 120, 121, 123]
        u, p = mw_boot(x, y, 10)
        self.assertTrue(u == 36.5 or u == 123.5)
        self.assertTrue(0 <= p <= 0.5)

    def test_kendall(self):
        """tests new kendall tau implamentation, returns tau, prob"""
        # test from pg. 594 Sokal and Rohlf, Box 15.7
        v1 = [8.7, 8.5, 9.4, 10, 6.3, 7.8, 11.9, 6.5, 6.6, 10.6, 10.2, 7.2,
              8.6, 11.1, 11.6]
        v2 = [5.95, 5.65, 6.00, 5.70, 4.70, 5.53, 6.40, 4.18, 6.15, 5.93, 5.70,
              5.68, 6.13, 6.30, 6.03]
        obs_tau = kendall(v1, v2)
        obs_prob = kendall_pval(obs_tau, len(v1))
        exp_tau = 0.49761335152811925
        exp_prob = 0.0097188572446995618
        assert_allclose(obs_tau, exp_tau)
        assert_allclose(obs_prob, exp_prob)
        # random vectors checked against scipy. v1 has 33 ties, v2 32
        v1 = array(
            [1.2, 9.7, 8.8, 1.7, 8.6, 9.9, 6.8, 7.3, 5.5, 5.4, 8.3,
             3.6, 7.5, 2., 9.3, 5.1, 8.4, 0.3, 8.2, 2.4, 9.8, 8.5,
             2.1, 6., 1.8, 3.7, 1.4, 4.6, 7.6, 5.2, 0.9, 5.2, 4.7,
             2.9, 5., 6.9, 1.3, 6.7, 5.2, 2.4, 6.9, 2., 7.4, 0.4,
             8.2, 9.5, 2.9, 5.7, 2.4, 8.8, 1.6, 3.5, 5.1, 3.6, 3.3,
             7.5, 0.9, 9.3, 5.4, 6.9, 9.3, 2.3, 1.9, 8.1, 3.2, 4.2,
             8.7, 3., 9.8, 5.3, 6.2, 4.8, 9., 2.8, 5.5, 8.4, 4.1,
             5.6, 5.4, 6.9, 3.8, 2.7, 0.3, 3.9, 8.2, 6.6, 1.9, 3.9,
             2., 4.4, 0.8, 6.5, 4.8, 1.5, 9.9, 9.1, 9.9, 6.2, 2.9,
             2.])
        v2 = array([6.6, 8.6, 3.9, 6.1, 0.9, 8.4, 10., 3.3, 0.4,
                    3.9, 7.6, 8.2, 8.6, 3., 6.9, 0.6, 8.4, 8.1,
                    6.3, 0.5, 5.2, 6.4, 8., 9.9, 1.2, 6.7, 8.4,
                    2.7, 8.4, 4.1, 4.6, 5.1, 5.2, 5.3, 2.2, 2.2,
                    4.3, 7.1, 1.4, 6.6, 7.6, 4.5, 7.8, 3.5, 7.1,
                    0.6, 4.6, 3.2, 2.2, 0.2, 3.9, 5.9, 7.7, 8.8,
                    1.3, 5.1, 5.6, 8.3, 8.8, 1.7, 5.2, 6.9, 1.3,
                    1.4, 4.9, 9.4, 2.3, 3.7, 9.1, 3.4, 1.6, 4.1,
                    9.7, 2.8, 9.9, 0.5, 2., 2.7, 3.3, 2.4, 3.6,
                    7.9, 6.5, 7., 4.2, 1.8, 1.6, 1.9, 5.5, 0.,
                    1.4, 2.2, 7.2, 8.2, 1.1, 2.5, 5.3, 0.2, 9., 0.2])
        exp_tau, exp_prob = (0.024867511238807951, 0.71392573687923555)
        obs_tau = kendall(v1, v2)
        obs_prob = kendall_pval(obs_tau, len(v1))
        assert_allclose(obs_tau, exp_tau)
        assert_allclose(obs_prob, exp_prob)


class TestDistMatrixPermutationTest(TestCase):

    """Tests of distance_matrix_permutation_test"""

    def setUp(self):
        """sets up variables for testing"""
        self.matrix = array(
            [[1, 2, 3, 4], [5, 6, 7, 8], [9, 10, 11, 12], [13, 14, 15, 16]])
        self.cells = [(0, 1), (1, 3)]
        self.cells2 = [(0, 2), (2, 3)]

    def test_ANOVA_one_way(self):
        """ANOVA one way returns same values as ANOVA on a stats package
        """
        g1 = array([10.0, 11.0, 10.0, 5.0, 6.0])
        g2 = array([1.0, 2.0, 3.0, 4.0, 1.0, 2.0])
        g3 = array([6.0, 7.0, 5.0, 6.0, 7.0])
        i = [g1, g2, g3]
        F, pval = ANOVA_one_way(i)

        assert_allclose(F, 18.565450643776831)
        assert_allclose(pval, 0.00015486238993089464)

    def test_kruskal_wallis(self):
        """Test kruskal_wallis on Sokal & Rohlf Box 13.6 dataset"""
        d_control = [75, 67, 70, 75, 65, 71, 67, 67, 76, 68]
        d_2_gluc = [57, 58, 60, 59, 62, 60, 60, 57, 59, 61]
        d_2_fruc = [58, 61, 56, 58, 57, 56, 61, 60, 57, 58]
        d_1_1 = [58, 59, 58, 61, 57, 56, 58, 57, 57, 59]
        d_2_sucr = [62, 66, 65, 63, 64, 62, 65, 65, 62, 67]
        data = [d_control, d_2_gluc, d_2_fruc, d_1_1, d_2_sucr]
        kw_stat, pval = kruskal_wallis(data)
        assert_allclose(kw_stat, 38.436807439)
        assert_allclose(pval, 9.105424085598766e-08)
        # test using a random data set against scipy
        x_0 = array([0, 0, 0, 31, 12, 0, 25, 26, 775, 13])
        x_1 = array([14, 15, 0, 15, 12, 13])
        x_2 = array([0, 0, 0, 55, 92, 11, 11, 11, 555])
        # kruskal(x_0, x_1, x_2) = (0.10761259465923653, 0.94761564440615031)
        exp = (0.10761259465923653, 0.94761564440615031)
        obs = kruskal_wallis([x_0, x_1, x_2])
        assert_allclose(obs, exp)


class PvalueTests(TestCase):

    '''Test that the methods for handling Pvalues return the results we expect.

    Note: eps is being set lower on some of these because Sokal and Rohlf
    provide only ~5 sig figs and our integrals diverge by that much or more.
    '''

    def setUp(self):
        '''Nothing needed for all tests.'''
        pass

    def test_fdr_correction(self):
        """Test that the fdr_correction works as anticipated."""
        pvals = array([.1, .7, .5, .3, .9])
        exp = array([.5, .7 * 5 / 4., .5 * 5 / 3., .3 * 5 / 2., .9])
        obs = fdr_correction(pvals)
        assert_allclose(obs, exp)

    def test_benjamini_hochberg_step_down(self):
        """Test that the BH step down procedure behaves as it does in R."""
        # r values
        # q = c(0.64771481,  0.93517796,  0.7169902 ,  0.18223457,  0.26918556,
        #  0.1450153 ,  0.22448242,  0.74723508,  0.89061034,  0.74007906)
        # p.adjust(q, method='BH')
        #  [1] 0.9340439 0.9351780 0.9340439 0.6729639 0.6729639 0.6729639
        #      0.6729639
        #  [8] 0.9340439 0.9351780 0.9340439
        pvals = array([0.64771481, 0.93517796, 0.7169902, 0.18223457,
                       0.26918556, 0.1450153, 0.22448242, 0.74723508,
                       0.89061034, 0.74007906])
        exp = array([0.9340439, 0.9351780, 0.9340439, 0.6729639, 0.6729639,
                     0.6729639, 0.6729639, 0.9340439, 0.9351780, 0.9340439])
        obs = benjamini_hochberg_step_down(pvals)
        assert_allclose(obs, exp)
        # example 2
        pvals = array([1.32305426, 1.9345059, 0.87129877, 1.89957702,
                       1.85712616, 0.68757988, 0.41248969, 0.20751712,
                       1.97658599, 1.06209437])
        exp = array([1., 1., 1., 1., 1., 1., 1., 1., 1., 1.])
        obs = benjamini_hochberg_step_down(pvals)
        assert_allclose(obs, exp)

    def test_bonferroni_correction(self):
        """Test that Bonferroni correction behaves correctly."""
        pvals = array([.1, .7, .5, .3, .9])
        exp = pvals * 5.
        obs = bonferroni_correction(pvals)
        assert_allclose(obs, exp)

    def test_fisher_z_transform(self):
        '''Test Fisher Z transform is correct.'''
        r = .657
        exp = .5 * log(1.657 / .343)
        obs = fisher_z_transform(r)
        assert_allclose(exp, obs)
        r = 1
        obs = fisher_z_transform(r)
        assert_allclose(obs, nan)
        r = -1
        obs = fisher_z_transform(r)
        assert_allclose(obs, nan)
        r = -5.6
        obs = fisher_z_transform(r)
        assert_allclose(obs, nan)
        # from sokal and rohlf pg 575
        r = .972
        obs = fisher_z_transform(r)
        exp = 2.12730
        assert_allclose(exp, obs, rtol=1e-4)

    def test_z_transform_pval(self):
        '''Test that pval associated with Fisher Z is correct.'''
        r = .6
        n = 100
        obs = z_transform_pval(r, n)
        exp = 3.4353390341723208e-09
        assert_allclose(exp, obs)
        r = .5
        n = 3
        obs = z_transform_pval(r, n)
        assert_allclose(obs, nan)

    def test_inverse_fisher_z_transform(self):
        '''Test that Fisher's Z transform is computed correctly.'''
        z = .65
        exp = 0.5716699660851171
        obs = inverse_fisher_z_transform(z)
        assert_allclose(exp, obs)

    def test_fisher_population_correlation(self):
        '''Test that the population rho and homogeneity coeff are correct.'''
        # note: the error tolerances are lower than they would normally be
        # because sokal and rolhf don't give many significant figures
        # example from Sokal and Rohlf Biometry pg. 580 - 582
        rs = array([.29, .7, .58, .56, .55, .67, .65, .61, .64, .56])
        ns = array([100, 46, 28, 74, 33, 27, 52, 26, 20, 17])
        zbar = .615268
        X2 = 15.26352
        pop_r = .547825
        hval = chi2prob(X2, len(ns) - 1)
        obs_p_rho, obs_hval = fisher_population_correlation(rs, ns)
        assert_allclose(obs_p_rho, pop_r, rtol=1e-5)
        assert_allclose(obs_hval, hval, rtol=1e-5)
        # test with nans
        rs = array(
            [.29, .7, nan, .58, .56, .55, .67, .65, .61, .64, .56])
        ns = array([100, 46, 400, 28, 74, 33, 27, 52, 26, 20, 17])
        obs_p_rho, obs_hval = fisher_population_correlation(rs, ns)
        assert_allclose(obs_p_rho, pop_r, rtol=1e-5)
        assert_allclose(obs_hval, hval, rtol=1e-5)
        # test with short vectors
        rs = [.6, .5, .4, .6, .7]
        ns = [10, 12, 42, 11, 3]
        obs_p_rho, obs_hval = fisher_population_correlation(rs, ns)
        assert_allclose(obs_p_rho, nan)
        assert_allclose(obs_hval, nan)
        # test with data with rs >1
        rs = [.6, .5, .4, 1.4]
        ns = [10, 50, 100, 10]
        self.assertRaises(ValueError, fisher_population_correlation, rs, ns)

    def test_assign_correlation_pval(self):
        '''Test that correlation pvalues are assigned correctly with each meth.
        '''
        # test with parametric t distribution, use example from Sokal and Rohlf
        # Biometry pg 576.
        r = .86519
        n = 12
        ts = 5.45618  # only 5 sig figs in sokal and rohlf
        exp = tprob(ts, n - 2, tails='two-sided')
        obs = assign_correlation_pval(r, n, 'parametric_t_distribution')
        assert_allclose(exp, obs, rtol=1e-5)
        # test with too few samples
        n = 3
        self.assertRaises(ValueError, assign_correlation_pval, r, n,
                          'parametric_t_distribution')
        # test with fisher_z_transform
        r = .29
        n = 100
        z = 0.29856626366017841  # .2981 in biometry
        exp = z_transform_pval(z, n)
        obs = assign_correlation_pval(r, n, 'fisher_z_transform')
        assert_allclose(exp, obs, rtol=1e-5)
        r = .61
        n = 26
        z = 0.70892135942740819  # .7089 in biometry
        exp = z_transform_pval(z, n)
        obs = assign_correlation_pval(r, n, 'fisher_z_transform')
        assert_allclose(exp, obs, rtol=1e-5)
        # prove that we can have specify the other options, and as long as we
        # dont have bootstrapped selected we are fine.
        v1 = array([10, 11, 12])
        v2 = array([10, 14, 15])
        obs = assign_correlation_pval(r, n, 'fisher_z_transform',
                                      permutations=1000, perm_test_fn=pearson,
                                      v1=v1, v2=v2)
        assert_allclose(exp, obs)
        # test with bootstrapping, seed for reproducibility.
        seed(0)
        v1 = array([54, 71, 60, 54, 42, 64, 43, 89, 96, 38])
        v2 = array([79, 52, 56, 92, 7, 8, 2, 83, 77, 87])
        # c = corrcoef(v1,v2)[0][1]
        exp = .357
        obs = assign_correlation_pval(0.33112494, 20000, 'bootstrapped',
                                      permutations=1000, perm_test_fn=pearson,
                                      v1=v1, v2=v2)
        assert_allclose(exp, obs)
        # make sure it throws an error
        self.assertRaises(ValueError, assign_correlation_pval, 7, 20000,
                          'bootstrapped', perm_test_fn=pearson, v1=None, v2=v2)
        # test that it does properly with kendall
        exp = kendall_pval(r, n)
        obs = assign_correlation_pval(r, n, 'kendall')
        assert_allclose(exp, obs)

    def test_cscore(self):
        '''Test cscore is calculated correctly.'''
        # test using example from Stone and Roberts pg 75
        v1 = array([1, 0, 0, 0, 1, 1, 0, 1, 0, 1])
        v2 = array([1, 1, 1, 0, 1, 0, 1, 1, 1, 0])
        obs = cscore(v1, v2)
        exp = 8
        self.assertEqual(obs, exp)
        # test using examples verified in ecosim
        v1 = array([4, 6, 12, 13, 14, 0, 0, 0, 14, 11, 9, 6, 0, 1, 1, 0, 0,
                    4])
        v2 = array([4, 0, 0, 113, 1, 2, 20, 0, 1, 0, 19, 16, 0, 13, 6, 0, 5,
                    4])
        # from R
        # library(vegan)
        # library(bipartite)
        # m = matrix(c(4,6,12,13,14,0,0,0,14,11,9,6,0,1,1,0,0,4,4,0,0,113,1,2,
        #              20,0,1,0,19,16,0,13,6,0,5,4), 18,2)
        # C.score(m, normalise=FALSE)
        exp = 9
        obs = cscore(v1, v2)
        self.assertEqual(obs, exp)


class DistributionTests(TestCase):
    '''Test that the distributions from scipy are perfoming as we expect.'''

    def setUp(self):
        '''Nothing needed for all tests.'''
        pass

    def test_normal_probability_distribution(self):
        '''Test that the normal probability distribution performs correctly.'''
        # test against R
        # library('stats')
        # pnorm(4.5, mean = 0, sd=1, lower.tail=TRUE)
        # 0.9999966
        p = normprob(4.5, direction='low', mean=0, std=1)
        assert_allclose(p, 0.9999966)
        # pnorm(-14.5, mean = -5, sd=20, lower.tail=FALSE)
        # 0.3173935
        p = normprob(-14.5, direction='two-sided', mean=-5, std=20)
        assert_allclose(p, 0.3173935*2)
        # > pnorm(4.5, mean = 0, sd=1, lower.tail=FALSE)
        # [1] 3.397673e-06
        p = normprob(4.5, direction='high', mean=0, std=1)
        assert_allclose(p, 3.397673e-06)
        p = normprob(4.5, direction='two-sided', mean=0, std=1)
        assert_allclose(p, 3.397673e-06*2)
        # test that a ValueError is correctly raised
        self.assertRaises(ValueError, normprob, 4.5, direction='dne')

    def test_chi2_probability_distribution(self):
        '''Test that chi2 probability distribution performs correctly.'''
        # test against R
        # library('stats')
        # pchisq(13.4, 4, lower.tail=TRUE)
        # 0.990522
        p = chi2prob(13.4, 4, direction='low')
        assert_allclose(p, 0.990522)
        # > pchisq(13.4, 4, lower.tail=FALSE)
        # [1] 0.009478022
        p = chi2prob(13.4, 4, direction='high')
        assert_allclose(p, 0.009478022)
        # test when we have a negative chi2 stat
        p = chi2prob(-10, 5, direction='high')
        assert_allclose(p, nan)
        # test another value
        # > pchisq(45, 35)
        # [1] 0.8800662
        p = chi2prob(45, 35, direction='low')
        assert_allclose(p, 0.8800662)
        # test that a ValueError is correctly raised
        self.assertRaises(ValueError, chi2prob, 4.5, 3, direction='dne')

    def test_t_probability_distribution(self):
        '''Test that the t probability distribution performs correctly.'''
        # test against R
        # library('stats')
        # pt(2.5, 10, lower.tail=TRUE)
        # 0.9842766
        t = tprob(2.5, 10, tails='low')
        assert_allclose(t, 0.9842766, atol=1e-7)
        # pt(2.5, 10, lower.tail=FALSE)
        # 0.01572342
        t = tprob(2.5, 10, tails='high')
        assert_allclose(t, 0.01572342, atol=1e-7)
        # both tails
        t = tprob(2.5, 10, tails='two-sided')
        assert_allclose(t, 2*0.01572342, atol=1e-7)
        # > pt(-6.7,2)
        # [1] 0.01077945
        t = tprob(-6.7, 2, tails='two-sided')
        assert_allclose(t, 2*0.01077945, atol=1e-7)
        # test that a ValueError is correctly raised
        self.assertRaises(ValueError, tprob, 4.5, 3, tails='dne')

    def test_f_probability_distribution(self):
        '''Test that the f probability distribution performs correctly.'''
        # test against R
        # library('stats')
        # pf(4.5, 3, 5)
        # 0.9305489
        p = fprob(4.5, 3, 5, direction='low')
        assert_allclose(p, 0.9305489, atol=1e-7)
        p = fprob(4.5, 3, 5, direction='high')
        assert_allclose(p, 1 - 0.9305489, atol=1e-7)
        # pf(33.5, 2, 5)
        # 0.9987292
        p = fprob(33.5, 2, 5, direction='low')
        assert_allclose(p, 0.9987292, atol=1e-7)
        # test when we have a negative f stat
        p = fprob(-10, 5, 6, direction='high')
        assert_allclose(p, nan)
        # test that a ValueError is correctly raised
        self.assertRaises(ValueError, fprob, 4.5, 3, 5, direction='dne')

if __name__ == "__main__":
    main()
