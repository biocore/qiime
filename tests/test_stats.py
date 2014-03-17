#!/usr/bin/env python
from __future__ import division

__author__ = "Michael Dwan"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Michael Dwan", "Logan Knecht",
               "Damien Coy", "Levi McCracken", "Andrew Cochran"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Test suite for classes, methods and functions of the stats module."""

from shutil import rmtree
from os.path import exists, join
from string import digits
from unittest import TestCase, main
from numpy.testing import assert_almost_equal
from cogent.util.misc import remove_files, create_dir
from numpy import array, asarray, roll, median, nan
from numpy.random import permutation, shuffle
import numpy as np
from itertools import izip
from types import StringType, ListType, FloatType, TupleType
from biom.parse import parse_biom_table
from skbio.core.distance import DistanceMatrix, SymmetricDistanceMatrix
from qiime.stats import (all_pairs_t_test, _perform_pairwise_tests,
                         Anosim, Best, CategoryStats, CorrelationStats,
                         DistanceMatrixStats, MantelCorrelogram, Mantel,
                         PartialMantel, Permanova, quantile, _quantile,
                         paired_difference_analyses)
from qiime.util import MetadataMap, get_qiime_temp_dir, get_tmp_filename


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
        elif observed is not None and isinstance(observed, (np.number, np.ndarray, FloatType)):
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
        self.overview_dm = SymmetricDistanceMatrix.from_file(
            self.overview_dm_str)

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
        self.single_ele_dm = SymmetricDistanceMatrix([[0]], ['s1'])

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
        np.random.seed(self.value_for_seed)
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
        exp = [['foo', 'bar', nan, nan, nan, nan, nan], ['foo', 'baz',
                                                         -
                                                         7.794228634059948, 0.008032650971672552, 0.016065301943345104,
                                                         nan, nan], ['bar', 'baz',
                                                                     -
                                                                     2.598076211353316, 0.060844967173160069, 0.12168993434632014,
                                                                     nan, nan]]
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
                          [DistanceMatrix(
                              array([[0, 2], [3, 0]]), ['foo', 'bar']),
                           DistanceMatrix(
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
                          [DistanceMatrix(
                              array([[0, 2], [3, 0]]), ['foo', 'bar']),
                           DistanceMatrix(
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
        mismatch = SymmetricDistanceMatrix(array([[0]]), ['s2'])

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


class CategoryStatsTests(TestHelper):
    """Tests for the CategoryStats class."""

    def setUp(self):
        """Define some useful data to use in testing."""
        super(CategoryStatsTests, self).setUp()
        self.cs_overview = CategoryStats(self.overview_map, [self.overview_dm],
                                         ["Treatment", "DOB"])

    def test_MetadataMap_setter(self):
        """Should set the mdmap property."""
        self.cs_overview.MetadataMap = self.overview_map
        self.assertEqual(self.cs_overview.MetadataMap, self.overview_map)

    def test_MetadataMap_setter_invalid_input(self):
        """Setter must receive the correct and compatible object types."""
        self.assertRaises(TypeError, setattr, self.cs_overview, 'MetadataMap',
                          "foo")
        self.assertRaises(TypeError, setattr, self.cs_overview, 'MetadataMap',
                          [])
        self.assertRaises(TypeError, setattr, self.cs_overview, 'MetadataMap',
                          {})
        self.assertRaises(TypeError, setattr, self.cs_overview, 'MetadataMap',
                          None)
        self.assertRaises(TypeError, setattr, self.cs_overview, 'MetadataMap',
                          self.overview_dm)

    def test_MetadataMap_getter(self):
        """Test valid return of MetadataMap property."""
        self.assertEqual(self.cs_overview.MetadataMap, self.overview_map)

    def test_Categories_setter_invalid_input(self):
        """Must receive a list of categories that are in the mapping file."""
        self.assertRaises(TypeError, setattr, self.cs_overview, 'Categories',
                          "Hello!")
        self.assertRaises(TypeError, setattr, self.cs_overview, 'Categories',
                          self.overview_dm)
        self.assertRaises(ValueError, setattr, self.cs_overview, 'Categories',
                          ["hehehe", 123, "hello"])
        self.assertRaises(ValueError, setattr, self.cs_overview, 'Categories',
                          ["foo"])

        # Test setting a unique category.
        self.assertRaises(ValueError, CategoryStats, self.test_map,
                          [self.overview_dm], ["Description"])
        cs_test = CategoryStats(self.test_map, [self.overview_dm], ["Foo"])
        self.assertRaises(ValueError, setattr, cs_test, 'Categories',
                          ["Description"])

        # Test setting a category with only a single value.
        self.assertRaises(ValueError, setattr, cs_test, 'Categories', ["Bar"])

        # Test setting a category that is non-numeric.
        self.assertRaises(TypeError, CategoryStats, self.overview_map,
                          [self.overview_dm], ["Treatment"],
                          suppress_numeric_category_check=False)

    def test_Categories_getter(self):
        """Test valid return of Categories property."""
        expected = ['Treatment', 'DOB']
        observed = self.cs_overview.Categories
        self.assertEqual(observed, expected)

    def test_RandomFunction_getter(self):
        """Test retrieval of a random function reference."""
        self.assertEqual(self.cs_overview.RandomFunction, permutation)

    def test_RandomFunction_setter(self):
        """Test setter for the random function to use in p-value calc."""
        self.assertEqual(self.cs_overview.RandomFunction, permutation)
        nrs = NonRandomShuffler()
        self.cs_overview.RandomFunction = nrs.permutation
        self.assertEqual(self.cs_overview.RandomFunction, nrs.permutation)

    def test_RandomFunction_setter_invalid_input(self):
        """Test setter for the random function with non-callable input."""
        self.assertRaises(TypeError, setattr, self.cs_overview,
                          'RandomFunction', 42)
        self.assertRaises(TypeError, setattr, self.cs_overview,
                          'RandomFunction', 42.0)
        self.assertRaises(TypeError, setattr, self.cs_overview,
                          'RandomFunction', "j")
        self.assertRaises(TypeError, setattr, self.cs_overview,
                          'RandomFunction', None)
        self.assertRaises(TypeError, setattr, self.cs_overview,
                          'RandomFunction', [])
        self.assertRaises(TypeError, setattr, self.cs_overview,
                          'RandomFunction', ())
        self.assertRaises(TypeError, setattr, self.cs_overview,
                          'RandomFunction', {})

    def test_validate_compatibility(self):
        """Test for compatible sample IDs between dms and mdmap."""
        self.assertEqual(self.cs_overview._validate_compatibility(), None)

        self.cs_overview.DistanceMatrices = [self.single_ele_dm]
        self.assertRaises(ValueError, self.cs_overview._validate_compatibility)

        self.cs_overview.DistanceMatrices = [self.overview_dm]
        self.cs_overview.MetadataMap = self.test_map
        self.assertRaises(ValueError, self.cs_overview._validate_compatibility)

    def test_call(self):
        """Test _call__() returns an empty result set."""
        self.assertEqual(self.cs_overview(), {})
        self.assertEqual(self.cs_overview(10), {})

    def test_call_bad_perms(self):
        """Test __call__() fails upon receiving invalid number of perms."""
        self.assertRaises(ValueError, self.cs_overview, -1)

    def test_call_incompatible_data(self):
        """Test __call__() fails after incompatible dms/mdmap pair is set."""
        self.cs_overview.DistanceMatrices = [self.single_ele_dm,
                                             self.single_ele_dm]
        self.assertRaises(ValueError, self.cs_overview)


class AnosimTests(TestHelper):
    """Tests for the Anosim class.

    This testing code is heavily based on Andrew Cochran's original suite of
    unit tests for ANOSIM.
    """

    def setUp(self):
        """Define some useful data to use in testing."""
        super(AnosimTests, self).setUp()

        # Define two small dms for easy testing. One has ties in the ranks.
        self.small_dm_str = ["\tsam1\tsam2\tsam3\tsam4",
                             "sam1\t0\t1\t5\t4",
                             "sam2\t1\t0\t3\t2",
                             "sam3\t5\t3\t0\t3",
                             "sam4\t4\t2\t3\t0"]
        self.small_dm = SymmetricDistanceMatrix.from_file(self.small_dm_str)

        self.small_dm_tie_str = ["\tsam1\tsam2\tsam3\tsam4",
                                 "sam1\t0\t1\t1\t4",
                                 "sam2\t1\t0\t3\t2",
                                 "sam3\t1\t3\t0\t3",
                                 "sam4\t4\t2\t3\t0"]
        self.small_dm_tie = SymmetricDistanceMatrix.from_file(
            self.small_dm_tie_str)

        self.small_map_str = ["#SampleID\tBarcodeSequence\
                              \tLinkerPrimerSequence\tTreatment\tDOB\
                              \tDescription",
                              "sam1\tAGCACGAGCCTA\tYATGCTGCCTCCCGTAGGAGT\
                              \tControl\t20061218\tControl_mouse_I.D._354",
                              "sam2\tAACTCGTCGATG\tYATGCTGCCTCCCGTAGGAGT\
                              \tControl\t20061218\tControl_mouse_I.D._355",
                              "sam3\tACAGACCACTCA\tYATGCTGCCTCCCGTAGGAGT\
                              \tFast\t20061126\tControl_mouse_I.D._356",
                              "sam4\tACCAGCGACTAG\tYATGCTGCCTCCCGTAGGAGT\
                              \tFast\t20070314\tControl_mouse_I.D._481"]
        self.small_map = MetadataMap.parseMetadataMap(self.small_map_str)

        # Create a group map, which maps sample ID to category value (e.g.
        # sample 1 to 'control' and sample 2 to 'fast'). This comes in handy
        # for testing some of the private methods in the Anosim class. This
        # group map can be used for testing both the small dm data and the
        # small dm with ties data.
        self.small_group_map = {}
        for samp_id in self.small_dm.ids:
            self.small_group_map[samp_id] = self.small_map.getCategoryValue(
                samp_id, 'Treatment')

        # Create three Anosim instances: one for the small dm, one for the
        # small dm with ties, and one for the overview tutorial dataset.
        self.anosim_small = Anosim(self.small_map, self.small_dm, 'Treatment')
        self.anosim_small_tie = Anosim(self.small_map, self.small_dm_tie,
                                       'Treatment')
        self.anosim_overview = Anosim(self.overview_map, self.overview_dm,
                                      'Treatment')

    def test_call_overview(self):
        """Test __call__() on overview data with Treatment category."""
        # These results were verified with R.
        exp = {'method_name': 'ANOSIM', 'p_value': 0.0080000000000000002,
               'r_value': 0.8125}
        obs = self.anosim_overview()
        self.assertEqual(obs['method_name'], exp['method_name'])
        assert_almost_equal(obs['r_value'], exp['r_value'])
        self.assertCorrectPValue(0, 0.06, self.anosim_overview)

    def test_call_small(self):
        """Test __call__() on small dm."""
        # These results were verified with R.
        exp = {'method_name': 'ANOSIM', 'p_value': 0.31, 'r_value': 0.625}
        obs = self.anosim_small()

        self.assertEqual(obs['method_name'], exp['method_name'])
        assert_almost_equal(obs['r_value'], exp['r_value'])
        self.assertCorrectPValue(0.28, 0.42, self.anosim_small)

    def test_call_small_ties(self):
        """Test __call__() on small dm with ties in ranks."""
        # These results were verified with R.
        exp = {'method_name': 'ANOSIM', 'p_value': 0.67600000000000005,
               'r_value': 0.25}
        obs = self.anosim_small_tie()

        self.assertEqual(obs['method_name'], exp['method_name'])
        assert_almost_equal(obs['r_value'], exp['r_value'])
        self.assertCorrectPValue(0.56, 0.75, self.anosim_small_tie)

    def test_call_no_perms(self):
        """Test __call__() on small dm with no permutations."""
        # These results were verified with R.
        exp = {'method_name': 'ANOSIM', 'p_value': 1.0, 'r_value': 0.625}
        obs = self.anosim_small(0)

        self.assertEqual(obs['method_name'], exp['method_name'])
        assert_almost_equal(obs['r_value'], exp['r_value'])
        assert_almost_equal(obs['p_value'], exp['p_value'])

    def test_call_incompatible_data(self):
        """Should fail on incompatible mdmap/dm combo and bad perms."""
        self.assertRaises(ValueError, self.anosim_small, -1)
        self.anosim_small.DistanceMatrices = [self.single_ele_dm]
        self.assertRaises(ValueError, self.anosim_small)

    def test_anosim_small(self):
        """Test _anosim() on small dm."""
        # These results were verified with R.
        exp = 0.625
        obs = self.anosim_small._anosim(self.small_group_map)
        assert_almost_equal(obs, exp)

    def test_anosim_small_ties(self):
        """Test _anosim() on small dm with ties."""
        # These results were verified with R.
        exp = 0.25
        obs = self.anosim_small_tie._anosim(self.small_group_map)
        assert_almost_equal(obs, exp)

    def test_remove_ties1(self):
        """Test removal of ties. Should return [1.5,1.5]."""
        result = self.anosim_small._remove_ties([1, 1], [1, 2])
        self.assertEqual(result, [1.5, 1.5])

    def test_remove_ties2(self):
        """Should return [3.5,3.5,3.5,3.5,3.5,3.5]."""
        result = self.anosim_small._remove_ties(
            [1, 1, 1, 1, 1, 1], [1, 2, 3, 4, 5, 6])
        self.assertEqual(result, [3.5, 3.5, 3.5, 3.5, 3.5, 3.5])

    def test_remove_ties3(self):
        """Should return [1,3.5,3.5,3.5,3.5,6]."""
        result = self.anosim_small._remove_ties(
            [1, 3, 3, 3, 3, 8], [1, 2, 3, 4, 5, 6])
        self.assertEqual(result, [1, 3.5, 3.5, 3.5, 3.5, 6])

    def test_remove_ties4(self):
        """Should return [1,2,3,4]."""
        result = self.anosim_small._remove_ties([1, 2, 3, 4], [1, 2, 3, 4])
        self.assertEqual(result, [1, 2, 3, 4])

    def test_remove_ties5(self):
        """Should return [1,3,3,3,5.5,5.5,7]."""
        result = self.anosim_small._remove_ties([1, 2, 2, 2, 3, 3, 5],
                                                [1, 2, 3, 4, 5, 6, 7])
        self.assertEqual(result, [1, 3, 3, 3, 5.5, 5.5, 7])

    def test_remove_ties6(self):
        """Should return [1.5,1.5,3.5,3.5]."""
        result = self.anosim_small._remove_ties([1, 1, 2, 2], [1, 2, 3, 4])
        self.assertEqual(result, [1.5, 1.5, 3.5, 3.5])

    def test_get_adjusted_vals(self):
        """Test computing adjusted ranks for ties."""
        exp = [4, 4, 4]
        obs = self.anosim_small._get_adjusted_vals([3, 4, 5], 0, 2)
        self.assertEqual(obs, exp)

    def test_compute_r1(self):
        """Should return .625 for the R statistic on the small dm."""
        sorted_rank = [1.0, 2.0, 3.5, 3.5, 5.0, 6.0]
        sorted_group = [1.0, 0.0, 0.0, 1.0, 0.0, 0.0]
        sorted_rank = array(sorted_rank)
        sorted_group = array(sorted_group)
        result = self.anosim_small._compute_r_value(
            sorted_rank,
            sorted_group,
            4)
        self.assertEqual(result, .625)

    def test_anosim_p_test(self):
        """p-value should be .5 for this test."""
        nrs = NonRandomShuffler()
        self.anosim_small.RandomFunction = nrs.permutation

        exp = {'method_name': 'ANOSIM', 'p_value': 0.5, 'r_value': 0.625}
        obs = self.anosim_small(3)

        self.assertEqual(obs['method_name'], exp['method_name'])
        assert_almost_equal(obs['r_value'], exp['r_value'])
        assert_almost_equal(obs['p_value'], exp['p_value'])


class PermanovaTests(TestHelper):

    def setUp(self):
        """Define some useful data to use in testing."""
        super(PermanovaTests, self).setUp()

        # Some distance matrices to help test Permanova.
        self.distmtx_str = ["\tsam1\tsam2\tsam3\tsam4",
                            "sam1\t0\t1\t5\t4",
                            "sam2\t1\t0\t3\t2",
                            "sam3\t5\t3\t0\t3",
                            "sam4\t4\t2\t3\t0"]

        self.distmtx = SymmetricDistanceMatrix.from_file(self.distmtx_str)
        self.distmtx_samples = self.distmtx.ids

        self.distmtx_tie_str = ["\tsam1\tsam2\tsam3\tsam4",
                                "sam1\t0\t1\t1\t4",
                                "sam2\t1\t0\t3\t2",
                                "sam3\t1\t3\t0\t3",
                                "sam4\t4\t2\t3\t0"]
        self.distmtx_tie = SymmetricDistanceMatrix.from_file(
            self.distmtx_tie_str)
        self.distmtx_tie_samples = self.distmtx_tie.ids

        # For testing with uneven group sizes.
        self.distmtx_uneven_str = ["\tsam1\tsam2\tsam3\tsam4\tsam5",
                                   "sam1\t0\t3\t7\t2\t1",
                                   "sam2\t3\t0\t5\t4\t1",
                                   "sam3\t7\t5\t0\t2\t6",
                                   "sam4\t2\t4\t2\t0\t2",
                                   "sam5\t1\t1\t6\t2\t0"]
        self.distmtx_uneven = SymmetricDistanceMatrix.from_file(
            self.distmtx_uneven_str)
        self.distmtx_uneven_samples = self.distmtx_uneven.ids

        # Some group maps to help test Permanova, data_map can be used with
        # distmtx and distmtx_tie while data_map_uneven can only be used
        # with distmtx_uneven.
        self.data_map_str = ["#SampleID\tBarcodeSequence\tLinkerPrimerSequence\
                \tTreatment\tDOB\tDescription",
                             "sam1\tAGCACGAGCCTA\tYATGCTGCCTCCCGTAGGAGT\tControl\t20061218\
                \tControl_mouse_I.D._354",
                             "sam2\tAACTCGTCGATG\tYATGCTGCCTCCCGTAGGAGT\tControl\t20061218\
                \tControl_mouse_I.D._355",
                             "sam3\tACAGACCACTCA\tYATGCTGCCTCCCGTAGGAGT\tFast\t20061126\
                \tControl_mouse_I.D._356",
                             "sam4\tACCAGCGACTAG\tYATGCTGCCTCCCGTAGGAGT\tFast\t20070314\
                \tControl_mouse_I.D._481"]
        self.data_map = MetadataMap.parseMetadataMap(self.data_map_str)

        # For testing with uneven group sizes.
        self.data_map_uneven_str = ["#SampleID\tBarcodeSequence\
                \tLinkerPrimerSequence\tTreatment\tDOB\tDescription",
                                    "sam1\tAGCACGAGCCTA\tYATGCTGCCTCCCGTAGGAGT\tControl\t20061218\
                \tControl_mouse_I.D._354",
                                    "sam2\tAACTCGTCGATG\tYATGCTGCCTCCCGTAGGAGT\tControl\t20061218\
                \tControl_mouse_I.D._355",
                                    "sam3\tACAGACCACTCA\tYATGCTGCCTCCCGTAGGAGT\tFast\t20061126\
                \tControl_mouse_I.D._356",
                                    "sam4\tACCAGCGACTAG\tYATGCTGCCTCCCGTAGGAGT\tAwesome\t20070314\
                \tControl_mouse_I.D._481",
                                    "sam5\tACCAGCGACTAG\tYATGCTGCCTCCCCTATADST\tAwesome\t202020\
                \tcontrolmouseid"]
        self.data_map_uneven = MetadataMap.parseMetadataMap(
            self.data_map_uneven_str)

        # Formatting the two data_maps to meet permanova requirements.
        self.map = {}
        for samp_id in self.data_map.SampleIds:
            self.map[samp_id] = self.data_map.getCategoryValue(samp_id,
                                                               'Treatment')
        self.map_uneven = {}
        for samp_id in self.data_map_uneven.SampleIds:
            self.map_uneven[samp_id] = self.data_map_uneven.getCategoryValue(
                samp_id, 'Treatment')

        # Creating instances of Permanova to run the tests on.
        self.permanova_plain = Permanova(self.data_map, self.distmtx,
                                         'Treatment')
        self.permanova_tie = Permanova(self.data_map, self.distmtx_tie,
                                       'Treatment')
        self.permanova_uneven = Permanova(self.data_map_uneven,
                                          self.distmtx_uneven, 'Treatment')
        self.permanova_overview = Permanova(self.overview_map,
                                            self.overview_dm, 'Treatment')

    def test_permanova1(self):
        """permanova() should return 4.4."""
        exp = 4.4
        obs = self.permanova_plain._permanova(self.map)
        self.assertEqual(obs, exp)

    def test_permanova2(self):
        """Should result in 2."""
        exp = 2
        obs = self.permanova_tie._permanova(self.map)
        self.assertEqual(obs, exp)

    def test_permanova3(self):
        """Should result in 3.58462."""
        exp = 3.58462
        obs = self.permanova_uneven._permanova(self.map_uneven)
        assert_almost_equal(obs, exp, decimal=4)

    def test_compute_f1(self):
        """Should return 4.4, testing just function."""
        distances = [1, 5, 4, 3, 2, 3]
        grouping = [0, -1, -1, -1, -1, 1]
        distances = array(distances)
        grouping = array(grouping)
        result = self.permanova_plain._compute_f_value(
            distances, grouping, 4, 2,
            [2, 2])
        self.assertEqual(result, 4.4)

    def test_call_plain(self):
        """Test __call__() on plain dm."""
        exp = {'method_name': 'PERMANOVA', 'p_value': "?", 'f_value': 4.4}
        obs = self.permanova_plain()

        self.assertEqual(obs['method_name'], exp['method_name'])
        assert_almost_equal(obs['f_value'], exp['f_value'])
        self.assertCorrectPValue(0.28, 0.42, self.permanova_plain)

    def test_call_tie(self):
        """Test __call__() on dm with ties in ranks."""
        exp = {'method_name': 'PERMANOVA', 'p_value': "?", 'f_value': 2}
        obs = self.permanova_tie()

        self.assertEqual(obs['method_name'], exp['method_name'])
        assert_almost_equal(obs['f_value'], exp['f_value'])
        self.assertCorrectPValue(0.56, 0.75, self.permanova_tie)

    def test_call_uneven(self):
        """Test __call__() on uneven group sizes with no permutations."""
        exp = {'method_name': 'PERMANOVA', 'p_value': 1.0, 'f_value': 3.58462}
        obs = self.permanova_uneven(0)

        self.assertEqual(obs['method_name'], exp['method_name'])
        assert_almost_equal(obs['f_value'], exp['f_value'], decimal=4)
        assert_almost_equal(obs['p_value'], exp['p_value'])

    def test_call_overview(self):
        """Test __call__() on the overview dataset."""
        exp = {'method_name': 'PERMANOVA', 'p_value': 0.039215686274509803,
               'f_value': 2.2966506517077487}
        obs = self.permanova_overview(50)

        self.assertEqual(obs['method_name'], exp['method_name'])
        assert_almost_equal(obs['f_value'], exp['f_value'])
        self.assertCorrectPValue(0.005, 0.07, self.permanova_overview, 50)

    def test_call_incompatible_data(self):
        """Should fail on incompatible mdmap/dm combo and bad perms."""
        self.assertRaises(ValueError, self.permanova_plain, -1)
        self.permanova_plain.DistanceMatrices = [self.single_ele_dm]
        self.assertRaises(ValueError, self.permanova_plain)


class BestTests(TestHelper):
    """Tests for the Best class."""

    def setUp(self):
        """Define some useful data to use in testing."""
        super(BestTests, self).setUp()

        self.bv_dm_88soils_str = ["\tMT2.141698\tCA1.141704\tBB2.141659\t"
                                  "CO2.141657\tTL3.141709\tSN3.141650", "MT2.141698\t0.0\t"
                                  "0.623818643706\t0.750015427505\t0.585201193913\t0.729023583672\t"
                                  "0.622135587669", "CA1.141704\t0.623818643706\t0.0\t0.774881224555"
                                  "\t0.649822398416\t0.777203137034\t0.629507320436", "BB2.141659\t"
                                  "0.750015427505\t0.774881224555\t0.0\t0.688845424001\t0.567470311282"
                                  "\t0.721707516043", "CO2.141657\t0.585201193913\t0.649822398416\t"
                                  "0.688845424001\t0.0\t0.658853575764\t0.661223617505", "TL3.141709\t"
                                  "0.729023583672\t0.777203137034\t0.567470311282\t0.658853575764\t0.0\t"
                                  "0.711173405838", "SN3.141650\t0.622135587669\t0.629507320436\t"
                                  "0.721707516043\t0.661223617505\t0.711173405838\t0.0"]
        self.bv_dm_88soils = SymmetricDistanceMatrix.from_file(
            self.bv_dm_88soils_str)

        self.bv_map_88soils_str = ["#SampleId\tTOT_ORG_CARB\tSILT_CLAY\t"
                                   "ELEVATION\tSOIL_MOISTURE_DEFICIT\tCARB_NITRO_RATIO\t"
                                   "ANNUAL_SEASON_TEMP\tANNUAL_SEASON_PRECPT\tPH\tCMIN_RATE\tLONGITUDE\t"
                                   "LATITUDE", "MT2.141698\t39.1\t35\t1000\t70\t23.087\t7\t450\t6.66\t"
                                   "19.7\t-114\t46.8", "CA1.141704\t16.7\t73\t2003\t198\t13\t10.3\t400\t"
                                   "7.27\t2.276\t-111.7666667\t36.05", "BB2.141659\t52.2\t44\t400\t-680\t"
                                   "21.4\t6.1\t1200\t4.6\t2.223\t-68.1\t44.86666667", "CO2.141657\t18.1\t"
                                   "24\t2400\t104\t31.8\t6.1\t350\t5.68\t9.223\t-105.3333333\t"
                                   "40.58333333", "TL3.141709\t53.9\t52\t894\t-212\t24.6\t-9.3\t400\t"
                                   "4.23\t16.456\t-149.5833333\t68.63333333", "SN3.141650\t16.6\t20\t"
                                   "3000\t-252\t13.9\t3.6\t600\t5.74\t6.289\t-118.1666667\t36.45"]
        self.bv_map_88soils = MetadataMap.parseMetadataMap(
            self.bv_map_88soils_str)

        self.cats = ['TOT_ORG_CARB', 'SILT_CLAY', 'ELEVATION',
                     'SOIL_MOISTURE_DEFICIT', 'CARB_NITRO_RATIO',
                     'ANNUAL_SEASON_TEMP', 'ANNUAL_SEASON_PRECPT', 'PH',
                     'CMIN_RATE', 'LONGITUDE', 'LATITUDE']
        self.best = Best(self.bv_dm_88soils, self.bv_map_88soils, self.cats)

    def test_vector_dist(self):
        """Test the _vector_dist helper method."""
        v1 = [1, 4, 2]
        v2 = [-1, 12, 4]

        exp = 8.48528137424
        obs = self.best._vector_dist(v1, v2)
        assert_almost_equal(exp, obs)

        v1 = [1, 2, 100, 4, 2]
        v2 = [-1, 12, 4, 12, 99]

        exp = 137.087563258
        obs = self.best._vector_dist(v1, v2)
        assert_almost_equal(exp, obs)

    def test_make_cat_mat(self):
        """Test the _make_cat_mat method."""
        exp = [('6.66', '46.8'),
               ('7.27', '36.05'),
               ('4.6', '44.86666667'),
               ('5.68', '40.58333333'),
               ('4.23', '68.63333333'),
               ('5.74', '36.45')]

        obs = self.best._make_cat_mat(['PH', 'LONGITUDE', 'LATITUDE'], [1, 3])
        self.assertEqual(exp, obs)

        exp = [('6.66', '-114', '46.8'),
               ('7.27', '-111.7666667', '36.05'),
               ('4.6', '-68.1', '44.86666667'),
               ('5.68', '-105.3333333', '40.58333333'),
               ('4.23', '-149.5833333', '68.63333333'),
               ('5.74', '-118.1666667', '36.45')]

        obs = self.best._make_cat_mat(['PH', 'LONGITUDE', 'LATITUDE'],
                                      [1, 2, 3])
        self.assertEqual(exp, obs)

    def test_derive_euclidean_dm(self):
        """Test the derive_euclidean_dm method."""
        cat_mat = [('6.66', '46.8'),
                   ('7.27', '36.05'),
                   ('4.6', '44.86666667'),
                   ('5.68', '40.58333333'),
                   ('4.23', '68.63333333'),
                   ('5.74', '36.45')]

        dm_lbls = ['MT2.141698', 'CA1.141704', 'BB2.141659',
                   'CO2.141657', 'TL3.141709', 'SN3.141650']

        mtx = [
            [0.0, 10.7672930674, 2.82513322958, 6.29343661968,
             21.9681438519, 10.3908084382],
            [10.7672930674, 0.0, 9.21208506093, 4.80408275125,
             32.7248408842, 1.58142340946],
            [2.82513322958, 9.21208506093, 0.0, 4.41739114202,
             23.7695465697, 8.49351975531],
            [6.29343661968, 4.80408275125, 4.41739114202, 0.0,
             28.0874527147, 4.13376879093],
            [21.9681438519, 32.7248408842, 23.7695465697,
             28.0874527147, 0.0, 32.2187374711],
            [10.3908084382, 1.58142340946, 8.49351975531,
             4.13376879093, 32.2187374711, 0.0]]

        exp = SymmetricDistanceMatrix(asarray(mtx), dm_lbls)
        obs = self.best._derive_euclidean_dm(cat_mat,
                                             self.bv_dm_88soils.shape[0])
        self.assertEqual(obs.ids, exp.ids)
        assert_almost_equal(obs.data, exp.data)

    def test_call(self):
        """Test the overall functionality of Best."""
        exp = {'method_name': 'BEST',
               'rho_vals': [
                   (0.75, '8'),
                   (0.5, '1,11'),
                   (0.5107142857142857, '1,8,11'),
                   (0.47857142857142854, '1,6,8,11'),
                   (0.46071428571428574, '1,5,6,8,9'),
                   (0.4357142857142857, '1,5,6,8,9,11'),
                   (0.38928571428571423, '1,2,4,5,6,7,8'),
                   (0.38928571428571423, '1,2,4,5,6,7,8,9'),
                   (0.38928571428571423, '1,2,4,5,6,7,8,9,10'),
                   (0.38928571428571423, '1,2,4,5,6,7,8,9,10,11'),
                   (0.16785714285714282, '1,2,3,4,5,6,7,8,9,10,11')],
               'num_vars': 11,
               'vars': ['TOT_ORG_CARB = 1', 'SILT_CLAY = 2',
                        'ELEVATION = 3', 'SOIL_MOISTURE_DEFICIT = 4',
                        'CARB_NITRO_RATIO = 5', 'ANNUAL_SEASON_TEMP = 6',
                        'ANNUAL_SEASON_PRECPT = 7', 'PH = 8',
                        'CMIN_RATE = 9', 'LONGITUDE = 10',
                        'LATITUDE = 11']}

        # the difference caused by the new spearmans_rho function is in the
        # final decimal place of the rho_vals. its not a significant difference,
        # but the test is failing because we are doing a dict comparison for
        # float values.
        obs = self.best()
        self.assertEqual(exp['method_name'], obs['method_name'])
        self.assertEqual(exp['num_vars'], obs['num_vars'])
        self.assertEqual(exp['vars'], obs['vars'])
        for i, j in zip(exp['rho_vals'], obs['rho_vals']):
            assert_almost_equal(i[0], j[0])
            self.assertEqual(i[1], j[1])
        # check that the keys are the same since we have checked all the values
        # we expect to be there are the same
        self.assertTrue(set(exp.keys()) == set(obs.keys()))


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
        dm1 = SymmetricDistanceMatrix(array([[0, 1, 2], [1, 0, 3], [2, 3, 0]]),
                                      ids)
        dm2 = SymmetricDistanceMatrix(array([[0, 2, 5], [2, 0, 8], [5, 8, 0]]),
                                      ids)

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
        m1_dm = SymmetricDistanceMatrix(m1, sample_ids)
        m2_dm = SymmetricDistanceMatrix(m2, sample_ids)

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
            SymmetricDistanceMatrix(array([[0, 1, 4], [1, 0, 3], [4, 3, 0]]),
                                    smpl_ids),
            SymmetricDistanceMatrix(array([[0, 2, 5], [2, 0, 8], [5, 8, 0]]),
                                    smpl_ids),
            SymmetricDistanceMatrix(array([[0, 9, 10], [9, 0, 2], [10, 2, 0]]),
                                    smpl_ids))

        self.small_pm_diff = PartialMantel(
            SymmetricDistanceMatrix(array([[0, 1, 4], [1, 0, 3], [4, 3, 0]]),
                                    smpl_ids),
            SymmetricDistanceMatrix(array([[0, 20, 51], [20, 0, 888],
                                           [51, 888, 0]]), smpl_ids),
            SymmetricDistanceMatrix(array([[0, 9, 10], [9, 0, 2], [10, 2, 0]]),
                                    smpl_ids))

        smpl_ids = ['s1', 's2', 's3', 's4', 's5']
        self.small_pm_diff2 = PartialMantel(
            SymmetricDistanceMatrix(array([[0, 1, 2, 3, 1.4],
                                           [1, 0, 1.5, 1.6, 1.7],
                                           [2, 1.5, 0, 0.8, 1.9],
                                           [3, 1.6, 0.8, 0, 1.0],
                                           [1.4, 1.7, 1.9, 1.0, 0]]),
                                    smpl_ids),
            SymmetricDistanceMatrix(array([[0, 1, 2, 3, 4.1],
                                           [1, 0, 5, 6, 7],
                                           [2, 5, 0, 8, 9],
                                           [3, 6, 8, 0, 10],
                                           [4.1, 7, 9, 10, 0]]), smpl_ids),
            SymmetricDistanceMatrix(array([[0, 1, 2, 3, 4],
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
        self.test_out = get_tmp_filename(tmp_dir=tmp_dir,
                                         prefix='qiime_paired_diff_tests_',
                                         suffix='',
                                         result_constructor=str)
        self.dirs_to_remove.append(self.test_out)
        create_dir(self.test_out)

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
        # expected t values returned, they should be less than (firmicutes) or greater (bacteroidetes) than 2 
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
        table = parse_biom_table(open(biom_table_fp, 'U'))
        self.assertItemsEqual(table.SampleIds, ['subject1', 'subject2'])
        self.assertItemsEqual(table.ObservationIds,
                              ['firmicutes-abundance', 'bacteroidetes-abundance'])
        assert_almost_equal(table
                              [table.getObservationIndex(
                                  'firmicutes-abundance')]
                              [table.getSampleIndex('subject1')],
                              0.1, 2)
        assert_almost_equal(table
                              [table.getObservationIndex(
                                  'bacteroidetes-abundance')]
                              [table.getSampleIndex('subject1')],
                              -0.07, 2)
        assert_almost_equal(table
                              [table.getObservationIndex(
                                  'firmicutes-abundance')]
                              [table.getSampleIndex('subject2')],
                              0.41, 2)
        assert_almost_equal(table
                              [table.getObservationIndex(
                                  'bacteroidetes-abundance')]
                              [table.getSampleIndex('subject2')],
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
        # expected t values returned, they should be less than (firmicutes) or greater (bacteroidetes) than 2
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

if __name__ == "__main__":
    main()
