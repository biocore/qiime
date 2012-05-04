#!/usr/bin/env python
from __future__ import division

__author__ = "Michael Dwan"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Michael Dwan", "Logan Knecht", "Jai Ram Rideout",
               "Andrew Cochran"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Michael Dwan"
__email__ = "mdwan.tgen@gmail.com"
__status__ = "Development"

"""Test suite for classes, methods and functions of the stats module."""

from cogent.util.unit_test import TestCase, main
from numpy import array, roll
from numpy.random import permutation

from qiime.stats import (Anosim, CategoryStats, CorrelationStats,
                         DistanceMatrixStats, MantelCorrelogram, Mantel,
                         PartialMantel)
from qiime.util import DistanceMatrix, MetadataMap


class TestHelper(TestCase):
    """Helper class that instantiates some commonly-used objects.

    This class should be subclassed by any test classes that want to use its
    members.
    """

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
        self.overview_dm = DistanceMatrix.parseDistanceMatrix(
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

        # A 1x1 dm.
        self.single_ele_dm = DistanceMatrix(array([[0]]), ['s1'], ['s1'])


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

        # Test constructor as well.
        self.assertRaises(TypeError, DistanceMatrixStats, None)
        self.assertRaises(TypeError, DistanceMatrixStats, 10)
        self.assertRaises(TypeError, DistanceMatrixStats, 20.0)
        self.assertRaises(TypeError, DistanceMatrixStats, "foo")
        self.assertRaises(TypeError, DistanceMatrixStats, {})
        self.assertRaises(TypeError, DistanceMatrixStats, self.overview_dm)
        self.assertRaises(TypeError, DistanceMatrixStats, [1])

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
        mismatch = DistanceMatrix(array([[0]]), ['s2'], ['s2'])

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
        """Must receive a list of strings that are in the mapping file."""
        self.assertRaises(TypeError, setattr, self.cs_overview, 'Categories',
                          "Hello!")
        self.assertRaises(TypeError, setattr, self.cs_overview, 'Categories',
                          self.overview_dm)
        self.assertRaises(ValueError, setattr, self.cs_overview, 'Categories',
                          ["hehehe", 123, "hello"])
        self.assertRaises(ValueError, setattr, self.cs_overview, 'Categories',
                          ["foo"])

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
        self.small_dm = DistanceMatrix.parseDistanceMatrix(self.small_dm_str)

        self.small_dm_tie_str = ["\tsam1\tsam2\tsam3\tsam4",
                                 "sam1\t0\t1\t1\t4",
                                 "sam2\t1\t0\t3\t2",
                                 "sam3\t1\t3\t0\t3",
                                 "sam4\t4\t2\t3\t0"]
        self.small_dm_tie = DistanceMatrix.parseDistanceMatrix(
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
        for samp_id in self.small_dm.SampleIds:
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
        self.assertFloatEqual(obs['r_value'], exp['r_value'])
        self.assertTrue(obs['p_value'] > 0 and obs['p_value'] < 0.06)

    def test_call_small(self):
        """Test __call__() on small dm."""
        # These results were verified with R.
        exp = {'method_name': 'ANOSIM', 'p_value': 0.31, 'r_value': 0.625}
        obs = self.anosim_small()

        self.assertEqual(obs['method_name'], exp['method_name'])
        self.assertFloatEqual(obs['r_value'], exp['r_value'])
        self.assertTrue(obs['p_value'] > 0.28 and obs['p_value'] < 0.42)

    def test_call_small_ties(self):
        """Test __call__() on small dm with ties in ranks."""
        # These results were verified with R.
        exp = {'method_name': 'ANOSIM', 'p_value': 0.67600000000000005,
               'r_value': 0.25}
        obs = self.anosim_small_tie()

        self.assertEqual(obs['method_name'], exp['method_name'])
        self.assertFloatEqual(obs['r_value'], exp['r_value'])
        self.assertTrue(obs['p_value'] > 0.56 and obs['p_value'] < 0.75)

    def test_call_no_perms(self):
        """Test __call__() on small dm with no permutations."""
        # These results were verified with R.
        exp = {'method_name': 'ANOSIM', 'p_value': 'NA', 'r_value': 0.625}
        obs = self.anosim_small(0)

        self.assertEqual(obs['method_name'], exp['method_name'])
        self.assertFloatEqual(obs['r_value'], exp['r_value'])
        self.assertEqual(obs['p_value'], exp['p_value'])

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
        self.assertFloatEqual(obs, exp)

    def test_anosim_small_ties(self):
        """Test _anosim() on small dm with ties."""
        # These results were verified with R.
        exp = 0.25
        obs = self.anosim_small_tie._anosim(self.small_group_map)
        self.assertFloatEqual(obs, exp)

    def test_remove_ties1(self):
        """Test removal of ties. Should return [1.5,1.5]."""
        result = self.anosim_small._remove_ties([1,1],[1,2])
        self.assertEqual(result, [1.5,1.5])

    def test_remove_ties2(self):
        """Should return [3.5,3.5,3.5,3.5,3.5,3.5]."""
        result = self.anosim_small._remove_ties([1,1,1,1,1,1],[1,2,3,4,5,6])
        self.assertEqual(result, [3.5,3.5,3.5,3.5,3.5,3.5])

    def test_remove_ties3(self):
        """Should return [1,3.5,3.5,3.5,3.5,6]."""
        result = self.anosim_small._remove_ties([1,3,3,3,3,8],[1,2,3,4,5,6])
        self.assertEqual(result, [1,3.5,3.5,3.5,3.5,6])

    def test_remove_ties4(self):
        """Should return [1,2,3,4]."""
        result = self.anosim_small._remove_ties([1,2,3,4],[1,2,3,4])
        self.assertEqual(result, [1,2,3,4])

    def test_remove_ties5(self):
        """Should return [1,3,3,3,5.5,5.5,7]."""
        result = self.anosim_small._remove_ties([1,2,2,2,3,3,5],
                                                [1,2,3,4,5,6,7])
        self.assertEqual(result, [1,3,3,3,5.5,5.5,7])

    def test_remove_ties6(self):
        """Should return [1.5,1.5,3.5,3.5]."""
        result = self.anosim_small._remove_ties([1,1,2,2],[1,2,3,4])
        self.assertEqual(result,[1.5,1.5,3.5,3.5])

    def test_get_adjusted_vals(self):
        """Test computing adjusted ranks for ties."""
        exp = [4, 4, 4]
        obs = self.anosim_small._get_adjusted_vals([3, 4, 5], 0, 2)
        self.assertEqual(obs, exp)

    def test_compute_r1(self):
        """Should return .625 for the R statistic on the small dm."""
        sorted_rank = [1.0,2.0,3.5,3.5,5.0,6.0]
        sorted_group = [1.0,0.0,0.0,1.0,0.0,0.0]
        sorted_rank = array(sorted_rank)
        sorted_group = array(sorted_group)
        result = self.anosim_small._compute_r_value(sorted_rank,sorted_group,4)
        self.assertEqual(result, .625)

    def test_anosim_p_test(self):
        """p-value should be .5 for this test."""
        nrs = NonRandomShuffler()
        self.anosim_small.RandomFunction = nrs.permutation

        exp = {'method_name': 'ANOSIM', 'p_value': 0.5, 'r_value': 0.625}
        obs = self.anosim_small(3)

        self.assertEqual(obs['method_name'], exp['method_name'])
        self.assertFloatEqual(obs['r_value'], exp['r_value'])
        self.assertFloatEqual(obs['p_value'], exp['p_value'])


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
        self.small_mc = MantelCorrelogram(
            DistanceMatrix(array([[0, 1, 2], [1, 0, 3], [2, 3, 0]]), ids, ids),
            DistanceMatrix(array([[0, 2, 5], [2, 0, 8], [5, 8, 0]]), ids, ids))

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
        self.assertFloatEqual(obs['class_index'], exp_class_index)

        exp_num_dist = [12, 6, 8, 10, 12, 16, 8]
        self.assertEqual(obs['num_dist'], exp_num_dist)

        exp_mantel_r = [0.73244729118260765, 0.31157641757444593,
            0.17627427296718071, None, None, None, None]
        self.assertFloatEqual(obs['mantel_r'], exp_mantel_r)

        # Test matplotlib Figure for a sane state.
        obs_fig = obs['correlogram_plot']
        obs_ax = obs_fig.get_axes()[0]
        self.assertEqual(obs_ax.get_title(), "Mantel Correlogram")
        self.assertEqual(obs_ax.get_xlabel(), "Distance class index")
        self.assertEqual(obs_ax.get_ylabel(), "Mantel correlation statistic")
        self.assertFloatEqual(obs_ax.get_xticks(), [0.57, 0.58, 0.59, 0.6,
            0.61, 0.62, 0.63, 0.64, 0.65])
        self.assertFloatEqual(obs_ax.get_yticks(), [0.1, 0.2, 0.3, 0.4, 0.5,
            0.6, 0.7, 0.8, 0.9])

        # Test p-values and corrected p-values.
        p_vals = obs['mantel_p']
        corr_p_vals = obs['mantel_p_corr']
        self.assertEqual(len(p_vals), 7)
        self.assertTrue(p_vals[0] >= 0 and p_vals[0] <= 0.01)
        self.assertTrue(p_vals[1] > 0.01 and p_vals[1] <= 0.1)
        self.assertTrue(p_vals[2] > 0.1 and p_vals[2] <= 0.5)
        self.assertEqual(p_vals[3:], [None, None, None, None])
        self.assertFloatEqual(corr_p_vals,
            [p_val * 3 if p_val is not None else None for p_val in p_vals])

    def test_call_small(self):
        """Test running a Mantel correlogram analysis on the smallest input."""
        # The expected output was verified with vegan's mantel correlogram
        # function.
        obs = self.small_mc()

        exp_method_name = 'Mantel Correlogram'
        self.assertEqual(obs['method_name'], exp_method_name)

        exp_class_index = [3.0, 5.0, 7.0]
        self.assertFloatEqual(obs['class_index'], exp_class_index)

        exp_num_dist = [2, 2, 2]
        self.assertEqual(obs['num_dist'], exp_num_dist)

        exp_mantel_r = [0.86602540378443871, None, None]
        self.assertFloatEqual(obs['mantel_r'], exp_mantel_r)

        # Test matplotlib Figure for a sane state.
        obs_fig = obs['correlogram_plot']
        obs_ax = obs_fig.get_axes()[0]
        self.assertEqual(obs_ax.get_title(), "Mantel Correlogram")
        self.assertEqual(obs_ax.get_xlabel(), "Distance class index")
        self.assertEqual(obs_ax.get_ylabel(), "Mantel correlation statistic")
        self.assertFloatEqual(obs_ax.get_xticks(), [2.85, 2.9, 2.95, 3., 3.05,
            3.1, 3.15, 3.2])
        self.assertFloatEqual(obs_ax.get_yticks(), [0.82, 0.83, 0.84, 0.85,
            0.86, 0.87, 0.88, 0.89, 0.9, 0.91])

        # Test p-values and corrected p-values.
        p_vals = obs['mantel_p']
        corr_p_vals = obs['mantel_p_corr']
        self.assertEqual(len(p_vals), 3)
        self.assertTrue(p_vals[0] >= 0 and p_vals[0] <= 0.5)
        self.assertEqual(p_vals[1:], [None, None])
        self.assertFloatEqual(corr_p_vals, p_vals)

    def test_find_distance_classes(self):
        """Test finding the distance classes a matrix's elements are in."""
        exp = (array([[-1,  0,  1], [ 0, -1,  2], [ 1,  2, -1]]),
               [3.0, 5.0, 7.0])
        obs = self.small_mc._find_distance_classes(
            self.small_mc.DistanceMatrices[1], 3)
        self.assertFloatEqual(obs, exp)

        exp = (array([[-1,  1,  2,  0,  0,  5,  7,  4,  6],
            [ 1, -1,  0,  2,  3,  6,  6,  6,  4],
            [ 2,  0, -1,  4,  5,  5,  7,  4,  6],
            [ 0,  2,  4, -1,  3,  3,  3,  3,  2],
            [ 0,  3,  5,  3, -1,  5,  7,  6,  6],
            [ 5,  6,  5,  3,  5, -1,  5,  2,  5],
            [ 7,  6,  7,  3,  7,  5, -1,  0,  0],
            [ 4,  6,  4,  3,  6,  2,  0, -1,  0],
            [ 6,  4,  6,  2,  6,  5,  0,  0, -1]]),
            [0.57381779, 0.60024231, 0.62666684, 0.65309137, 0.67951589,
             0.70594042, 0.73236494, 0.75878947])
        obs = self.mc._find_distance_classes(
            self.mc.DistanceMatrices[1], 8)
        self.assertFloatEqual(obs, exp)

    def test_find_distance_classes_invalid_num_classes(self):
        """Test finding the distance classes for a bad number of classes."""
        self.assertRaises(ValueError, self.mc._find_distance_classes,
                self.mc.DistanceMatrices[1], 0)
        self.assertRaises(ValueError, self.mc._find_distance_classes,
                self.mc.DistanceMatrices[1], -1)

    def test_find_break_points(self):
        """Test finding equal-spaced breakpoints in a range."""
        exp = [-2.2204460492503131e-16, 1.0, 2.0, 3.0, 4.0, 5.0]
        obs = self.mc._find_break_points(0, 5, 5)
        self.assertFloatEqual(obs, exp)

        exp = [-2.0, -1.66666666667, -1.33333333333, -1.0]
        obs = self.mc._find_break_points(-2, -1, 3)
        self.assertFloatEqual(obs, exp)

        exp = [-1.0, -0.5, 0.0, 0.5, 1.0]
        obs = self.mc._find_break_points(-1, 1, 4)
        self.assertFloatEqual(obs, exp)

        exp = [-1.0, 1.0]
        obs = self.mc._find_break_points(-1, 1, 1)
        self.assertFloatEqual(obs, exp)

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
        self.assertFloatEqual(obs, exp)

    def test_correct_p_values_all_None(self):
        """Test p-value correction for all None p-values."""
        exp = [None, None]
        obs = self.mc._correct_p_values([None, None])
        self.assertEqual(obs, exp)

    def test_correct_p_values_mixed(self):
        """Test p-value correction for mixture of None and valid p-values."""
        exp = [None, 0.008, 0.01, None]
        obs = self.mc._correct_p_values([None, 0.004, 0.005, None])
        self.assertFloatEqual(obs, exp)

    def test_correct_p_values_no_change(self):
        """Test p-value correction where none is needed."""
        exp = [None, 0.008]
        obs = self.mc._correct_p_values([None, 0.008])
        self.assertFloatEqual(obs, exp)
        exp = [0.007]
        obs = self.mc._correct_p_values([0.007])
        self.assertFloatEqual(obs, exp)

    def test_correct_p_values_large_correction(self):
        """Test p-value correction that exceeds 1.0."""
        exp = [1, None, 0.03, 0.03]
        obs = self.mc._correct_p_values([0.5, None, 0.01, 0.01])
        self.assertFloatEqual(obs, exp)

    def test_correct_p_values_empty(self):
        """Test p-value correction on empty list."""
        exp = []
        obs = self.mc._correct_p_values([])
        self.assertFloatEqual(obs, exp)

    def test_generate_correlogram(self):
        """Test creating a correlogram plot."""
        obs_fig = self.mc._generate_correlogram([0, 1, 2], [-0.9, 0, 0.9],
                [0.001, 0.1, 0.9])
        obs_ax = obs_fig.get_axes()[0]
        self.assertEqual(obs_ax.get_title(), "Mantel Correlogram")
        self.assertEqual(obs_ax.get_xlabel(), "Distance class index")
        self.assertEqual(obs_ax.get_ylabel(), "Mantel correlation statistic")
        self.assertFloatEqual(obs_ax.get_xticks(), [0., 0.5, 1., 1.5, 2.])
        self.assertFloatEqual(obs_ax.get_yticks(), [-1., -0.5, 0., 0.5, 1.])

    def test_generate_correlogram_empty(self):
        """Test creating a correlogram plot with no data."""
        obs_fig = self.mc._generate_correlogram([], [], [])
        obs_ax = obs_fig.get_axes()[0]
        self.assertEqual(obs_ax.get_title(), "Mantel Correlogram")
        self.assertEqual(obs_ax.get_xlabel(), "Distance class index")
        self.assertEqual(obs_ax.get_ylabel(), "Mantel correlation statistic")
        self.assertFloatEqual(obs_ax.get_xticks(),
            [0., 0.2, 0.4, 0.6, 0.8, 1.0])
        self.assertFloatEqual(obs_ax.get_yticks(),
            [0., 0.2, 0.4, 0.6, 0.8, 1.0])


class MantelTests(TestHelper):
    """Tests for the Mantel class."""

    def setUp(self):
        """Set up Mantel instances for use in tests."""
        super(MantelTests, self).setUp()

        # Used to test that the constuctor sets the default permutations
        # correctly.
        self.defaultPermutations = 999
        self.overview_mantel = Mantel(self.overview_dm, self.overview_dm,
                                      'greater')

        # Create three small test distance matrices. These match the ones used
        # in PyCogent's mantel unit tests.
        m1 = array([[0, 1, 2], [1, 0, 3], [2, 3, 0]])
        m2 = array([[0, 2, 7], [2, 0, 6], [7, 6, 0]])
        m3 = array([[0, 0.5, 0.25], [0.5, 0, 0.1], [0.25, 0.1, 0]])
        sample_ids = ["S1", "S2", "S3"]

        self.m1_dm = DistanceMatrix(m1, sample_ids, sample_ids)
        self.m2_dm = DistanceMatrix(m2, sample_ids, sample_ids)
        self.m3_dm = DistanceMatrix(m3, sample_ids, sample_ids)

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

    def test_call(self):
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

        overview_mantel = Mantel(self.overview_dm, self.overview_dm, 'greater')
        overview_mantel_output = overview_mantel(999)

        obs_method_name = overview_mantel_output['method_name']
        obs_num_permutations = overview_mantel_output['num_perms']
        obs_p_value = overview_mantel_output['p_value']
        obs_r_value = overview_mantel_output['r_value']
        obs_perm_stats_len = len(overview_mantel_output['perm_stats'])
        obs_tail_type = overview_mantel_output['tail_type']

        self.assertEqual(expected_method_name, obs_method_name)
        self.assertTrue(obs_p_value > 0 and obs_p_value < 0.006)
        self.assertFloatEqual(expected_r_value, obs_r_value)
        self.assertFloatEqual(expected_perm_stats_len, obs_perm_stats_len)
        self.assertEqual(expected_number_of_permutations, obs_num_permutations)
        self.assertEqual(expected_tail_type, obs_tail_type)

    # The remaining tests in this class were grabbed from PyCogent's mantel
    # unit tests. They should be removed once we start using PyCogent's version
    # of mantel in a future release. They have only been modified slightly to
    # use the new interface provided in the Mantel class, but the inputs and
    # outputs should all be the same.
    def test_mantel_test_one_sided_greater(self):
        """Test one-sided mantel test (greater)."""
        # This test output was verified by R (their mantel function does a
        # one-sided greater test).
        mantel = Mantel(self.m1_dm, self.m1_dm, 'greater')
        p, stat, perms = mantel._mantel_test(999)

        self.assertTrue(p > 0.09 and p < 0.25)
        self.assertFloatEqual(stat, 1.0)
        self.assertEqual(len(perms), 999)

        mantel = Mantel(self.m1_dm, self.m2_dm, 'greater')
        p, stat, perms = mantel._mantel_test(999)

        self.assertTrue(p > 0.2 and p < 0.5)
        self.assertFloatEqual(stat, 0.755928946018)
        self.assertEqual(len(perms), 999)

    def test_mantel_test_one_sided_less(self):
        """Test one-sided mantel test (less)."""
        # This test output was verified by R (their mantel function does a
        # one-sided greater test, but I modified their output to do a one-sided
        # less test).
        mantel = Mantel(self.m1_dm, self.m1_dm, 'less')
        p, stat, perms = mantel._mantel_test(999)
        self.assertFloatEqual(p, 1.0)
        self.assertFloatEqual(stat, 1.0)
        self.assertEqual(len(perms), 999)

        mantel = Mantel(self.m1_dm, self.m2_dm, 'less')
        p, stat, perms = mantel._mantel_test(999)
        self.assertTrue(p > 0.6 and p < 1.0)
        self.assertFloatEqual(stat, 0.755928946018)
        self.assertEqual(len(perms), 999)

        mantel = Mantel(self.m1_dm, self.m3_dm, 'less')
        p, stat, perms = mantel._mantel_test(999)
        self.assertTrue(p > 0.1 and p < 2.5)
        self.assertFloatEqual(stat, -0.989743318611)
        self.assertEqual(len(perms), 999)

    def test_mantel_test_two_sided(self):
        """Test two-sided mantel test."""
        # This test output was verified by R (their mantel function does a
        # one-sided greater test, but I modified their output to do a two-sided
        # test).
        mantel = Mantel(self.m1_dm, self.m1_dm, 'two sided')
        p, stat, perms = mantel._mantel_test(999)
        self.assertTrue(p > 0.20 and p < 0.45)
        self.assertFloatEqual(stat, 1.0)
        self.assertEqual(len(perms), 999)

        mantel = Mantel(self.m1_dm, self.m2_dm, 'two sided')
        p, stat, perms = mantel._mantel_test(999)
        self.assertTrue(p > 0.6 and p < 0.75)
        self.assertFloatEqual(stat, 0.755928946018)
        self.assertEqual(len(perms), 999)

        mantel = Mantel(self.m1_dm, self.m3_dm, 'two sided')
        p, stat, perms = mantel._mantel_test(999)
        self.assertTrue(p > 0.2 and p < 0.45)
        self.assertFloatEqual(stat, -0.989743318611)
        self.assertEqual(len(perms), 999)


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
        self.small_pm = PartialMantel(DistanceMatrix(array([[1, 3, 2],
            [1, 1, 3], [4, 3, 1]]), smpl_ids, smpl_ids),
            DistanceMatrix(array([[0, 2, 5], [2, 0, 8], [5, 8, 0]]), smpl_ids,
            smpl_ids), DistanceMatrix(array([[10, 7, 13], [9, 7, 0],
            [10, 2, 8]]), smpl_ids, smpl_ids))

        self.small_pm_diff = PartialMantel(DistanceMatrix(array([[1, 3, 2],
            [1, 1, 3], [4, 3, 1]]), smpl_ids, smpl_ids),
            DistanceMatrix(array([[100, 25, 53], [20, 30, 87], [51, 888, 0]]),
            smpl_ids, smpl_ids), DistanceMatrix(array([[10, 7, 13], [9, 7, 0],
            [10, 2, 8]]), smpl_ids, smpl_ids))

        smpl_ids = ['s1', 's2', 's3', 's4', 's5']
        self.small_pm_diff2 = PartialMantel(
            DistanceMatrix(array([[0,1,2,3,1.4],
                                  [1,0,1.5,1.6,1.7],
                                  [2,1.5,0,0.8,1.9],
                                  [3,1.6,0.8,0,1.0],
                                  [1.4,1.7,1.9,1.0,0]]), smpl_ids, smpl_ids),
            DistanceMatrix(array([[0,1,2,3,4.1],
                                  [1,0,5,6,7],
                                  [2,5,0,8,9],
                                  [3,6,8,0,10],
                                  [4.1,7,9,10,0]]), smpl_ids, smpl_ids),
            DistanceMatrix(array([[0,1,2,3,4],
                                  [1,0,5,6,7],
                                  [2,5,0,8,9.1],
                                  [3,6,8,0,10],
                                  [4,7,9.1,10,0]]), smpl_ids, smpl_ids))

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

    def test_call(self):
        """Test running partial Mantel analysis on valid input."""
        obs = self.pm()
        exp_method_name = 'Partial Mantel'
        exp_mantel_r = 0.49999999999999989

        self.assertEqual(obs['method_name'], exp_method_name)
        self.assertFloatEqual(obs['mantel_r'], exp_mantel_r)
        self.assertTrue(obs['mantel_p'] >= 0.001 and obs['mantel_p'] < 0.01)

    def test_call_small(self):
        """Test the running of partial Mantel analysis on small input."""
        obs = self.small_pm()
        exp_method_name = 'Partial Mantel'
        self.assertEqual(obs['method_name'], exp_method_name)

        exp_mantel_r = 0.99999999999999944
        self.assertFloatEqual(obs['mantel_r'], exp_mantel_r)
        self.assertTrue(obs['mantel_p'] > 0.40 and obs['mantel_p'] < 0.60)

        obs = self.small_pm_diff()
        exp_method_name = 'Partial Mantel'
        self.assertEqual(obs['method_name'], exp_method_name)

        exp_mantel_r = 0.99999999999999734
        self.assertFloatEqual(obs['mantel_r'], exp_mantel_r)
        self.assertTrue(obs['mantel_p'] > 0.25 and obs['mantel_p'] < 0.4)

        obs = self.small_pm_diff2()
        exp_method_name = 'Partial Mantel'
        self.assertEqual(obs['method_name'], exp_method_name)

        exp_mantel_r = -0.350624881409
        self.assertFloatEqual(obs['mantel_r'], exp_mantel_r)
        self.assertTrue(obs['mantel_p'] > 0.8)


if __name__ == "__main__":
    main()
