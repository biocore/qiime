#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"
__status__ = "Development"

"""Test suite for the estimate_observation_richness.py module."""

from biom.parse import parse_biom_table
from biom.table import Table
from cogent.util.unit_test import TestCase, main
from numpy import array

from qiime.estimate_observation_richness import (AbstractFullRichnessEstimator,
        AbstractPointEstimator, Chao1FullRichnessEstimator, EmptySampleError,
        EmptyTableError, MultinomialPointEstimator,
        ObservationRichnessEstimator, RichnessEstimatesResults)

class ObservationRichnessEstimatorTests(TestCase):
    """Tests for the ObservationRichnessEstimator class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        # Single sample, 6 observations, one of which isn't observed in sample.
        self.biom_table1 = parse_biom_table(biom_table_str1)
        self.chao1_estimator = Chao1FullRichnessEstimator()
        self.multi_estimator = MultinomialPointEstimator()

        self.estimator1 = ObservationRichnessEstimator(self.biom_table1,
                                                       self.chao1_estimator,
                                                       self.multi_estimator)

    def test_constructor(self):
        """Test instantiating an AbstractObservationRichnessEstimator."""
        self.assertTrue(isinstance(self.estimator1,
                                   ObservationRichnessEstimator))

    def test_constructor_empty_table(self):
        """Test instantiating an estimator with an empty table."""
        empty_table = Table(array([]), [], [])
        self.assertRaises(EmptyTableError,
                          ObservationRichnessEstimator, empty_table,
                          self.chao1_estimator, self.multi_estimator)

    def test_constructor_empty_sample(self):
        """Test instantiating an estimator with a sample that has no obs."""
        empty_sample_table = parse_biom_table(empty_sample_table_str)
        self.assertRaises(EmptySampleError,
                          ObservationRichnessEstimator, empty_sample_table,
                          self.chao1_estimator, self.multi_estimator)

    def test_getSampleCount(self):
        """Test estimator returns correct number of samples."""
        self.assertEqual(self.estimator1.getSampleCount(), 1)

    def test_getTotalIndividualCounts(self):
        """Returns correct total number of observed individuals per sample."""
        # Verified with iNEXT.
        self.assertEqual(self.estimator1.getTotalIndividualCounts(), [15])

    def test_getObservationCounts(self):
        """Returns correct number of (observed) observations per sample."""
        # Verified with iNEXT.
        self.assertEqual(self.estimator1.getObservationCounts(), [5])

    def test_getAbundanceFrequencyCounts(self):
        """Returns correct abundance frequency counts for each sample."""
        # Verified with iNEXT.
        obs = list(self.estimator1.getAbundanceFrequencyCounts())
        self.assertEqual(obs, [[1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])

    def test_call_interpolate(self):
        """Test __call__ computes correct estimates (interpolation)."""
        # Verified with iNEXT (http://glimmer.rstudio.com/tchsieh/inext/).
        # SE estimates differ because they use a different technique. SE
        # estimates have been verified against values in Colwell 2012 instead
        # (in a separate unit test).

        # Just reference.
        obs = self.estimator1(start=15, stop=15, step_size=1)
        #self.assertFloatEqual(obs.RawEstimatesData,
        #                      [[('S1', 15, 5, 0.674199862463)]])
        self.assertFloatEqual(obs,
                              [[('S1', 15, 5, 0.674199862463)]])

        # start=1 and reference.
        obs = self.estimator1(start=1, stop=1, step_size=1)
        self.assertFloatEqual(obs,
                              [[('S1', 1, 1.0, 0.250252397843),
                                ('S1', 15, 5, 0.674199862463)]])

        # Points in between start=1 and reference.
        obs = self.estimator1(start=1, stop=15, step_size=5)
        self.assertFloatEqual(obs,
                              [[('S1', 1, 1.0, 0.250252397843),
                                ('S1', 6, 3.7382617382617385, 0.676462867498),
                                ('S1', 11, 4.666666666666667, 0.669471144282),
                                ('S1', 15, 5, 0.674199862463)]])

    def test_call_extrapolate(self):
        """Test __call__ computes correct estimates (extrapolation)."""
        # Verified with iNEXT. Differs slightly from their output because
        # they've slightly modified Colwell 2012 equation 9, and we're using
        # the original one.

        #obs = self.estimator1(start=15, stop=30, step_size=15)
        # TODO fix expected variance
        #self.assertFloatEqual(obs.RawEstimatesData,
        #                      [[('S1', 15, 5, 0.674199862463),
        #                        ('S1', 30, 5.4415544562981095, float('inf'))]])

        #obs = self.estimator1(start=20, stop=30, step_size=5)
        # TODO fix expected variance
        #self.assertFloatEqual(obs.RawEstimatesData,
        #                      [[('S1', 15, 5, 0.674199862463),
        #                        ('S1', 20, 5.2555272427983537, float('inf')),
        #                        ('S1', 25, 5.38046614197245, float('inf')),
        #                        ('S1', 30, 5.4415544562981095, float('inf'))]])

    def test_call_full_range(self):
        """Test __call__ computes correct estimates (inter/extrapolation)."""
        # TODO test me!
        pass

    def test_get_points_to_estimate_invalid_input(self):
        """Raises an error on invalid input."""
        # Invalid min.
        self.assertRaises(ValueError, self.estimator1._get_points_to_estimate,
                          0, 10, 1, 5)

        # Invalid step_size.
        self.assertRaises(ValueError, self.estimator1._get_points_to_estimate,
                          1, 10, 0, 5)

        # max < min.
        self.assertRaises(ValueError, self.estimator1._get_points_to_estimate,
                          1, -1, 1, 5)

    def test_get_points_to_estimate(self):
        """Correctly calculates estimation points given range parameters."""
        # Ref in range.
        obs = self.estimator1._get_points_to_estimate(1, 5, 1, 4)
        self.assertEqual(obs, [1, 2, 3, 4, 5])

        # Ref not in range.
        obs = self.estimator1._get_points_to_estimate(5, 10, 2, 4)
        self.assertEqual(obs, [4, 5, 7, 9])


class AbstractFullRichnessEstimatorTests(TestCase):
    """Tests for the AbstractFullRichnessEstimator class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.abstract_estimator = AbstractFullRichnessEstimator()

    def test_estimateFullRichness(self):
        """Test should raise error."""
        self.assertRaises(NotImplementedError,
                          self.abstract_estimator.estimateFullRichness,
                          [1, 2, 3], 6)

    def test_estimateUnobservedObservationCount(self):
        """Test should raise error."""
        self.assertRaises(NotImplementedError,
                self.abstract_estimator.estimateUnobservedObservationCount,
                [1, 2, 3])


class Chao1FullRichnessEstimatorTests(TestCase):
    """Tests for the Chao1FullRichnessEstimator class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.abundance_frequency_counts1 = [1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                                            0, 0, 0]
        self.abundance_frequency_counts2 = [1, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0,
                                            0, 0, 0]
        self.abundance_frequency_counts3 = self.abundance_frequency_counts2[:]
        self.abundance_frequency_counts3[1] = -1

        self.chao1_estimator = Chao1FullRichnessEstimator()

    def test_estimateFullRichness(self):
        """Test returns correct Chao1 full observation richness estimate."""
        # Verified with iNEXT.

        # f2 > 0
        obs = self.chao1_estimator.estimateFullRichness(
                self.abundance_frequency_counts1, 5)
        self.assertFloatEqual(obs, 5.5)

        # f2 == 0
        obs = self.chao1_estimator.estimateFullRichness(
                self.abundance_frequency_counts2, 4)
        self.assertFloatEqual(obs, 4)

        # f2 < 0
        self.assertRaises(ValueError,
                          self.chao1_estimator.estimateFullRichness,
                          self.abundance_frequency_counts3, 4)

    def test_estimateUnobservedObservationCount(self):
        """Test returns correct Chao1 estimate of num unobserved obs."""
        # Verified with iNEXT.
        obs = self.chao1_estimator.estimateUnobservedObservationCount(
                self.abundance_frequency_counts1)
        self.assertFloatEqual(obs, 0.5)


class AbstractPointEstimatorTests(TestCase):
    """Tests for the AbstractPointEstimator class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.abstract_estimator = AbstractPointEstimator()

    def test_estimateExpectedObservationCount(self):
        """Test should raise error."""
        self.assertRaises(NotImplementedError,
                self.abstract_estimator.estimateExpectedObservationCount,
                1, 2, 3, [4], 5)

    def test_estimateExpectedObservationCountStdErr(self):
        """Test should raise error."""
        self.assertRaises(NotImplementedError,
                self.abstract_estimator.estimateExpectedObservationCountStdErr,
                1, 2, 3, [4], 5, 6)


class MultinomialPointEstimatorTests(TestCase):
    """Tests for the MultinomialPointEstimator class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.colwell_fk = colwell_abundance_freq_counts
        self.chao1_estimator = Chao1FullRichnessEstimator()
        self.estimator1 = MultinomialPointEstimator()

    def test_estimateExpectedObservationCount_interpolate(self):
        """Test computing S(m) using data from Colwell 2012 paper."""
        # Verified against results in Colwell 2012 paper.

        # m = 1 (min)
        obs = self.estimator1.estimateExpectedObservationCount(1, 237,
                self.colwell_fk, 112, self.chao1_estimator)
        self.assertFloatEqual(obs, 1.0)

        # m = 20
        obs = self.estimator1.estimateExpectedObservationCount(20, 237,
                self.colwell_fk, 112, self.chao1_estimator)
        self.assertFloatEqual(obs, 15.891665207609165)

        # m = 200
        obs = self.estimator1.estimateExpectedObservationCount(200, 237,
                self.colwell_fk, 112, self.chao1_estimator)
        self.assertFloatEqual(obs, 98.63181822376555)

        # m = 237 (max)
        obs = self.estimator1.estimateExpectedObservationCount(237, 237,
                self.colwell_fk, 112, self.chao1_estimator)
        self.assertFloatEqual(obs, 112.00)

    def test_estimateExpectedObservationCount_extrapolate(self):
        """Test computing S(n+m*) using data from Colwell 2012 paper."""
        # Verified against results in Colwell 2012 paper.

        # m = 237 (min)
        obs = self.estimator1.estimateExpectedObservationCount(237, 237,
                self.colwell_fk, 112, self.chao1_estimator)
        self.assertFloatEqual(obs, 112)

        # m = 337 (n+100)
        obs = self.estimator1.estimateExpectedObservationCount(337, 237,
                self.colwell_fk, 112, self.chao1_estimator)
        self.assertFloatEqual(obs, 145.7369598336187)

        # m = 1237 (n+1000)
        obs = self.estimator1.estimateExpectedObservationCount(1237, 237,
                self.colwell_fk, 112, self.chao1_estimator)
        self.assertFloatEqual(obs, 335.67575295919767)

    def test_estimateExpectedObservationCountStdErr_interpolate(self):
        """Test computing std err of S_m using data from Colwell 2012."""
        # Verified against results in Colwell 2012 paper.

        # m = 1 (min)
        # Note: Colwell 2012 list 0.00 in their table, but after extensive
        # searching I'm not sure why. All other values match theirs, so I'm
        # guessing they're treating 1 as a special case.
        obs = self.estimator1.estimateExpectedObservationCountStdErr(
                1, 237, self.colwell_fk, 112, 1.0, self.chao1_estimator)
        self.assertFloatEqual(obs, 0.20541870170521284)

        # m = 20
        obs = self.estimator1.estimateExpectedObservationCountStdErr(
                20, 237, self.colwell_fk, 112, 15.891665207609165,
                self.chao1_estimator)
        self.assertFloatEqual(obs, 1.9486745986194465)

        # m = 200
        obs = self.estimator1.estimateExpectedObservationCountStdErr(
                200, 237, self.colwell_fk, 112, 98.63181822376555,
                self.chao1_estimator)
        self.assertFloatEqual(obs, 8.147805938386115)

        # m = 237 (max)
        obs = self.estimator1.estimateExpectedObservationCountStdErr(
                237, 237, self.colwell_fk, 112, 112, self.chao1_estimator)
        self.assertFloatEqual(obs, 9.22019783913399)

    def test_estimateExpectedObservationCountStdErr_extrapolate(self):
        """Test computing std err of S(n+m*) using data from Colwell 2012."""
        # Verified against results in Colwell 2012 paper.

        # m = 337 (n+100)
        obs = self.estimator1.estimateExpectedObservationCountStdErr(
                337, 237, self.colwell_fk, 112, 145.7369598336187,
                self.chao1_estimator)
        self.assertFloatEqual(obs, 12.2033650407)

        # m = 437 (n+200)
        obs = self.estimator1.estimateExpectedObservationCountStdErr(
                437, 237, self.colwell_fk, 112, 176.25, self.chao1_estimator)
        self.assertFloatEqual(obs, 15.38155289184887)

        # m = 1237 (n+1000)
        obs = self.estimator1.estimateExpectedObservationCountStdErr(
                1237, 237, self.colwell_fk, 112, 335.68, self.chao1_estimator)
        # Paper shows 48.96, so we're off a little here. Not by a lot, likely
        # just due to rounding differences.
        self.assertFloatEqual(obs, 48.951306831638895)


class RichnessEstimatesResultsTests(TestCase):
    """Tests for the RichnessEstimatesResults class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        #self.res1 = RichnessEstimatesResults()
        pass

    def test_estimateExpectedObservationCount_interpolate(self):
        """Test computing S(m) using data from Colwell 2012 paper."""
        pass


# OTU ID S1 taxonomy
# OTU0   0  foo;bar;baz
# OTU1   1  foo;bar;bazz
# OTU2   2  foo;bar;bazzz
# OTU3   3  foo;bar;bazzzz
# OTU4   4  foo;bar;bazzzzz
# OTU5   5  foo;bar;bazzzzzz
biom_table_str1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "BIOM-Format 1.1.2","date": "2013-04-11T11:39:44.032365","matrix_type": "sparse","matrix_element_type": "float","shape": [6, 1],"data": [[1,0,1.0],[2,0,2.0],[3,0,3.0],[4,0,4.0],[5,0,5.0]],"rows": [{"id": "OTU0", "metadata": {"taxonomy": ["foo", "bar", "baz"]}},{"id": "OTU1", "metadata": {"taxonomy": ["foo", "bar", "bazz"]}},{"id": "OTU2", "metadata": {"taxonomy": ["foo", "bar", "bazzz"]}},{"id": "OTU3", "metadata": {"taxonomy": ["foo", "bar", "bazzzz"]}},{"id": "OTU4", "metadata": {"taxonomy": ["foo", "bar", "bazzzzz"]}},{"id": "OTU5", "metadata": {"taxonomy": ["foo", "bar", "bazzzzzz"]}}],"columns": [{"id": "S1", "metadata": null}]}"""

empty_sample_table_str = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "BIOM-Format 1.1.2","date": "2013-04-11T13:02:56.774981","matrix_type": "dense","matrix_element_type": "float","shape": [1, 1],"data": [[0]],"rows": [{"id": "OTU0", "metadata": null}],"columns": [{"id": "S1", "metadata": null}]}"""

# Taken from Colwell 2012 Osa old growth sample (Table 1b).
colwell_abundance_freq_counts = [84.0, 10.0, 4.0, 3.0, 5.0, 1.0, 2.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 1.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0]


if __name__ == "__main__":
    main()
