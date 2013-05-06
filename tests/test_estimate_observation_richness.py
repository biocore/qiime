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

from collections import defaultdict

from biom.parse import parse_biom_table
from biom.table import Table

from cogent.util.unit_test import TestCase, main

from numpy import array

from qiime.estimate_observation_richness import (AbstractPointEstimator,
        Chao1MultinomialPointEstimator, EmptySampleError,
        EmptyTableError, ObservationRichnessEstimator,
        RichnessEstimatesResults)

class ObservationRichnessEstimatorTests(TestCase):
    """Tests for the ObservationRichnessEstimator class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        # Single sample, 6 observations, one of which isn't observed in sample.
        self.biom_table1 = parse_biom_table(biom_table_str1)
        self.multi_estimator = Chao1MultinomialPointEstimator()

        self.estimator1 = ObservationRichnessEstimator(self.multi_estimator,
                                                       self.biom_table1)

    def test_constructor(self):
        """Test instantiating an ObservationRichnessEstimator."""
        self.assertTrue(isinstance(self.estimator1,
                                   ObservationRichnessEstimator))

    def test_constructor_empty_table(self):
        """Test instantiating an estimator with an empty table."""
        empty_table = Table(array([]), [], [])
        self.assertRaises(EmptyTableError,
                          ObservationRichnessEstimator, self.multi_estimator,
                          empty_table)

    def test_constructor_empty_sample(self):
        """Test instantiating an estimator with a sample that has no obs."""
        empty_sample_table = parse_biom_table(empty_sample_table_str)
        self.assertRaises(EmptySampleError,
                          ObservationRichnessEstimator, self.multi_estimator,
                          empty_sample_table)

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
        exp = [defaultdict(int, {1: 1, 2: 1, 3: 1, 4: 1, 5: 1})]
        obs = list(self.estimator1.getAbundanceFrequencyCounts())
        self.assertEqual(obs, exp)

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


class AbstractPointEstimatorTests(TestCase):
    """Tests for the AbstractPointEstimator class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.abstract_estimator = AbstractPointEstimator()

    def test_call(self):
        """Test should raise error."""
        with self.assertRaises(NotImplementedError):
            self.abstract_estimator(1, 2, 3, [4])


class Chao1MultinomialPointEstimatorTests(TestCase):
    """Tests for the Chao1MultinomialPointEstimator class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.estimator1 = Chao1MultinomialPointEstimator()
        self.colwell_fk = colwell_abundance_freq_counts
        self.abundance_fk1 = defaultdict(int, {1: 1, 2: 1, 3: 1, 4: 1, 5: 1})
        self.abundance_fk2 = defaultdict(int, {1: 1, 3: 1, 4: 1, 5: 1})
        self.abundance_fk3 = self.abundance_fk2.copy()
        self.abundance_fk3[2] = -1

    def test_call_interpolate(self):
        """Test computing S(m) using data from Colwell 2012 paper."""
        # Verified against results in Colwell 2012 paper.

        # m = 1 (min)
        # Note: Colwell 2012 list the std err as 0.00 in their table, but after
        # extensive searching I'm not sure why. All other values match theirs,
        # so I'm guessing they're treating 1 as a special case (since you can't
        # have an observation count of less than one if you have exactly one
        # individual).
        obs = self.estimator1(1, 237, self.colwell_fk, 112)
        self.assertFloatEqual(obs, (1.0, 0.20541870170521284))

        # m = 20
        obs = self.estimator1(20, 237, self.colwell_fk, 112)
        self.assertFloatEqual(obs, (15.891665207609165, 1.9486745986194465))

        # m = 200
        obs = self.estimator1(200, 237, self.colwell_fk, 112)
        self.assertFloatEqual(obs, (98.63181822376555, 8.147805938386115))

        # m = 237 (max)
        obs = self.estimator1(237, 237, self.colwell_fk, 112)
        self.assertFloatEqual(obs, (112.00, 9.22019783913399))

    def test_call_extrapolate(self):
        """Test computing S(n+m*) using data from Colwell 2012 paper."""
        # Verified against results in Colwell 2012 paper.

        # m = 337 (n+100)
        obs = self.estimator1(337, 237, self.colwell_fk, 112)
        self.assertFloatEqual(obs, (145.7369598336187, 12.2033650407))

        # m = 437 (n+200)
        obs = self.estimator1(437, 237, self.colwell_fk, 112)
        self.assertFloatEqual(obs, (176.24777891095846, 15.38155289184887))

        # m = 1237 (n+1000)
        # Paper shows the std err as 48.96, so we're off a little here. Not by
        # a lot, likely just due to rounding differences.
        obs = self.estimator1(1237, 237, self.colwell_fk, 112)
        self.assertFloatEqual(obs, (335.67575295919767, 48.951306831638895))

    def test_estimateFullRichness(self):
        """Test returns correct Chao1 full observation richness estimate."""
        # Verified with iNEXT.

        # f2 > 0
        obs = self.estimator1.estimateFullRichness(self.abundance_fk1, 5)
        self.assertFloatEqual(obs, 5.5)

        # f2 == 0
        obs = self.estimator1.estimateFullRichness(self.abundance_fk2, 4)
        self.assertFloatEqual(obs, 4)

        # f2 < 0
        self.assertRaises(ValueError, self.estimator1.estimateFullRichness,
                          self.abundance_fk3, 4)

    def test_estimateUnobservedObservationCount(self):
        """Test returns correct Chao1 estimate of num unobserved obs."""
        # Verified with iNEXT.
        obs = self.estimator1.estimateUnobservedObservationCount(
                self.abundance_fk1)
        self.assertFloatEqual(obs, 0.5)


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
colwell_abundance_freq_counts = defaultdict(int, {1: 84, 2: 10, 3: 4, 4: 3, 5: 5, 6: 1, 7: 2, 8: 1, 14: 1, 42: 1})


if __name__ == "__main__":
    main()
