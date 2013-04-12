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
        AbstractObservationRichnessEstimator, Chao1FullRichnessEstimator,
        EmptySampleError, EmptyTableError, ObservationRichnessInterpolator)

class AbstractObservationRichnessEstimatorTests(TestCase):
    """Tests for the AbstractObservationRichnessEstimator class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        # Single sample, 6 observations, one of which isn't observed in sample.
        self.biom_table1 = parse_biom_table(biom_table_str1)
        self.chao1_estimator = Chao1FullRichnessEstimator()
        self.abstract_estimator1 = AbstractObservationRichnessEstimator(
                self.biom_table1, self.chao1_estimator)

    def test_constructor(self):
        """Test instantiating an AbstractObservationRichnessEstimator."""
        self.assertTrue(isinstance(self.abstract_estimator1,
                                   AbstractObservationRichnessEstimator))

    def test_constructor_empty_table(self):
        """Test instantiating an estimator with an empty table."""
        empty_table = Table(array([]), [], [])
        self.assertRaises(EmptyTableError,
                          AbstractObservationRichnessEstimator, empty_table,
                          self.chao1_estimator)

    def test_constructor_empty_sample(self):
        """Test instantiating an estimator with a sample that has no obs."""
        empty_sample_table = parse_biom_table(empty_sample_table_str)
        self.assertRaises(EmptySampleError,
                          AbstractObservationRichnessEstimator,
                          empty_sample_table, self.chao1_estimator)

    def test_getSampleCount(self):
        """Test estimator returns correct number of samples."""
        self.assertEqual(self.abstract_estimator1.getSampleCount(), 1)

    def test_getTotalIndividualCounts(self):
        """Returns correct total number of observed individuals per sample."""
        self.assertEqual(self.abstract_estimator1.getTotalIndividualCounts(),
                         [15])

    def test_getObservationCounts(self):
        """Returns correct number of (observed) observations per sample."""
        self.assertEqual(self.abstract_estimator1.getObservationCounts(), [5])

    def test_getAbundanceFrequencyCounts(self):
        """Returns correct abundance frequency counts for each sample."""
        obs = list(self.abstract_estimator1.getAbundanceFrequencyCounts())
        self.assertEqual(obs, [[1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0]])

    def test_call(self):
        """Test call raises error."""
        self.assertRaises(NotImplementedError, self.abstract_estimator1)


class ObservationRichnessInterpolatorTests(TestCase):
    """Tests for the ObservationRichnessInterpolator class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.biom_table1 = parse_biom_table(biom_table_str1)
        self.chao1_estimator = Chao1FullRichnessEstimator()
        self.interpolator1 = ObservationRichnessInterpolator(self.biom_table1,
                self.chao1_estimator)
        self.colwell_fk = colwell_abundance_freq_counts

    def test_constructor(self):
        """Test instantiating an ObservationRichnessInterpolator."""
        self.assertTrue(isinstance(self.interpolator1,
                                   AbstractObservationRichnessEstimator))
        self.assertTrue(isinstance(self.interpolator1,
                                   ObservationRichnessInterpolator))

    def test_call(self):
        """Test __call__ computes correct interpolation data."""
        # Verified against iNEXT (http://glimmer.rstudio.com/tchsieh/inext/).
        # SE estimates differ because they use a different technique. SE
        # estimates have been verified against values in Colwell 2012 instead
        # (in a separate unit test).
        obs = self.interpolator1(point_count=1)
        self.assertFloatEqual(obs, [[(15, 5, 0.674199862463)]])

        obs = self.interpolator1(point_count=2)
        self.assertFloatEqual(obs, [[(1, 1.0, 0.250252397843),
                                     (15, 5, 0.674199862463)]])

        obs = self.interpolator1(point_count=4)
        self.assertFloatEqual(obs, [[(1, 1.0, 0.250252397843),
                                     (6, 3.7382617382617385, 0.676462867498),
                                     (11, 4.666666666666667, 0.669471144282),
                                     (15, 5, 0.674199862463)]])

    def test_estimate_expected_observation_count(self):
        """Test computing S_m using data from Colwell 2012 paper."""
        # Verified against results in Colwell 2012 paper.

        # m = 1 (min)
        obs = self.interpolator1._estimate_expected_observation_count(1, 237,
                self.colwell_fk, 112)
        self.assertFloatEqual(obs, 1.0)

        # m = 20
        obs = self.interpolator1._estimate_expected_observation_count(20, 237,
                self.colwell_fk, 112)
        self.assertFloatEqual(obs, 15.891665207609165)

        # m = 200
        obs = self.interpolator1._estimate_expected_observation_count(200, 237,
                self.colwell_fk, 112)
        self.assertFloatEqual(obs, 98.63181822376555)

        # m = 237 (max)
        obs = self.interpolator1._estimate_expected_observation_count(237, 237,
                self.colwell_fk, 112)
        self.assertFloatEqual(obs, 112.00)

    def test_estimate_expected_observation_count_std_err(self):
        """Test computing std err of S_m using data from Colwell 2012."""
        # Verified against results in Colwell 2012 paper.

        # m = 1 (min)
        # Note: Colwell 2012 list 0.00 in their table, but after extensive
        # searching I'm not sure why. All other values match theirs, so I'm
        # guessing they're treating 1 as a special case.
        obs = self.interpolator1._estimate_expected_observation_count_std_err(
                1, 237, self.colwell_fk, 112, 1.0)
        self.assertFloatEqual(obs, 0.20541870170521284)

        # m = 20
        obs = self.interpolator1._estimate_expected_observation_count_std_err(
                20, 237, self.colwell_fk, 112, 15.891665207609165)
        self.assertFloatEqual(obs, 1.9486745986194465)

        # m = 200
        obs = self.interpolator1._estimate_expected_observation_count_std_err(
                200, 237, self.colwell_fk, 112, 98.63181822376555)
        self.assertFloatEqual(obs, 8.147805938386115)

        # m = 237 (max)
        obs = self.interpolator1._estimate_expected_observation_count_std_err(
                237, 237, self.colwell_fk, 112, 112)
        self.assertFloatEqual(obs, 9.22019783913399)


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
        self.chao1_estimator = Chao1FullRichnessEstimator()

    def test_estimateFullRichness(self):
        """Test returns correct Chao1 full observation richness estimate."""
        obs = self.chao1_estimator.estimateFullRichness(
                self.abundance_frequency_counts1, 5)
        self.assertFloatEqual(obs, 5.5)

    def test_estimateUnobservedObservationCount(self):
        """Test returns correct Chao1 estimate of num unobserved obs."""
        obs = self.chao1_estimator.estimateUnobservedObservationCount(
                self.abundance_frequency_counts1)
        self.assertFloatEqual(obs, 0.5)


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
