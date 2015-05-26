#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Test suite for the estimate_observation_richness.py module."""

from collections import defaultdict
from StringIO import StringIO

from biom.parse import parse_biom_table
from biom.table import Table
from unittest import TestCase, main
from numpy.testing import assert_almost_equal
from numpy import asarray, array

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
        self.estimator1 = ObservationRichnessEstimator(self.biom_table1,
                                                       Chao1MultinomialPointEstimator)

    def test_constructor(self):
        """Test instantiating an ObservationRichnessEstimator."""
        self.assertTrue(isinstance(self.estimator1,
                                   ObservationRichnessEstimator))

    def test_constructor_empty_table(self):
        """Test instantiating an estimator with an empty table."""
        empty_table = Table(array([]), [], [])
        self.assertRaises(EmptyTableError, ObservationRichnessEstimator,
                          empty_table, Chao1MultinomialPointEstimator)

    def test_getSampleCount(self):
        """Test estimator returns correct number of samples."""
        self.assertEqual(self.estimator1.getSampleCount(), 1)

    def test_call_interpolate(self):
        """Test __call__ computes correct estimates (interpolation)."""
        # Verified with iNEXT (http://glimmer.rstudio.com/tchsieh/inext/).
        # SE estimates differ because they use a different technique. SE
        # estimates have been verified against values in Colwell 2012 instead
        # (in separate unit tests).

        # Just reference.
        obs = self.estimator1(start=15, stop=15, num_steps=1)
        self.assertEqual(obs.getSampleCount(), 1)
        assert_almost_equal(obs.getEstimates('S1'),
                              [(15, 5, 0.674199862463, 3.67859255119, 6.32140744881)])

        # start=1 and reference.
        obs = self.estimator1(start=1, stop=1, num_steps=1)
        self.assertEqual(obs.getSampleCount(), 1)
        assert_almost_equal(obs.getEstimates('S1'),
                              [(1, 1.0, 0.250252397843, 0.509514313183, 1.49048568682),
                               (15, 5, 0.674199862463, 3.67859255119, 6.32140744881)])

        # Points in between start=1 and reference.
        obs = self.estimator1(start=1, stop=15, num_steps=3)
        self.assertEqual(obs.getSampleCount(), 1)
        assert_almost_equal(obs.getEstimates('S1'),
                              [(1, 1.0, 0.250252397843, 0.509514313183, 1.49048568682),
                               (5, 3.40326340326, 0.655024590447, 2.119438797,
                                4.68708800953),
                                  (9, 4.4001998002, 0.680106580075,
                                   3.0672153976, 5.7331842028),
                                  (13, 4.85714285714, 0.665379090563, 3.55302380357,
                                   6.16126191071),
                                  (15, 5, 0.674199862463, 3.67859255119, 6.32140744881)])

    def test_call_extrapolate(self):
        """Test __call__ computes correct estimates (extrapolation)."""
        # Verified with iNEXT. Differs slightly from their output because
        # they've slightly modified Colwell 2012 equation 9, and we're using
        # the original one. SE estimates differ because they use a different
        # technique. SE estimates have been verified against values in Colwell
        # 2012 instead (in separate unit tests).

        obs = self.estimator1(start=15, stop=30, num_steps=1)
        self.assertEqual(obs.getSampleCount(), 1)
        assert_almost_equal(obs.getEstimates('S1'),
                              [(15, 5, 0.674199862463, 3.67859255119, 6.32140744881),
                               (30, 5.4415544562981095, 1.073911829557642, 3.33672594779,
                                7.5463829648)])

        obs = self.estimator1(start=20, stop=30, num_steps=2)
        self.assertEqual(obs.getSampleCount(), 1)
        assert_almost_equal(obs.getEstimates('S1'),
                              [(15, 5, 0.674199862463, 3.67859255119, 6.32140744881),
                               (20, 5.2555272427983537, 0.77331345626875192, 3.73986071975,
                                6.77119376585),
                                  (25, 5.38046614197245, 0.93220670591157662, 3.55337457224,
                                   7.20755771171),
                                  (30, 5.4415544562981095, 1.073911829557642, 3.33672594779,
                                   7.5463829648)])

    def test_get_points_to_estimate_invalid_input(self):
        """Raises an error on invalid input."""
        # Invalid min.
        self.assertRaises(ValueError, self.estimator1._get_points_to_estimate,
                          5, 0, 10, 1)

        # Invalid num_steps.
        self.assertRaises(ValueError, self.estimator1._get_points_to_estimate,
                          5, 1, 10, 0)

        # max < min.
        self.assertRaises(ValueError, self.estimator1._get_points_to_estimate,
                          5, 1, -1, 1)

    def test_get_points_to_estimate(self):
        """Correctly calculates estimation points given range parameters."""
        # Ref in range.
        obs = self.estimator1._get_points_to_estimate(4, 1, 5, 4)
        self.assertEqual(obs, [1, 2, 3, 4, 5])

        # Ref not in range.
        obs = self.estimator1._get_points_to_estimate(4, 5, 10, 2)
        self.assertEqual(obs, [4, 5, 7, 9])

        # stop not supplied.
        obs = self.estimator1._get_points_to_estimate(5, 5, num_steps=2)
        self.assertEqual(obs, [5, 17, 29])


class AbstractPointEstimatorTests(TestCase):

    """Tests for the AbstractPointEstimator class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.colwell_data1 = asarray(colwell_data1)
        self.colwell_data2 = asarray(colwell_data2)

        self.est1 = AbstractPointEstimator(asarray([0, 1, 2, 3, 4, 5]))
        self.est2 = AbstractPointEstimator(self.colwell_data1)
        self.est3 = AbstractPointEstimator(self.colwell_data2)

    def test_constructor(self):
        """Test instantiating an AbstractPointEstimator instance."""
        self.assertTrue(isinstance(self.est1, AbstractPointEstimator))

    def test_constructor_empty_sample(self):
        """Test instantiating an estimator with a sample that has no obs."""
        with self.assertRaises(EmptySampleError):
            _ = AbstractPointEstimator(asarray([0, 0, 0, 0.0, 0, 0.0]))

    def test_getTotalIndividualCount(self):
        """Returns correct total number of observed individuals."""
        # Verified with iNEXT.
        self.assertEqual(self.est1.getTotalIndividualCount(), 15)

        # Verified against results in Colwell 2012 paper.
        self.assertEqual(self.est2.getTotalIndividualCount(), 976)
        self.assertEqual(self.est3.getTotalIndividualCount(), 237)

    def test_getObservationCount(self):
        """Returns correct number of (observed) observations."""
        # Verified with iNEXT.
        self.assertEqual(self.est1.getObservationCount(), 5)

        # Verified against results in Colwell 2012 paper.
        self.assertEqual(self.est2.getObservationCount(), 140)
        self.assertEqual(self.est3.getObservationCount(), 112)

    def test_getAbundanceFrequencyCounts(self):
        """Returns correct abundance frequency counts."""
        # Verified with iNEXT.
        exp = defaultdict(int, {1: 1, 2: 1, 3: 1, 4: 1, 5: 1})
        obs = self.est1.getAbundanceFrequencyCounts()
        self.assertEqual(obs, exp)

        # Verified against results in Colwell 2012 paper.
        self.assertEqual(self.est2.getAbundanceFrequencyCounts(), colwell_fk1)
        self.assertEqual(self.est3.getAbundanceFrequencyCounts(), colwell_fk2)

    def test_call(self):
        """Test should raise error."""
        with self.assertRaises(NotImplementedError):
            self.est1(1)


class Chao1MultinomialPointEstimatorTests(TestCase):

    """Tests for the Chao1MultinomialPointEstimator class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.colwell_data1 = asarray(colwell_data1)
        self.colwell_data2 = asarray(colwell_data2)

        self.samp_data1 = asarray([1, 2, 3, 4, 5])
        self.samp_data2 = asarray([1, 3, 4, 5])

        self.estimator1 = Chao1MultinomialPointEstimator(self.colwell_data1)
        self.estimator2 = Chao1MultinomialPointEstimator(self.colwell_data2)
        self.estimator3 = Chao1MultinomialPointEstimator(self.samp_data1)
        self.estimator4 = Chao1MultinomialPointEstimator(self.samp_data2)

    def test_estimateUnobservedObservationCount(self):
        """Test returns correct Chao1 estimate of num unobserved obs."""
        # Verified with iNEXT.
        obs = self.estimator3.estimateUnobservedObservationCount()
        assert_almost_equal(obs, 0.5)

    def test_estimateFullRichness(self):
        """Test returns correct Chao1 full observation richness estimate."""
        # Verified with iNEXT.

        # f2 > 0
        obs = self.estimator3.estimateFullRichness()
        assert_almost_equal(obs, 5.5)

        # f2 == 0
        obs = self.estimator4.estimateFullRichness()
        assert_almost_equal(obs, 4)

    def test_call_interpolate(self):
        """Test computing S(m) using data from Colwell 2012 paper."""
        # Verified against results in Colwell 2012 paper.

        # Second-growth data.

        # m = 1 (min)
        # Note: Colwell 2012 list the std err as 0.00 in their table, but after
        # extensive searching I'm not sure why. All other values match theirs,
        # so I'm guessing they're treating 1 as a special case (since you can't
        # have an observation count of less than one if you have exactly one
        # individual).
        obs = self.estimator1(1)
        assert_almost_equal(obs, (1.0, 0.17638208235509734, 0.654297471066,
                                    1.34570252893))

        # m = 100
        obs = self.estimator1(100)
        assert_almost_equal(obs, (44.295771605749465, 4.3560838094150975,
                                    35.7580042257, 52.8335389858))

        # m = 800
        obs = self.estimator1(800)
        assert_almost_equal(obs, (126.7974481741264, 7.7007346056227375,
                                    111.704285693, 141.890610656))

        # m = 976 (max)
        obs = self.estimator1(976)
        assert_almost_equal(obs, (140, 8.4270097160038446, 123.483364459,
                                    156.516635541))

        # Old-growth data.

        # m = 1 (min)
        obs = self.estimator2(1)
        assert_almost_equal(obs, (1.0, 0.20541870170521284, 0.597386742907,
                                    1.40261325709))

        # m = 20
        obs = self.estimator2(20)
        assert_almost_equal(obs, (15.891665207609165, 1.9486745986194465,
                                    12.0723331767, 19.7109972385))

        # m = 200
        obs = self.estimator2(200)
        assert_almost_equal(obs, (98.63181822376555, 8.147805938386115,
                                    82.6624120315, 114.601224416))

        # m = 237 (max)
        obs = self.estimator2(237)
        assert_almost_equal(obs, (112.00, 9.22019783913399, 93.928744305,
                                    130.071255695))

    def test_call_extrapolate(self):
        """Test computing S(n+m*) using data from Colwell 2012 paper."""
        # Verified against results in Colwell 2012 paper.

        # Second-growth data.

        # m = 1076 (n+100)
        obs = self.estimator1(1076)
        assert_almost_equal(obs, (146.99829023479796, 8.8700520745653257,
                                    129.613307628, 164.383272842))

        # m = 1176 (n+200)
        obs = self.estimator1(1176)
        assert_almost_equal(obs, (153.6567465407886, 9.3364370482687296,
                                    135.357666182, 171.955826899))

        # m = 1976 (n+1000)
        obs = self.estimator1(1976)
        assert_almost_equal(obs, (196.51177687081162, 13.989113717395064,
                                    169.093617809, 223.929935933))

        # Old-growth data.

        # m = 337 (n+100)
        obs = self.estimator2(337)
        assert_almost_equal(obs, (145.7369598336187, 12.20489285355208,
                                    121.815809405, 169.658110262))

        # m = 437 (n+200)
        obs = self.estimator2(437)
        assert_almost_equal(obs, (176.24777891095846, 15.382655350552035,
                                    146.098328437, 206.397229385))

        # m = 1237 (n+1000)
        obs = self.estimator2(1237)
        assert_almost_equal(obs, (335.67575295919767, 48.962273606327834,
                                    239.71146009, 431.640045829))

    def test_call_invalid_input(self):
        """Test error is raised on invalid input."""
        with self.assertRaises(ValueError):
            self.estimator1(42, confidence_level=0)

    def test_call_na_samples(self):
        """Test on sample without any singletons or doubletons."""
        est = Chao1MultinomialPointEstimator(asarray([4, 3, 4, 5]))
        obs = est(42)
        self.assertEqual(obs, (None, None, None, None))

    def test_partial_derivative_f1(self):
        """Test computes correct partial derivative wrt f1."""
        # Verified with Wolfram Alpha.

        # f2 > 0
        obs = self.estimator1._partial_derivative_f1(2, 3, 10, 42)
        assert_almost_equal(obs, 1.22672908818)

        # f2 == 0
        obs = self.estimator1._partial_derivative_f1(2, 0, 10, 42)
        assert_almost_equal(obs, 1.272173492918482)

        # f1 == 0, f2 == 0
        obs = self.estimator1._partial_derivative_f1(0, 0, 10, 42)
        assert_almost_equal(obs, 1.2961664362634027)

    def test_partial_derivative_f2(self):
        """Test computes correct partial derivative wrt f2."""
        # Verified with Wolfram Alpha.

        # f2 > 0
        obs = self.estimator1._partial_derivative_f2(2, 3, 10, 42)
        assert_almost_equal(obs, 0.9651585982441183)

        # f2 == 0
        obs = self.estimator1._partial_derivative_f2(2, 0, 10, 42)
        assert_almost_equal(obs, 0.9208698803111386)

        # f1 ==0, f2 == 0
        obs = self.estimator1._partial_derivative_f2(0, 0, 10, 42)
        assert_almost_equal(obs, 1.0)


class RichnessEstimatesResultsTests(TestCase):

    """Tests for the RichnessEstimatesResults class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.res1 = RichnessEstimatesResults()

        self.res2 = RichnessEstimatesResults()
        self.res2.addSample('S2', 52)
        self.res2.addSampleEstimate('S2', 1, 3, 0.4, 2.5, 3.5)
        self.res2.addSample('S1', 42)
        self.res2.addSampleEstimate('S1', 10, 20, 2.5, 2.5, 3.5)
        self.res2.addSampleEstimate('S1', 20, 30, 3.5, 2.5, 3.5)
        self.res2.addSampleEstimate('S1', 5, 21, 1.5, 2.5, 3.5)

    def test_constructor(self):
        """Test instantiating a RichnessEstimatesResults instance."""
        self.assertTrue(isinstance(self.res1, RichnessEstimatesResults))
        self.assertEqual(self.res1.getSampleCount(), 0)

    def test_getSampleCount(self):
        """Test getting the number of samples in the results container."""
        self.assertEqual(self.res1.getSampleCount(), 0)

        self.res1.addSample('S1', 42)
        self.assertEqual(self.res1.getSampleCount(), 1)

        self.res1.addSample('S2', 43)
        self.assertEqual(self.res1.getSampleCount(), 2)

    def test_getReferenceIndividualCount(self):
        """Test getting the original number of individuals in a sample."""
        with self.assertRaises(ValueError):
            self.res1.getReferenceIndividualCount('S1')

        self.res1.addSample('S1', 42)
        self.assertEqual(self.res1.getReferenceIndividualCount('S1'), 42)

    def test_getEstimates(self):
        """Test getting the estimates for a sample."""
        with self.assertRaises(ValueError):
            self.res1.getEstimates('S1')

        self.res1.addSample('S1', 42)
        self.res1.addSampleEstimate('S1', 15, 30, 4.75, 2.5, 3.5)
        self.res1.addSampleEstimate('S1', 10, 20, 2.5, 2.5, 3.5)
        assert_almost_equal(self.res1.getEstimates('S1'),
                              [(10, 20, 2.5, 2.5, 3.5), (15, 30, 4.75, 2.5, 3.5)])

    def test_addSample(self):
        """Test adding a new sample to the results container."""
        self.res1.addSample('S1', 42)
        self.assertEqual(self.res1.getSampleCount(), 1)
        self.assertEqual(self.res1.getReferenceIndividualCount('S1'), 42)

        with self.assertRaises(ValueError):
            self.res1.addSample('S1', 45)

    def test_addSampleEstimate(self):
        """Test adding a new estimate for a sample."""
        with self.assertRaises(ValueError):
            self.res1.addSampleEstimate('S1', 10, 20, 2.5, 2.5, 3.5)

        self.res1.addSample('S1', 42)
        self.res1.addSampleEstimate('S1', 10, 20, 2.5, 2.5, 3.5)
        assert_almost_equal(self.res1.getEstimates('S1'),
                              [(10, 20, 2.5, 2.5, 3.5)])

        with self.assertRaises(ValueError):
            self.res1.addSampleEstimate('S1', 10, 35, 0.002, 2.5, 3.5)

    def test_toTable(self):
        """Test writing results container to a table."""
        # Empty results.
        out_f = StringIO()
        self.res1.toTable(out_f)
        self.assertEqual(out_f.getvalue(),
                         "SampleID\tSize\tEstimate\tStd Err\tCI (lower)\tCI (upper)\n")
        out_f.close()

        # Results with multiple samples.
        exp = """SampleID\tSize\tEstimate\tStd Err\tCI (lower)\tCI (upper)
S1\t5\t21\t1.5\t2.5\t3.5
S1\t10\t20\t2.5\t2.5\t3.5
S1\t20\t30\t3.5\t2.5\t3.5
S2\t1\t3\t0.4\t2.5\t3.5
"""
        out_f = StringIO()
        self.res2.toTable(out_f)
        self.assertEqual(out_f.getvalue(), exp)
        out_f.close()

        # Custom header.
        exp = """foo\tbar\tbaz\tbazaar\tbazaaar\tbazaaaar
S1\t5\t21\t1.5\t2.5\t3.5
"""
        out_f = StringIO()
        self.res1.addSample('S1', 42)
        self.res1.addSampleEstimate('S1', 5, 21, 1.5, 2.5, 3.5)
        self.res1.toTable(out_f,
                          header=['foo', 'bar', 'baz', 'bazaar', 'bazaaar', 'bazaaaar'])
        self.assertEqual(out_f.getvalue(), exp)
        out_f.close()

        # Invalid header.
        with self.assertRaises(ValueError):
            out_f = StringIO()
            self.res1.toTable(out_f, header=['foo'])

        # Cells with None as their value.
        exp = """SampleID\tSize\tEstimate\tStd Err\tCI (lower)\tCI (upper)
S1\t43\tN/A\tN/A\tN/A\tN/A
"""
        out_f = StringIO()
        res = RichnessEstimatesResults()
        res.addSample('S1', 42)
        res.addSampleEstimate('S1', 43, None, None, None, None)
        res.toTable(out_f)
        self.assertEqual(out_f.getvalue(), exp)
        out_f.close()


# OTU ID S1 taxonomy
# OTU0   0  foo;bar;baz
# OTU1   1  foo;bar;bazz
# OTU2   2  foo;bar;bazzz
# OTU3   3  foo;bar;bazzzz
# OTU4   4  foo;bar;bazzzzz
# OTU5   5  foo;bar;bazzzzzz
biom_table_str1 = """{"id": "None","format": "Biological Observation Matrix 1.0.0","format_url": "http://biom-format.org","type": "OTU table","generated_by": "BIOM-Format 1.1.2","date": "2013-04-11T11:39:44.032365","matrix_type": "sparse","matrix_element_type": "float","shape": [6, 1],"data": [[1,0,1.0],[2,0,2.0],[3,0,3.0],[4,0,4.0],[5,0,5.0]],"rows": [{"id": "OTU0", "metadata": {"taxonomy": ["foo", "bar", "baz"]}},{"id": "OTU1", "metadata": {"taxonomy": ["foo", "bar", "bazz"]}},{"id": "OTU2", "metadata": {"taxonomy": ["foo", "bar", "bazzz"]}},{"id": "OTU3", "metadata": {"taxonomy": ["foo", "bar", "bazzzz"]}},{"id": "OTU4", "metadata": {"taxonomy": ["foo", "bar", "bazzzzz"]}},{"id": "OTU5", "metadata": {"taxonomy": ["foo", "bar", "bazzzzzz"]}}],"columns": [{"id": "S1", "metadata": null}]}"""

# Taken from Colwell 2012 Osa second growth sample (Table 1a). Added some zeros
# as these should be ignored.
colwell_data1 = [64,
                 1,
                 1,
                 1,
                 1,
                 0.0,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 2,
                 2,
                 2,
                 2,
                 0,
                 2,
                 2,
                 2,
                 2,
                 2,
                 2,
                 2,
                 2,
                 2,
                 2,
                 2,
                 2,
                 2,
                 3,
                 3,
                 3,
                 3,
                 4,
                 4,
                 4,
                 4,
                 4,
                 5,
                 5,
                 5,
                 5,
                 5,
                 6,
                 6,
                 6,
                 6,
                 6,
                 7,
                 7,
                 7,
                 7,
                 7,
                 8,
                 8,
                 8,
                 9,
                 10,
                 10,
                 11,
                 11,
                 11,
                 12,
                 12,
                 77,
                 14,
                 14,
                 17,
                 19,
                 19,
                 20,
                 20,
                 20,
                 21,
                 24,
                 26,
                 40,
                 71,
                 57,
                 57,
                 60,
                 0]
colwell_fk1 = defaultdict(
    int,
    {1: 70,
     2: 17,
     3: 4,
     4: 5,
     5: 5,
     6: 5,
     7: 5,
     8: 3,
     9: 1,
     10: 2,
     11: 3,
     12: 2,
     14: 2,
     17: 1,
     19: 2,
     20: 3,
     21: 1,
     24: 1,
     26: 1,
     40: 1,
     57: 2,
     60: 1,
     64: 1,
     71: 1,
     77: 1})

# Taken from Colwell 2012 Osa old growth sample (Table 1b). Added some zeros as
# these should be ignored.
colwell_data2 = [0,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 0.0,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 0,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 1,
                 2,
                 2,
                 2,
                 2,
                 2,
                 2,
                 2,
                 2,
                 2,
                 2,
                 3,
                 3,
                 3,
                 3,
                 4,
                 4,
                 4,
                 5,
                 5,
                 5,
                 5,
                 5,
                 6,
                 7,
                 7,
                 8,
                 42,
                 14]
colwell_fk2 = defaultdict(
    int,
    {1: 84,
     2: 10,
     3: 4,
     4: 3,
     5: 5,
     6: 1,
     7: 2,
     8: 1,
     14: 1,
     42: 1})


if __name__ == "__main__":
    main()
