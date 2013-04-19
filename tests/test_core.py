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

"""Test suite for the core.py module."""

from cogent.util.unit_test import TestCase, main
from numpy import array

from qiime.core import DistanceMatrix

class DistanceMatrixTests(TestCase):
    """Tests for the DistanceMatrix class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.data1 = [[0, 1], [1, 0]]
        self.data2 = [[0, 2], [2, 0]]
        self.sids1 = ['a', 'b']
        self.dm1 = DistanceMatrix(self.data1)
        self.dm2 = DistanceMatrix(self.data1, self.sids1)
        self.dm3 = DistanceMatrix(self.data2, self.sids1)

    def test_constructor(self):
        """Correctly constructs DistanceMatrix instances."""
        self.assertTrue(isinstance(self.dm1, DistanceMatrix))
        self.assertEqual(self.dm1.SampleIds, None)
        self.assertEqual(self.dm1.shape, (2, 2))

    def test_view_casting(self):
        """Correctly constructs DistanceMatrix instances using view casting."""
        # DistanceMatrix -> DistanceMatrix
        dm = self.dm1.view(DistanceMatrix)
        self.assertEqual(dm, self.dm1)
        dm.SampleIds = [1,2,3]
        self.assertEqual(dm, self.dm1)

        # TODO finish me

    def test_equals(self):
        """Correctly identifies instances that are equal (or not)."""
        eq_dm = DistanceMatrix(self.data1)
        self.assertTrue(self.dm1.equals(self.dm1))
        self.assertTrue(self.dm1.equals(eq_dm))
        self.assertTrue(self.dm2.equals(self.dm2))

        self.assertFalse(self.dm1.equals(array(self.data1)))
        self.assertFalse(self.dm1.equals(self.dm2))
        self.assertFalse(self.dm2.equals(self.dm3))


if __name__ == "__main__":
    main()
