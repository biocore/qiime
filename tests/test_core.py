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

from StringIO import StringIO

from cogent.util.unit_test import TestCase, main
from numpy import array, ndarray

from qiime.core import (DistanceMatrix, InvalidDistanceMatrixError,
                        InvalidDistanceMatrixFormatError,
                        SampleIdMismatchError)

class DistanceMatrixTests(TestCase):
    """Tests for the DistanceMatrix class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.data1 = [[0, 1], [1, 0]]
        self.data2 = [[0, 2], [2, 0]]
        self.data3 = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]
        self.data4 = [[0, 1], [1.5, 0]]

        self.sids1 = ['a', 'b']
        self.sids2 = ['a', 'b', 'c']

        self.data_f1 = '\ta\tb\na\t0\t1\nb\t1.5\t0'.split('\n')
        self.bad_data_f1 = 'a\tb\na\t0\t1\nb\t1'.split('\n')
        self.bad_data_f2 = '\ta\tb\nb\t0\t1\na\t1\t0'.split('\n')
        self.bad_data_f3 = '\ta\tb\na\t0\t1\nb\t1\t0\nfoo'.split('\n')

        self.dm1 = DistanceMatrix(self.data1)
        self.dm2 = DistanceMatrix(self.data1, self.sids1)
        self.dm3 = DistanceMatrix(self.data2, self.sids1)
        self.dm4 = DistanceMatrix(self.data3, self.sids2)
        self.dm5 = DistanceMatrix(self.data4, self.sids1)

    def test_fromFile(self):
        """Test parsing distance matrix file into a DistanceMatrix instance."""
        # Header with leading tab.
        obs = DistanceMatrix.fromFile(self.data_f1)
        self.assertTrue(self.dm5.equals(obs))

        # Header without leading tab.
        data = self.data_f1[:]
        data[0] = 'a\tb'
        obs = DistanceMatrix.fromFile(data)
        self.assertTrue(self.dm5.equals(obs))

        # Extra newlines at end.
        data = ('\n'.join(self.data_f1) + '\n\n').split('\n')
        obs = DistanceMatrix.fromFile(data)
        self.assertTrue(self.dm5.equals(obs))

    def test_fromFile_invalid_input(self):
        """Raises error on ill-formatted distance matrix file."""
        # Empty dm.
        self.assertRaises(InvalidDistanceMatrixError, DistanceMatrix.fromFile,
                          [])

        # Number of values don't match number of sample IDs.
        self.assertRaises(InvalidDistanceMatrixFormatError,
                          DistanceMatrix.fromFile, self.bad_data_f1)

        # Mismatched sample IDs.
        self.assertRaises(SampleIdMismatchError, DistanceMatrix.fromFile,
                          self.bad_data_f2)

        # Extra data at end.
        self.assertRaises(InvalidDistanceMatrixFormatError,
                          DistanceMatrix.fromFile, self.bad_data_f3)

    def test_constructor(self):
        """Correctly constructs DistanceMatrix instances."""
        # Without sample IDs.
        self.assertTrue(isinstance(self.dm1, DistanceMatrix))
        self.assertEqual(self.dm1.SampleIds, None)
        self.assertEqual(self.dm1.shape, (2, 2))

        # With sample IDs.
        self.assertTrue(isinstance(self.dm2, DistanceMatrix))
        self.assertEqual(self.dm2.SampleIds, self.sids1)
        self.assertEqual(self.dm2.shape, (2, 2))

    def test_constructor_invalid_input(self):
        """Raises error on invalid distance matrix data."""
        # Empty data.
        self.assertRaises(InvalidDistanceMatrixError, DistanceMatrix, [])

        # Invalid number of dimensions.
        self.assertRaises(InvalidDistanceMatrixError, DistanceMatrix,
                          [1, 2, 3])

        # Dimensions don't match.
        self.assertRaises(InvalidDistanceMatrixError, DistanceMatrix,
                          [[1, 2, 3]])

        # Number of sample IDs don't match dimensions.
        self.assertRaises(InvalidDistanceMatrixError, DistanceMatrix,
                          self.data1, ['a', 'b', 'c'])

    def test_view_casting(self):
        """Correctly constructs DistanceMatrix instances using view casting."""
        # DistanceMatrix -> DistanceMatrix
        dm = self.dm1.view(DistanceMatrix)
        self.assertTrue(dm.equals(self.dm1))

        # DistanceMatrix -> ndarray
        arr = self.dm2.view(ndarray)
        self.assertTrue(isinstance(arr, ndarray))
        self.assertEqual(arr.shape, (2, 2))
        # We shouldn't be able to access SampleIds or write to the array.
        self.assertRaises(AttributeError, getattr, arr, 'SampleIds')
        self.assertRaises(RuntimeError, arr.__setitem__, (0, 0), 42)

        # ndarray -> DistanceMatrix
        arr = array(self.data1)
        dm = arr.view(DistanceMatrix)
        self.assertTrue(isinstance(dm, DistanceMatrix))
        self.assertEqual(dm.shape, (2, 2))
        self.assertTrue(dm.SampleIds is None)
        # We shouldn't be able to write to the dm.
        self.assertRaises(RuntimeError, dm.__setitem__, (0, 0), 42)

    def test_view_casting_invalid_input(self):
        """Raises error on invalid distance matrix data."""
        # Empty data.
        arr = array([])
        self.assertRaises(InvalidDistanceMatrixError, arr.view, DistanceMatrix)

        # Invalid number of dimensions.
        arr = array([1, 2, 3])
        self.assertRaises(InvalidDistanceMatrixError, arr.view, DistanceMatrix)

        # Dimensions don't match.
        arr = array([[1, 2, 3]])
        self.assertRaises(InvalidDistanceMatrixError, arr.view, DistanceMatrix)

    def test_new_from_template(self):
        """Correctly constructs DistanceMatrix instances new-from-template."""
        # We're not creating a deep copy, but rather a new instance that points
        # to the old data (including sample IDs).
        dm = self.dm2[:]
        self.assertTrue(dm.equals(self.dm2))
        self.assertRaises(RuntimeError, dm.__setitem__, (0, 0), 42)
        # Not a deep copy.
        self.assertTrue(dm.SampleIds is self.dm2.SampleIds)

    def test_new_from_template_invalid_input(self):
        """Raises error on invalid new-from-template distance matrices."""
        # TODO
        pass

    def test_getslice(self):
        """Test that __getslice__ defers to __getitem__."""
        # Slice of first dimension only.
        obs = self.dm2[1:]
        self.assertEqual(obs, array([[1, 0]]))
        with self.assertRaises(RuntimeError):
            obs[0][1] = 42

    def test_getitem(self):
        """Test __getitem__ works correctly and performs validation checks."""
        # Single element access.
        obs = self.dm1[0,0]
        self.assertEqual(obs, 0)

        # Single element access (via two __getitem__ calls).
        obs = self.dm2[0][0]
        self.assertEqual(obs, 0)

        # Row access. Should get an ndarray that isn't writeable.
        obs = self.dm2[1]
        self.assertEqual(obs, array([1, 0]))
        with self.assertRaises(RuntimeError):
            obs[0] = 42

        # Get the entire distance matrix (should stay DistanceMatrix type).
        obs = self.dm2[0:,0:]
        self.assertTrue(self.dm2.equals(obs))
        with self.assertRaises(RuntimeError):
            obs[0,0] = 42

        #obs = self.dm2[1:,1:]

    def test_copy(self):
        """Correctly copies DistanceMatrix instances, including SampleIds."""
        dm = self.dm2.copy()
        self.assertTrue(dm.equals(self.dm2))
        self.assertRaises(RuntimeError, dm.__setitem__, (0, 0), 42)
        # Check for correct deep copy.
        self.assertFalse(dm.SampleIds is self.dm2.SampleIds)

    def test_equals(self):
        """Correctly identifies instances that are equal (or not)."""
        eq_dm = DistanceMatrix(self.data1)
        self.assertTrue(self.dm1.equals(self.dm1))
        self.assertTrue(self.dm1.equals(eq_dm))
        self.assertTrue(self.dm2.equals(self.dm2))

        self.assertFalse(self.dm1.equals(array(self.data1)))
        self.assertFalse(self.dm1.equals(self.dm2))
        self.assertFalse(self.dm2.equals(self.dm3))

    def test_toFile(self):
        """Correctly formats and writes distance matrix to file."""
        # Ints.
        f = StringIO()
        self.dm2.toFile(f)
        obs = f.getvalue()
        f.close()
        self.assertEqual(obs, '\ta\tb\na\t0\t1\nb\t1\t0\n')

        # Floats.
        f = StringIO()
        self.dm5.toFile(f)
        obs = f.getvalue()
        f.close()
        self.assertEqual(obs, '\ta\tb\na\t0.0\t1.0\nb\t1.5\t0.0\n')

    def test_format(self):
        """Correctly formats distance matrix for writing to file."""
        # Without header.
        obs = self.dm1._format()
        self.assertEqual(obs, [[0, 1], [1, 0]])

        # With header.
        obs = self.dm2._format()
        self.assertEqual(obs, [['', 'a', 'b'], ['a', 0, 1], ['b', 1, 0]])

        # Without header, including ints and floats.
        obs = self.dm5._format()
        self.assertEqual(obs,
                         [['', 'a', 'b'], ['a', 0.0, 1.0], ['b', 1.5, 0.0]])


if __name__ == "__main__":
    main()
