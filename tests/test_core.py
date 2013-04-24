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
from numpy import ceil, add, array, matrix, ndarray, power, ones_like

from qiime.core import (DistanceMatrix, InvalidDistanceMatrixError,
                        InvalidDistanceMatrixFormatError,
                        SampleIdMismatchError)

class DistanceMatrixTests(TestCase):
    """Tests for the DistanceMatrix class."""

    def setUp(self):
        """Define some sample data that will be used by the tests."""
        self.data1 = [[42.42]]
        self.data2 = [[0, 1], [1, 0]]
        self.data3 = [[0, 2], [2, 0]]
        self.data4 = [[0, 1], [1.5, 0]]
        self.data5 = [[1, 2, 3], [4, 5, 6], [7, 8, 9]]

        self.sids1 = ['a']
        self.sids2 = ['a', 'b']
        self.sids3 = ['a', 'b', 'c']

        self.data_f1 = '\ta\tb\na\t0\t1\nb\t1.5\t0'.split('\n')
        self.bad_data_f1 = 'a\tb\na\t0\t1\nb\t1'.split('\n')
        self.bad_data_f2 = '\ta\tb\nb\t0\t1\na\t1\t0'.split('\n')
        self.bad_data_f3 = '\ta\tb\na\t0\t1\nb\t1\t0\nfoo'.split('\n')

        self.dm1 = DistanceMatrix(self.data1, self.sids1)
        self.dm2 = DistanceMatrix(self.data2, self.sids2)
        self.dm3 = DistanceMatrix(self.data3, self.sids2)
        self.dm4 = DistanceMatrix(self.data4, self.sids2)
        self.dm5 = DistanceMatrix(self.data5, self.sids3)

    def test_round_trip_read_write(self):
        """Test reading, writing, and reading again works as expected."""
        # Read.
        dm1 = DistanceMatrix.fromFile(self.data_f1)

        # Write.
        f = StringIO()
        dm1.toFile(f)
        f.seek(0)

        # Read.
        dm2 = DistanceMatrix.fromFile(f)
        self.assertTrue(dm1.equals(dm2))

    def test_fromFile(self):
        """Test parsing distance matrix file into a DistanceMatrix instance."""
        # Header with leading tab.
        obs = DistanceMatrix.fromFile(self.data_f1)
        self.assertTrue(self.dm4.equals(obs))

        # Header without leading tab.
        data = self.data_f1[:]
        data[0] = 'a\tb'
        obs = DistanceMatrix.fromFile(data)
        self.assertTrue(self.dm4.equals(obs))

        # Extra newlines at end.
        data = ('\n'.join(self.data_f1) + '\n\n').split('\n')
        obs = DistanceMatrix.fromFile(data)
        self.assertTrue(self.dm4.equals(obs))

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
        self.assertTrue(isinstance(self.dm1, DistanceMatrix))
        self.assertEqual(self.dm1.SampleIds, tuple(self.sids1))
        self.assertEqual(self.dm1.shape, (1, 1))
        self.assertEqual(self.dm1[0][0], 42.42)

    def test_constructor_invalid_input(self):
        """Raises error on invalid distance matrix data."""
        # Empty data.
        self.assertRaises(InvalidDistanceMatrixError, DistanceMatrix, [], [])

        # Invalid number of dimensions.
        self.assertRaises(InvalidDistanceMatrixError, DistanceMatrix,
                          [1, 2, 3], ['a'])

        # Dimensions don't match.
        self.assertRaises(InvalidDistanceMatrixError, DistanceMatrix,
                          [[1, 2, 3]], ['a'])

        # Duplicate sample IDs.
        self.assertRaises(InvalidDistanceMatrixError, DistanceMatrix,
                          self.data2, ['a', 'a'])

        # Number of sample IDs don't match dimensions.
        self.assertRaises(InvalidDistanceMatrixError, DistanceMatrix,
                          self.data2, ['a', 'b', 'c'])

    def test_view_casting(self):
        """Correctly supports view-casting."""
        # DistanceMatrix -> DistanceMatrix
        with self.assertRaises(TypeError):
            _ = self.dm1.view(DistanceMatrix)

        # ndarray -> DistanceMatrix
        with self.assertRaises(TypeError):
            _ = array(self.data1).view(DistanceMatrix)

        # DistanceMatrix -> ndarray
        arr = self.dm1.view(ndarray)
        self.assertTrue(isinstance(arr, ndarray))
        self.assertEqual(arr, array(self.data1))
        # We shouldn't be able to access SampleIds.
        self.assertRaises(AttributeError, getattr, arr, 'SampleIds')

        # DistanceMatrix -> matrix
        mat = self.dm1.view(matrix)
        self.assertTrue(isinstance(mat, matrix))
        self.assertEqual(mat, matrix(self.data1))
        # We shouldn't be able to access SampleIds.
        self.assertRaises(AttributeError, getattr, mat, 'SampleIds')

    def test_ufuncs(self):
        """Test that ufuncs work correctly with DistanceMatrix instances."""
        # DistanceMatrix + ndarray = ndarray
        data = array([[2, 3], [4, 5]])
        exp = array([[2, 4], [5, 5]])
        obs = add(self.dm2, data)
        self.assertEqual(obs, exp)
        self.assertEqual(type(obs), ndarray)

        # DistanceMatrix + DistanceMatrix = ndarray
        dm = DistanceMatrix(data, ['c', 'd'])
        obs = add(dm, self.dm2)
        self.assertEqual(obs, exp)
        self.assertEqual(type(obs), ndarray)

        # Shape mismatch (2x2 + 3x3).
        with self.assertRaises(ValueError):
            _ = add(self.dm2, self.dm5)

    def test_common_ndarray_methods(self):
        """Test that common ndarray methods work as expected with dms."""
        # ndarray.max
        self.assertFloatEqual(self.dm1.max(), 42.42)

        obs = self.dm4.max(axis=0)
        self.assertFloatEqual(obs, array([1.5, 1]))
        self.assertEqual(type(obs), ndarray)

        obs = self.dm4.max(axis=1)
        self.assertFloatEqual(obs, array([1, 1.5]))
        self.assertEqual(type(obs), ndarray)

        # ndarray.min
        self.assertFloatEqual(self.dm1.min(), 42.42)
        
        obs = self.dm4.min(axis=0)
        self.assertEqual(obs, array([0, 0]))
        self.assertEqual(type(obs), ndarray)

        obs = self.dm4.min(axis=1)
        self.assertEqual(obs, array([0, 0]))
        self.assertEqual(type(obs), ndarray)

        # ndarray.sum
        self.assertEqual(self.dm2.sum(), 2)

        obs = self.dm4.sum(axis=0)
        self.assertEqual(obs, array([1.5, 1]))
        self.assertEqual(type(obs), ndarray)

        obs = self.dm4.sum(axis=1)
        self.assertEqual(obs, array([1, 1.5]))
        self.assertEqual(type(obs), ndarray)

        # ndarray.all
        self.assertFalse(self.dm2.all())

        true_dm = DistanceMatrix([[1, 1], [1, 1]], self.sids2)
        self.assertTrue(true_dm.all())
        self.assertEqual(true_dm.all(axis=0), array([True, True]))
        self.assertEqual(true_dm.all(axis=1), array([True, True]))

    def test_getslice(self):
        """Test that __getslice__ defers to __getitem__."""
        # Slice of first dimension only.
        obs = self.dm3[1:]
        self.assertEqual(obs, array([[2, 0]]))
        self.assertEqual(type(obs), ndarray)

    def test_getitem(self):
        """Test __getitem__ delegates to underlying ndarray."""
        # Single element access.
        obs = self.dm2[0,1]
        self.assertEqual(obs, 1)

        # Single element access (via two __getitem__ calls).
        obs = self.dm2[0][1]
        self.assertEqual(obs, 1)

        # Row access.
        obs = self.dm2[1]
        self.assertEqual(obs, array([1, 0]))
        self.assertEqual(type(obs), ndarray)

        # Grab all data.
        obs = self.dm2[0:,0:]
        self.assertEqual(obs, array(self.data2))
        self.assertEqual(type(obs), ndarray)

    def test_mul(self):
        """Test __mul__ delegates to underlying ndarray."""
        obs = self.dm2 * self.dm3
        self.assertEqual(obs, array([[0, 2], [2, 0]]))
        self.assertEqual(type(obs), ndarray)

    def test_str(self):
        """Test getting string representation of DistanceMatrix instances."""
        obs = str(self.dm1)
        self.assertEqual(obs, '[[ 42.42]]\nSample IDs: a')

        obs = str(self.dm4)
        self.assertEqual(obs, '[[ 0.   1. ]\n [ 1.5  0. ]]\nSample IDs: a, b')

    def test_copy(self):
        """Correctly copies DistanceMatrix instances, including SampleIds."""
        # Copies should be equal.
        dm = self.dm1.copy()
        self.assertTrue(self.dm1.equals(dm))

        # After modifying, shouldn't be equal.
        dm[0][0] = 10
        self.assertFalse(self.dm1.equals(dm))

    def test_SampleIds(self):
        """Test getting/setting the sample IDs."""
        # Getter.
        obs = self.dm1.SampleIds
        self.assertEqual(obs, tuple(self.sids1))

        # Setter.
        self.dm1.SampleIds = ['foo']
        obs = self.dm1.SampleIds
        self.assertEqual(obs, tuple(['foo']))

        # Invalid new sample IDs.
        with self.assertRaises(InvalidDistanceMatrixError):
            self.dm1.SampleIds = ['foo', 'bar']
        # Make sure the original object's state hasn't been corrupted.
        obs = self.dm1.SampleIds
        self.assertEqual(obs, tuple(['foo']))

    def test_NumSamples(self):
        """Test getting the number of samples."""
        self.assertEqual(self.dm1.NumSamples, 1)
        self.assertEqual(self.dm2.NumSamples, 2)
        self.assertEqual(self.dm5.NumSamples, 3)

    def test_equals(self):
        """Correctly identifies dm instances that are equal (or not)."""
        self.assertTrue(self.dm2.equals(self.dm2))

        eq_dm = DistanceMatrix(self.data2, self.sids2[:])
        self.assertTrue(self.dm2.equals(eq_dm))
        self.assertTrue(eq_dm.equals(self.dm2))

        # Different class.
        self.assertFalse(self.dm2.equals(array(self.data2)))

        # Different sample IDs.
        ne_dm = DistanceMatrix(self.data2, ['c', 'd'])
        self.assertFalse(self.dm2.equals(ne_dm))

        # Different data.
        self.assertFalse(self.dm2.equals(self.dm3))

    def test_extractTriangle(self):
        """Test extracting upper and lower triangle."""
        # 1x1
        self.assertEqual(self.dm1.extractTriangle(), [])
        self.assertEqual(self.dm1.extractTriangle(upper=True), [])

        # 2x2
        self.assertEqual(self.dm2.extractTriangle(), [1])
        self.assertEqual(self.dm2.extractTriangle(upper=True), [1])

        # 3x3
        self.assertEqual(self.dm5.extractTriangle(), [4, 7, 8])
        self.assertEqual(self.dm5.extractTriangle(upper=True), [2, 3, 6])

    def test_isSymmetricAndHollow(self):
        """Test for symmetry and hollowness of various dms."""
        # 1x1
        self.assertFalse(self.dm1.isSymmetricAndHollow())

        # 2x2
        self.assertTrue(self.dm2.isSymmetricAndHollow())

    def test_toFile(self):
        """Correctly formats and writes distance matrix to file."""
        # Ints.
        f = StringIO()
        self.dm3.toFile(f)
        obs = f.getvalue()
        f.close()
        self.assertEqual(obs, '\ta\tb\na\t0\t2\nb\t2\t0\n')

        # Floats.
        f = StringIO()
        self.dm4.toFile(f)
        obs = f.getvalue()
        f.close()
        self.assertEqual(obs, '\ta\tb\na\t0.0\t1.0\nb\t1.5\t0.0\n')

    def test_validate(self):
        """Empty stub: DistanceMatrix._validate already tested elsewhere."""
        pass

    def test_format_for_writing(self):
        """Correctly formats distance matrix for writing to file."""
        # Without header.
        obs = self.dm2._format_for_writing(include_header=False)
        self.assertEqual(obs, [[0, 1], [1, 0]])

        # With header.
        obs = self.dm2._format_for_writing()
        self.assertEqual(obs, [['', 'a', 'b'], ['a', 0, 1], ['b', 1, 0]])

        # With header, including ints and floats.
        obs = self.dm4._format_for_writing()
        self.assertEqual(obs,
                         [['', 'a', 'b'], ['a', 0.0, 1.0], ['b', 1.5, 0.0]])


if __name__ == "__main__":
    main()
