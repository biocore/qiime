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

"""Contains core data classes used in QIIME."""

from csv import writer

from numpy import array_equal, asarray, copy, ndarray, trace

from qiime.parse import parse_distmat
from qiime.pycogent_backports.test import is_symmetric_and_hollow

class InvalidDistanceMatrixError(Exception):
    pass

class InvalidDistanceMatrixFormatError(Exception):
    pass

class SampleIdMismatchError(Exception):
    pass

class DistanceMatrix(ndarray):
    """

    Implementation is based on numpy subclassing guide:
        http://docs.scipy.org/doc/numpy/user/basics.subclassing.html
    """

    # TODO:
    #
    # - make SampleIds a tuple
    # - make getter/setter properties for SampleIds
    # - check for unique SampleIds
    # - test round-trip read-write-read

    __array_priority__ = -1000.0

    @classmethod
    def fromFile(cls, dm_f):
        """Parses a QIIME distance matrix file into a DistanceMatrix object.

        Arguments:
            dm_f - a file(-like) object containing a QIIME-formatted distance
                matrix
        """
        sids = None
        matrix_data = []

        for line_idx, line in enumerate(dm_f):
            tokens = map(lambda e: e.strip(), line.strip().split('\t'))

            if line_idx == 0:
                # We're at the header (sample IDs).
                sids = tokens
            elif line_idx <= len(sids):
                if tokens[0] == sids[line_idx - 1]:
                    row_data = map(float, tokens[1:])

                    if len(row_data) == len(sids):
                        matrix_data.append(row_data)
                    else:
                        raise InvalidDistanceMatrixFormatError("The number of "
                                "values in row number %d doesn't match the "
                                "number of sample IDs in the header." %
                                line_idx)
                else:
                    raise SampleIdMismatchError("Encountered mismatched "
                            "sample IDs while parsing the distance matrix "
                            "file. Please ensure the sample IDs match between "
                            "the distance matrix header (first row) and the "
                            "row labels (first column).")
            else:
                if ''.join(tokens):
                    # If it isn't a blank line, raise an error because we
                    # shouldn't ignore extra data.
                    raise InvalidDistanceMatrixFormatError("Encountered extra "
                            "rows without corresponding sample IDs in the "
                            "header.")

        return cls(matrix_data, SampleIds=sids)

    @staticmethod
    def _validate_data(data, sids):
        """

        This method is static because it is called from within __new__ (where
        self doesn't exist yet) and also within normal methods.
        """
        if 0 in data.shape:
            raise InvalidDistanceMatrixError("Data must be at least 1x1 in "
                                             "size.")
        elif len(data.shape) != 2:
            raise InvalidDistanceMatrixError("Data must have exactly two "
                                             "dimensions.")
        elif data.shape[0] != data.shape[1]:
            raise InvalidDistanceMatrixError("Data must be square (i.e. have "
                                             "the same number of rows and "
                                             "columns).")
        elif sids is not None and len(sids) != data.shape[0]:
            raise InvalidDistanceMatrixError("The number of sample IDs must "
                                             "match the number of "
                                             "rows/columns in the data.")

    def __new__(cls, input_array, SampleIds=None):
        dm = asarray(input_array)
        cls._validate_data(dm, SampleIds)
        dm = dm.view(cls)

        dm.SampleIds = SampleIds
        dm.flags.writeable = False
        return dm

    def __array_finalize__(self, obj):
        if obj is None:
            return
        else:
            sids = getattr(obj, 'SampleIds', None)
            self._validate_data(obj, sids)

            self.SampleIds = sids
            self.flags.writeable = False

#    def __array_wrap__(self, out_arr, context=None):
#        return super(DistanceMatrix, self).__array_wrap__(out_arr, context)

    def __getslice__(self, start, stop):
        """

        Taken from http://stackoverflow.com/a/14555197 (including docstring):

        This solves a subtle bug, where __getitem__ is not called, and all the
        dimensional checking not done, when a slice of only the first dimension
        is taken, e.g. a[1:3]. From the Python docs:
           Deprecated since version 2.0: Support slice objects as parameters to
           the __getitem__() method. (However, built-in types in CPython
           currently still implement __getslice__(). Therefore, you have to
           override it in derived classes when implementing slicing.)
        """
        return self.__getitem__(slice(start, stop))

    def __getitem__(self, index):
        # We need to override __getitem__ because the __array_finalize__ call
        # is delayed for certain slicing operations, e.g. dm[1:,1:]. Thus, we
        # need to also have a validation check here. See numpy.matrix's
        # __array_finalize__ and __getitem__ implementations for more details.
        out = super(DistanceMatrix, self).__getitem__(index)

        if isinstance(out, ndarray) and self.shape != out.shape:
            out = out.view(ndarray)

        return out

    def __str__(self):
        result = super(DistanceMatrix, self).__str__()

        if self.SampleIds is None:
            result += '\nNo sample IDs'
        else:
            result += '\nSample IDs: %s' % ', '.join(self.SampleIds)

        return result

    def copy(self):
        # We use numpy.copy instead of calling the superclass copy because that
        # one doesn't work with immutable arrays (and changing
        # self.flags.writeable to True temporarily doesn't fix the issue
        # either). numpy.copy returns an ndarray, so we have to view-cast back
        # to DistanceMatrix.
        clone = copy(self).view(DistanceMatrix)

        if self.SampleIds is not None:
            # Deep copy.
            clone.SampleIds = self.SampleIds[:]
        else:
            clone.SampleIds = None

        return clone

    def max(self, axis=None, out=None):
        return self.view(ndarray).max(axis=axis, out=out)

    def min(self, axis=None, out=None):
        return self.view(ndarray).min(axis=axis, out=out)

    def sum(self, axis=None, dtype=None, out=None):
        return self.view(ndarray).sum(axis=axis, dtype=dtype, out=out)

    def all(self, axis=None, out=None):
        return self.view(ndarray).all(axis=axis, out=out)

    @property
    def NumSamples(self):
        """Returns the number of samples (i.e. number of rows or columns)."""
        return self.shape[0]

    def equals(self, other):
        # Use array_equal instead of (a == b).all() because of this issue:
        # http://stackoverflow.com/a/10582030
        # Shape is also checked for us in array_equal.
        if (isinstance(other, self.__class__) and
            self.SampleIds == other.SampleIds and
            array_equal(self, other)):
            return True
        else:
            return False

    def toFile(self, out_f, include_header=True, delimiter='\t'):
        dm_rows = self._format_for_writing(include_header=include_header)
        dm_writer = writer(out_f, delimiter=delimiter, lineterminator='\n')
        dm_writer.writerows(dm_rows)

    def extractTriangle(self, upper=False):
        # Naive implementation...
        result = []

        for col_idx in range(self.shape[0]):
            for row_idx in range(self.shape[0]):
                if upper:
                    if col_idx > row_idx:
                        result.append(self[row_idx][col_idx])
                else:
                    if col_idx < row_idx:
                        result.append(self[row_idx][col_idx])

        return result

    def isSymmetricAndHollow(self):
        """Returns True if the distance matrix is symmetric and hollow."""
        #return (self.T == self).all() and (trace(self) == 0)
        return is_symmetric_and_hollow(self)

    def _format_for_writing(self, include_header=True):
        if include_header and self.SampleIds is not None:
            rows = [[''] + self.SampleIds]

            for sid, dm_row in zip(self.SampleIds, self):
                row = [sid]
                for dm_col in dm_row:
                    row.append(dm_col)

                rows.append(row)
        else:
            rows = [[dm_col for dm_col in dm_row] for dm_row in self]

        return rows
