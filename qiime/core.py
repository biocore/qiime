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

class DistanceMatrix(object):
    """

    numpy ndarray delegation idea taken from:
        http://stackoverflow.com/a/4759140
    """

    # TODO:
    #
    # - make SampleIds a tuple
    # - make getter/setter properties for SampleIds
    # - check for unique SampleIds
    # - test round-trip read-write-read

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

        return cls(matrix_data, sids)

    def __init__(self, data, SampleIds):
        """"""
        self._data = asarray(data)
        self.SampleIds = SampleIds
        self._validate()

    def __getattr__(self, attr):
        """"""
        return getattr(self._data, attr)

    # Need to forward python "special" methods since we're a new-style class.
    # Will define others as necessary. Details:
    #     http://docs.python.org/2/reference/datamodel.html#new-style-special-lookup

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
        return self._data.__getitem__(index)

    def __mul__(self, other):
        return self._data.__mul__(other)

    # End "special" method forwarding.

    def __str__(self):
        """"""
        return str(self._data) + '\nSample IDs: %s' % ', '.join(self.SampleIds)

    def copy(self):
        """"""
        return DistanceMatrix(self._data.copy(), self.SampleIds[:])

    @property
    def NumSamples(self):
        """Returns the number of samples (i.e. number of rows or columns)."""
        return len(self.SampleIds)

    def equals(self, other):
        """"""
        # Use array_equal instead of (a == b).all() because of this issue:
        # http://stackoverflow.com/a/10582030
        # Shape is also checked for us in array_equal.
        if (isinstance(other, self.__class__) and
            self.SampleIds == other.SampleIds and
            array_equal(self, other)):
            return True
        else:
            return False

    def extractTriangle(self, upper=False):
        """"""
        # Naive implementation... need to fix.
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
        return is_symmetric_and_hollow(self)

    def toFile(self, out_f, include_header=True, delimiter='\t'):
        """"""
        dm_rows = self._format_for_writing(include_header=include_header)
        dm_writer = writer(out_f, delimiter=delimiter, lineterminator='\n')
        dm_writer.writerows(dm_rows)

    def _validate(self):
        """"""
        if 0 in self.shape:
            raise InvalidDistanceMatrixError("Data must be at least 1x1 in "
                                             "size.")
        elif len(self.shape) != 2:
            raise InvalidDistanceMatrixError("Data must have exactly two "
                                             "dimensions.")
        elif self.shape[0] != self.shape[1]:
            raise InvalidDistanceMatrixError("Data must be square (i.e. have "
                                             "the same number of rows and "
                                             "columns).")
        elif len(self.SampleIds) != self.shape[0]:
            raise InvalidDistanceMatrixError("The number of sample IDs must "
                                             "match the number of "
                                             "rows/columns in the data.")

    def _format_for_writing(self, include_header=True):
        """"""
        if include_header:
            rows = [[''] + self.SampleIds]

            for sid, dm_row in zip(self.SampleIds, self):
                row = [sid]
                for dm_col in dm_row:
                    row.append(dm_col)
                rows.append(row)
        else:
            rows = [[dm_col for dm_col in dm_row] for dm_row in self]

        return rows
