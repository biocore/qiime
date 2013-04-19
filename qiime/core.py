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

from numpy import array_equal, asarray, copy, ndarray

class InvalidDistanceMatrixError(Exception):
    pass

class DistanceMatrix(ndarray):
    """

    Implementation is based on numpy subclassing guide:
        http://docs.scipy.org/doc/numpy/user/basics.subclassing.html
    """

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

        return clone

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
