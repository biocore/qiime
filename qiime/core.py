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

    def __new__(cls, input_array, SampleIds=None):
        dm = asarray(input_array)

        if 0 in dm.shape:
            raise InvalidDistanceMatrixError("The input array must be at "
                                             "least 1x1 in size.")
        elif len(dm.shape) != 2:
            raise InvalidDistanceMatrixError("The input array must have "
                                             "exactly two dimensions.")
        elif dm.shape[0] != dm.shape[1]:
            raise InvalidDistanceMatrixError("The input array must be square "
                                             "(i.e. have the same number of "
                                             "rows and columns).")
        elif SampleIds is not None and len(SampleIds) != dm.shape[0]:
            raise InvalidDistanceMatrixError("The number of sample IDs must "
                                             "match the number of "
                                             "rows/columns in the distance "
                                             "matrix.")

        dm = dm.view(cls)
        dm.SampleIds = SampleIds
        dm.flags.writeable = False
        return dm

    def __array_finalize__(self, obj):
        if obj is None:
            return
        else:
            if 0 in obj.shape:
                raise InvalidDistanceMatrixError("The input array must be at "
                                                 "least 1x1 in size.")
            elif len(obj.shape) != 2:
                raise InvalidDistanceMatrixError("The input array must have "
                                                 "exactly two dimensions.")
            elif obj.shape[0] != obj.shape[1]:
                raise InvalidDistanceMatrixError("The input array must be square "
                                                 "(i.e. have the same number of "
                                                 "rows and columns).")

            sids = getattr(obj, 'SampleIds', None)
            if sids is not None and len(sids) != obj.shape[0]:
                raise InvalidDistanceMatrixError("The number of sample IDs must "
                                                 "match the number of "
                                                 "rows/columns in the distance "
                                                 "matrix.")

            self.SampleIds = sids
            self.flags.writeable = False

    def copy(self):
        clone = copy(self).view(DistanceMatrix)

        if self.SampleIds is not None:
            clone.SampleIds = self.SampleIds[:]

        return clone

    def equals(self, other):
        # Use array_equal instead of (a == b).all() because of this issue:
        # http://stackoverflow.com/a/10582030
        if (isinstance(other, self.__class__) and
            self.SampleIds == other.SampleIds and
            array_equal(self, other)):
            return True
        else:
            return False
