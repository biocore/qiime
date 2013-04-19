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

from numpy import array_equal, asarray, ndarray

class InvalidDistanceMatrixError(Exception):
    pass

class DistanceMatrix(ndarray):
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
        if SampleIds is not None and len(SampleIds) != dm.shape[0]:
            raise InvalidDistanceMatrixError("The number of sample IDs must "
                                             "match the number of "
                                             "rows/columns.")

        dm = dm.view(cls)
        dm.SampleIds = SampleIds
        dm.flags.writeable = False
        return dm

    def __array_finalize__(self, obj):
        if obj is None:
            return
        else:
            self.SampleIds = getattr(obj, 'SampleIds', None)

    def equals(self, other):
        if (isinstance(other, self.__class__) and
            self.SampleIds == other.SampleIds and
            array_equal(self, other)):
            return True
        else:
            return False
