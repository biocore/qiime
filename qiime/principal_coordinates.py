#!/usr/bin/env python
import numpy
import sys
from optparse import OptionParser
import cogent.cluster.metric_scaling as ms
from qiime.format import format_coords
from qiime.parse import parse_distmat

__author__ = "Justin Kuzynski"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Rob Knight", "Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Development"


def pcoa(file):
    samples, distmtx = parse_distmat(file)
    # coords, each row is an axis
    coords, eigvals = ms.principal_coordinates_analysis(distmtx)
    
    pcnts = (numpy.abs(eigvals)/sum(numpy.abs(eigvals)))*100
    idxs_descending = pcnts.argsort()[::-1]
    coords = coords[idxs_descending]
    eigvals = eigvals[idxs_descending]
    pcnts = pcnts[idxs_descending]
    
    return format_coords(samples, coords.T, eigvals, pcnts)

