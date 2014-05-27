#!/usr/bin/env python
from __future__ import division

__author__ = "Jose Antonio Navas Molina"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jose Antonio Navas Molina", "Antonio Gonzalez Pena",
               "Yoshiki Vazquez Baeza"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jose Antonio Navas Molina"
__email__ = "josenavasmolina@gmail.com"

from skbio.maths.gradient import (AverageVectors, TrajectoryVectors,
                                  DifferenceVectors, WindowDifferenceVectors)


TRAJECTORY_ALGORITHMS = ['avg', 'trajectory', 'diff', 'wdiff']


def run_trajectory_analysis(ord_res, metamap, vector_category,
                            sort_category=None, algorithm='avg',
                            vectors_axes=3, weighted=False,
                            window_size=None):
    """
    Parameters
    ----------
    ord_res : skbio.maths.stats.ordination.OrdinationResults
        The ordination results
    metamap : pandas.DataFrame
        The metadata map where index are samples ids and columns are the
        metadata categories
    vector_category : str
    sort_category : str, optional
    algorithm : str, optional
    vectors_axes : int, optional
    weighted : bool, optional
    window_size : int, optional

    Returns
    -------
    """
    if algorithm == 'avg':
        vecs = AverageVectors(ord_res, metamap, vector_category, sort_category,
                              vectors_axes, weighted)
    elif algorithm == 'trajectory':
        vecs = TrajectoryVectors(ord_res, metamap, vector_category,
                                 sort_category, vectors_axes, weighted)
    elif algorithm == 'diff':
        vecs = DifferenceVectors(ord_res, metamap, vector_category,
                                 sort_category, vectors_axes, weighted)
    elif algorithm == 'wdiff':
        vecs = WindowDifferenceVectors(ord_res, metamap, vector_category,
                                       window_size, sort_category,
                                       vectors_axes, weighted)
    else:
        raise ValueError("Algorithm %s not recognized" % algorithm)

    return vecs
