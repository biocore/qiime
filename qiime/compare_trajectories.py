from __future__ import division

__author__ = "Jose Antonio Navas Molina"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jose Antonio Navas Molina", "Antonio Gonzalez Pena",
               "Yoshiki Vazquez Baeza"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jose Antonio Navas Molina"
__email__ = "josenavasmolina@gmail.com"

from itertools import izip

import pandas as pd
from skbio.stats.gradient import (AverageGradientANOVA,
                                   TrajectoryGradientANOVA,
                                   FirstDifferenceGradientANOVA,
                                   WindowDifferenceGradientANOVA)


TRAJECTORY_ALGORITHMS = ['avg', 'trajectory', 'diff', 'wdiff']


def run_trajectory_analysis(ord_res, metadata_map, trajectory_categories=None,
                            sort_category=None, algorithm='avg', axes=3,
                            weighted=False, window_size=None):
    """Executes the volatility analysis using 'algorithm' on the input data

    Parameters
    ----------
    ord_res : skbio.stats.ordination.OrdinationResults
        The ordination results
    metadata_map : pandas.DataFrame
        The metadata map where index are samples ids and columns are the
        metadata categories
    trajectory_categories : list of str, optional
        A list of metadata categories to use to create the trajectories. If
        None is passed, the trajectories for all metadata categories are
        computed. Default: None, compute all of them.
    sort_category : str, optional
        The metadata category to use to sort the trajectories. Default: None,
        sort by sample id.
    algorithm : str, optional
        The algorithm used to perform the volatility analysis. Default: avg.
    axes : int, optional
        The number of axes to account while doing the trajectory specific
        calculations. Pass 0 to compute all of them. Default: 3
    weighted : bool, optional
        If true, the output is weighted by the space between samples in the
        `sort_category` column
    window_size : int or long, optional
        If `algorithm` is 'wdiff', the window size to use while computing the
        difference

    Returns
    -------
    skbio.stats.gradient.GradientANOVAResults
        The results of running the volatility analysis over the input data

    Raises
    ------
    ValueError
        If `algorithm` is not a recognized.
    """
    coords_dict = {sid: coord for sid, coord in izip(ord_res.site_ids,
                                                      ord_res.site)}
    coords = pd.DataFrame.from_dict(coords_dict, orient='index')

    if algorithm == 'avg':
        alg = AverageGradientANOVA(coords, ord_res.proportion_explained,
                                   metadata_map, trajectory_categories,
                                   sort_category, axes, weighted)
    elif algorithm == 'trajectory':
        alg = TrajectoryGradientANOVA(coords, ord_res.proportion_explained,
                                      metadata_map, trajectory_categories,
                                      sort_category, axes, weighted)
    elif algorithm == 'diff':
        alg = FirstDifferenceGradientANOVA(
            coords, ord_res.proportion_explained, metadata_map,
            trajectory_categories, sort_category, axes, weighted)
    elif algorithm == 'wdiff':
        alg = WindowDifferenceGradientANOVA(
            coords, ord_res.proportion_explained, metadata_map, window_size,
            trajectory_categories=trajectory_categories,
            sort_category=sort_category, axes=axes, weighted=weighted)
    else:
        raise ValueError("Algorithm %s not recognized" % algorithm)

    return alg.get_trajectories()
