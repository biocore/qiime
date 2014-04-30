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

import numpy as np
from copy import deepcopy

from skbio.maths.stats.ordination import OrdinationResults
from skbio.maths.gradient import make_groups


TRAJECTORY_ALGORITHMS = ['avg', 'trajectory', 'diff', 'wdiff', None]


def normalize_samples(ord_res, metamap):
    """Creates a copy of `ord_res` including only those samples present in
    metamap and a copy of  `metamap` including only those samples present in
    `ord_res`

    Parameters
    ----------
    ord_res : skbio.maths.stats.ordination.OrdinationResults
        The ordination results
    metamap : pandas.DataFrame
        The metadata map where index are samples ids and columns are the
        metadata categories

    Returns
    -------
    OrdinationResults
        A copy of `ord_res` including only the samples present in `metamap`
    DataFrame
        A copy of `metamap` including only the samples present in `ord_res`

    Raises
    ------
    ValueError
        If `ord_res` and `metamap` does not have samples in common
    """
    ord_sample_ids = set(ord_res.site_ids)
    mm_sample_ids = set(metamap.index)
    sample_ids = ord_sample_ids.intersection(mm_sample_ids)

    if not sample_ids:
        raise ValueError("Ordination results and metadata map had no samples "
                         "in common")

    if ord_sample_ids == sample_ids:
        # No changes needed in ord_res
        n_ord_res = deepcopy(ord_res)
    else:
        # Need to take a subset of ord_res
        boolean_idx = np.in1d(ord_res.site_ids, list(sample_ids))
        idx = np.where(boolean_idx)[0]
        n_ord_res = OrdinationResults(eigvals=ord_res.eigvals,
                                      proportion_explained=
                                      ord_res.proportion_explained,
                                      site=ord_res.site[idx, :],
                                      site_ids=
                                      np.asarray(ord_res.site_ids)[idx])

    if mm_sample_ids == sample_ids:
        # No changes needed in metamap
        n_metamap = deepcopy(metamap)
    else:
        # Need to take a subset of metamap
        n_metamap = metamap.ix[sample_ids]

    return n_ord_res, n_metamap


def run_trajectory_analysis(ord_res, metamap, vector_category,
                            sort_category=None, algorithm=None,
                            vectors_axes=3, weighted=False,
                            window_size=None):
    """
    Parameters
    ----------
    ord_res :
    metamap :
    vector_category :
    sort_category :
    algorithm :
    vectors_axes :
    weighted :
    window_size :

    Returns
    -------
    """
    # Remove any samples from ord_res not present in mapping file
    # And remove any samples from metamap not present in ord_res
    ord_res, metamap = normalize_samples(ord_res, metamap)

    # Create groups
    groups = make_groups(ord_res, metamap, vector_category, sort_category)

    # add_vectors = {}
    # add_vectors['vectors'] = [vector_category, sort_category]
    # add_vectors['weight_by_vector'] = weighted

    # if algorithm:
    #     add_vectors['vectors_output'] = {}
    #     add_vectors['vectors_algorithm'] = algorithm
    #     add_vectors['eigvals'] = coord_data[3]
    #     add_vectors['window_size'] = window_size

    #     # Checks specific for the modified first difference algorithm
    #     if add_vectors['vectors_algorithm'] == 'wdiff':
    #         # Sanity check as the value can only be greater or equal to one
    #         if add_vectors['window_size'] < 1:
    #             option_parser.error("The value of window_size is invalid,the"
    #                                 "value must be greater than zero, not %d"
    #                                 % add_vectors['window_size'])
    # else:
    #     add_vectors['vectors_algorithm'] = None
    # add_vectors['vectors_path'] = output_fp
