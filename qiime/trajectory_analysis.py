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

from skbio.maths.stats.ordination import OrdinationResults
from skbio.maths.gradient import make_groups


TRAJECTORY_ALGORITHMS = ['avg', 'trajectory', 'diff', 'wdiff', None]


def filter_unmapped_samples(ord_res, metamap):
    """Creates a copy of `ord_res` including only those samples present in
    metamap

    Parameters
    ----------
    ord_res : skbio.maths.stats.ordination.OrdinationResults
        The ordination results to filter
    metamap : qiime.util.MetadataMap
        The metadata map with the sample ids to include

    Returns
    -------
    OrdinationResults
        A copy of `ord_res` including only the samples present in metamap

    Raises
    ------
    ValueError
        If `ord_res` and `metamap` does not have samples in common
    """
    site_ids = []
    site_idx = []
    mapped_samples = metamap.SampleIds
    for idx, sid in enumerate(ord_res.site_ids):
        if sid in mapped_samples:
            site_ids.append(sid)
            site_idx.append(idx)

    if not site_idx:
        raise ValueError("Ordination results and mapping file had no samples "
                         "in common")

    site = ord_res.site[site_idx, :]

    return OrdinationResults(eigvals=ord_res.eigvals,
                             proportion_explained=ord_res.proportion_explained,
                             site=site,
                             site_ids=site_ids)


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
    # Remove any samples not present in mapping file
    ord_res = filter_unmapped_samples(ord_res, metamap)

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
