#!/usr/bin/env python
# File created on 08 Jun 2012
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Greg Caporaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from numpy import array
from biom.exception import TableException


def get_filter_to_core_f(table,
                         sample_ids=None,
                         fraction_for_core=1.):
    """ return function that filters a table to its core observations

        table: the biom-format table object to filter
        sample_ids: list of sample ids of interest for the core
         computation (default: all samples are of interest)
        fraction_for_core: the fraction of the sample_ids that
         an observation must have a non-zero count for to be
         considered a core observation

    """
    if not (0. <= fraction_for_core <= 1.):
        raise ValueError(
            "invalid fraction_for_core passed to core filter: %1.2f is outside of range [0,1]." %
            fraction_for_core)
    # generate the position mask, which contains True at SampleIds
    # positions that contain an id in sample_ids
    if sample_ids is None:
        position_mask = array([True] * len(table.ids()))
    else:
        position_mask = array([s in sample_ids for s in table.ids()])
    # determine the number of sample_ids that must have a non-zero
    # value for an OTU to be considered part of the core
    min_count = fraction_for_core * position_mask.sum()

    def f(values, obs_ids, obs_md):
        # count the sample ids with non-zero observation
        # counts that are in sample_ids. if that is greater than
        # the minimum required count, return True
        return ((values != 0) & position_mask).sum() >= min_count
    return f


def filter_table_to_core(table,
                         sample_ids=None,
                         fraction_for_core=1.):
    """ filter a table to it's core observations

        table: the biom-format table object to filter
        sample_ids: list of sample ids of interest for the core
         computation (default: all samples are of interest)
        fraction_for_core: the fraction of the sample_ids that
         an observation must have a non-zero count for to be
         considered a core observation

    """
    filter_f = get_filter_to_core_f(table, sample_ids, fraction_for_core)
    return table.filter(filter_f, axis='observation', inplace=False)


def core_observations_across_sample_ids(table,
                                        sample_ids=None,
                                        fraction_for_core=1.):
    """ get the list of core observations in table

        table: the biom-format table object to filter
        sample_ids: list of sample ids of interest for the core
         computation (default: all samples are of interest)
        fraction_for_core: the fraction of the sample_ids that
         an observation must have a non-zero count for to be
         considered a core observation

    """
    try:
        result = list(filter_table_to_core(
            table, sample_ids, fraction_for_core).ids(axis='observation'))

    except TableException:
        result = []
    return result
