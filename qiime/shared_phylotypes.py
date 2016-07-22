#!/usr/bin/env python
# File created on 12 Aug 2010
from __future__ import division

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jens Reeder, Justin Kuczynski", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jose Clemente"
__email__ = "jose.clemente@gmail.com"

"""Computes shared phylotypes between samples"""

from numpy import logical_and, zeros, ones
from qiime.format import format_distance_matrix


def _calc_shared_phylotypes_pairwise(otu_table, i, j):
    """Calculate shared otus between two samples in column i and j.

    otu_table: OTU tables as a OTUtable subclass

    i: a sample id in the OTU table

    j: a sample id in the OTU table
    """
    shared_phylos = logical_and(
        otu_table.data(i, 'sample'),
        otu_table.data(j, 'sample'))

    return shared_phylos.sum()


def _calc_shared_phylotypes_multiple(otu_table, idxs):
    """Calculate shared otus between several samples indexed by values in idxes.

    otu_table: OTU table as a OTUtable subclass
    idxs: list of sample ids in the OTU table
    """

    if len(idxs) < 2:
        raise ValueError("calc_shared_phylotypes_multiple needs at least two "
                         "sampleIDs to comapre")

    shared_phylos = ones(len(otu_table.ids(axis='observation')))

    for id_ in idxs:
        shared_phylos = logical_and(shared_phylos, otu_table.data(id_, 'sample'))

    return shared_phylos.sum()


def calc_shared_phylotypes(otu_table, reference_sample=None):
    """Calculates number of shared phylotypes for each pair of sample.

    infile: otu table filehandle

    reference_sample: if set, will use this sample name to calculate shared
        OTUs between reference sample, and pair of samples. Useful, e.g. when
        the reference sample is the Donor in a transplant study
    """
    if reference_sample:
        ref_idx = reference_sample

    sample_ids = otu_table.ids()
    num_samples = len(sample_ids)
    result_array = zeros((num_samples, num_samples), dtype=int)
    for i, samp1_id in enumerate(sample_ids):
        for j, samp2_id in enumerate(sample_ids[:i + 1]):
            if reference_sample:
                result_array[i, j] = result_array[j, i] = \
                    _calc_shared_phylotypes_multiple(otu_table,
                                                     [samp1_id, samp2_id,
                                                      ref_idx])
            else:
                result_array[i, j] = result_array[j, i] = \
                    _calc_shared_phylotypes_pairwise(otu_table, samp1_id,
                                                     samp2_id)

    return format_distance_matrix(sample_ids, result_array) + "\n"
