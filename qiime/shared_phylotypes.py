#!/usr/bin/env python
# File created on 12 Aug 2010
from __future__ import division

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Jens Reeder, Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.2.0"
__maintainer__ = "Jose Clemente"
__email__ = "jose.clemente@gmail.com"
__status__ = "Release"
 
"""Computes shared phylotypes between samples"""

from qiime.parse import parse_otu_table
from numpy import logical_and, zeros, ones
from qiime.format import format_distance_matrix

def _calc_shared_phylotypes_pairwise(otu_table,i,j):
    """Calculate shared otus between two samples in column i and j.

    otu_table: OTU tables as 2D-array
    
    i: integer index
    
    j: integer index
    """

    shared_phylos = logical_and(otu_table[:,i], otu_table[:,j])
    
    return shared_phylos.sum()

def _calc_shared_phylotypes_multiple(otu_table, idxs):
    """Calculate shared otus between several samples indexed by values in idxes.

    otu_table: OTU table as 2D-array
    idxs: list of integer indexes
    """

    if len(idxs)< 2:
        raise ValueError, "calc_shared_phylotypes_multiple needs at least two sampleIDs to comapre"
    shared_phylos = ones(len(otu_table[:,1]))

    for idx in idxs:
        shared_phylos = logical_and(shared_phylos, otu_table[:,idx])

    return shared_phylos.sum()

def calc_shared_phylotypes(infile, reference_sample=None):
    """Calculates number of shared phylotypes for each pair of sample.

    infile: otu table filehandle

    reference_sample: if set, will use this sample name to calculate shared OTUs
                      between reference sample, and pair of samples. Useful, 
                      e.g. when the reference sample is the Donor in a transplant study
    """

    sample_ids, otu_ids, otu_table, lineages = parse_otu_table(infile)
 
    if reference_sample:
        ref_idx = sample_ids.index(reference_sample)
    (n,m) = otu_table.shape
    result_array = zeros((m,m), dtype=int)
    for i in range(m):
        for j in range (i+1):
            if reference_sample:
                result_array[i,j] = result_array[j,i] = \
                    _calc_shared_phylotypes_multiple(otu_table, [i, j, ref_idx])
            else:  
                result_array[i,j] = result_array[j,i] = \
                    _calc_shared_phylotypes_pairwise(otu_table, i, j)
                
    return format_distance_matrix(sample_ids, result_array)+"\n"
