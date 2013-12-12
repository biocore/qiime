#!/usr/bin/env python
# File created on 12 Aug 2010
from __future__ import division

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jens Reeder, Justin Kuczynski","Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jose Clemente"
__email__ = "jose.clemente@gmail.com"
 
"""Computes shared phylotypes between samples"""

from biom.parse import parse_biom_table
from numpy import logical_and, zeros, ones
from qiime.format import format_distance_matrix

def _calc_shared_phylotypes_pairwise(otu_table,i,j):
    """Calculate shared otus between two samples in column i and j.

    otu_table: OTU tables as a OTUtable subclass
    
    i: a sample id in the OTU table
    
    j: a sample id in the OTU table
    """
    shared_phylos = logical_and(otu_table.sampleData(i), otu_table.sampleData(j))
    #shared_phylos = logical_and(otu_table[:,i], otu_table[:,j])
    
    return shared_phylos.sum()

def _calc_shared_phylotypes_multiple(otu_table, idxs):
    """Calculate shared otus between several samples indexed by values in idxes.

    otu_table: OTU table as a OTUtable subclass
    idxs: list of sample ids in the OTU table
    """

    if len(idxs)< 2:
        raise ValueError, "calc_shared_phylotypes_multiple needs at least two sampleIDs to comapre"
    #shared_phylos = ones(len(otu_table[:,1]))
    shared_phylos = ones(len(otu_table.ObservationIds))
    #for idx in idxs:
    for id_ in idxs:
        #shared_phylos = logical_and(shared_phylos, otu_table[:,idx])
        shared_phylos = logical_and(shared_phylos, otu_table.sampleData(id_))

    return shared_phylos.sum()

def calc_shared_phylotypes(infile, reference_sample=None):
    """Calculates number of shared phylotypes for each pair of sample.

    infile: otu table filehandle

    reference_sample: if set, will use this sample name to calculate shared OTUs
                      between reference sample, and pair of samples. Useful, 
                      e.g. when the reference sample is the Donor in a transplant study
    """

    otu_table = parse_biom_table(infile)

    if reference_sample:
        #ref_idx = sample_ids.index(reference_sample)
        ref_idx = reference_sample
    
    num_samples = len(otu_table.SampleIds)
    result_array = zeros((num_samples, num_samples), dtype=int)
    for i,samp1_id in enumerate(otu_table.SampleIds):
        for j,samp2_id in enumerate(otu_table.SampleIds[:i+1]):
            if reference_sample:
                result_array[i,j] = result_array[j,i] = \
                    _calc_shared_phylotypes_multiple(otu_table, 
                                                 [samp1_id, samp2_id, ref_idx])
            else:  
                result_array[i,j] = result_array[j,i] = \
                    _calc_shared_phylotypes_pairwise(otu_table, samp1_id, 
                                                      samp2_id)
                
    return format_distance_matrix(otu_table.SampleIds, result_array)+"\n"
