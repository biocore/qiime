#!/usr/bin/env python
# File created on 07 Oct 2009.
from __future__ import division

__author__ = "Catherine Lozupone"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Catherine Lozupone", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.8.0"
__maintainer__ = "Catherine Lozupone"
__email__ = "lozupone@colorado.edu"

from collections import defaultdict
from numpy import zeros
from string import strip
from qiime.format import format_otu_table
from biom.table import table_factory

"""Prepare otu table for performing otu significance studies that account
for longitudinal study designs.
"""

def get_sample_individual_info(mapping_data, header, individual_column, \
    timepoint_zero_column):
    """make dicts  mapping:
        samples_from_subject: sample names mapped to samples from the same sub
        sample_to_subtract: sample names mapped to timepoint zero sample names
    """
    samples_from_subject = defaultdict(list)
    sample_to_subtract = {}
    if timepoint_zero_column in header:
        zero_timepoint_index = header.index(timepoint_zero_column)
    else:
        raise ValueError("The category mapping file must have a column called reference_sample. Please see the script documentation for more details.")
    if individual_column in header:
        individual_index = header.index(individual_column)
    else:
        raise ValueError("The category mapping file must have a column called individual. Please see the script documentation for more details.")
    for i in mapping_data:
        sample_name = i[0]
        SUBJECT = i[individual_index]
        for j in mapping_data:
            if j[individual_index] == SUBJECT and j[zero_timepoint_index] == '1':
                sample_to_subtract[sample_name] = j[0]
            if j[individual_index] == SUBJECT:
                samples_from_subject[sample_name].append(j[0])
    return samples_from_subject, sample_to_subtract

def make_new_otu_counts(otu_table, sample_to_subtract, samples_from_subject):
    """make the converted otu table
    """
    new_sample_ids = sample_to_subtract.keys()
    new_sample_ids.sort()
    new_otu_counts = zeros([len(otu_table.ObservationIds),
                            len(new_sample_ids)])
    for index1, otu in enumerate(otu_table.ObservationIds):
        for index2, sample in enumerate(new_sample_ids):
            tpz_sample = sample_to_subtract[sample]
            if tpz_sample in otu_table.SampleIds:
                tpz_sample_index = otu_table.SampleIds.index(tpz_sample)
            else:
                raise ValueError("There are samples in the category mapping file that are not in the otu table, such as sample: " + tpz_sample + ". Removing these samples from the category mapping file will allow you to proceed.")
            #get the new count as the relative abundance of the otu at
            #the later timepoint minus the relative abundance at timepoint zero
            old_sample_index = otu_table.SampleIds.index(sample)
            new_count = otu_table[index1, old_sample_index] - \
                otu_table[index1, tpz_sample_index]
            #make sure that the count is not zero across all of the subject's
            #samples
            has_nonzeros = False
            subject_sample_ids = samples_from_subject[sample]
            for i in subject_sample_ids:
                sample_index = otu_table.SampleIds.index(i)
                if otu_table[index1, sample_index] > 0:
                    has_nonzeros = True
            if has_nonzeros:
                new_otu_counts[index1, index2] = new_count
            else:
                new_otu_counts[index1, index2] = 999999999
    return table_factory(new_otu_counts, new_sample_ids,
                         otu_table.ObservationIds,
                         observation_metadata=otu_table.ObservationMetadata)

def longitudinal_otu_table_conversion_wrapper(otu_table, category_mapping,\
    individual_column, timepoint_zero_column):
    """returns the modified otu_table"""
    otu_table = otu_table.normObservationBySample()
    mapping_data, header, comments = category_mapping
    samples_from_subject, sample_to_subtract = \
        get_sample_individual_info(mapping_data, header, individual_column, \
        timepoint_zero_column)
    return make_new_otu_counts(otu_table, sample_to_subtract,
                               samples_from_subject)
