#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"
__status__ = "Development"

"""Contains functionality to estimate the observation richness of samples."""

class EmptyTableError(Exception):
    pass

class EmptySampleError(Exception):
    pass

class AbstractObservationRichnessEstimator(object):
    def __init__(self, biom_table):
        if biom_table.isEmpty():
            raise EmptyTableError("The input BIOM table cannot be empty.")
        self._biom_table = biom_table

        for c in self.getTotalIndividualCounts():
            if c < 1:
                raise EmptySampleError("Encountered a sample without any "
                                       "recorded observations.")

    def getSampleCount(self):
        return len(self._biom_table.SampleIds)

    def getTotalIndividualCounts(self):
        return self._biom_table.sum(axis='sample')

    def getObservationCounts(self):
        return [(e > 0).sum(0) for e in self._biom_table.iterSampleData()]

    def getAbundanceFrequencyCounts(self):
        for samp_data, num_individuals in zip(
                self._biom_table.iterSampleData(),
                self.getTotalIndividualCounts()):
            samp_abundance_freq_count = []
            for i in range(1, int(num_individuals + 1)):
                samp_abundance_freq_count.append((samp_data == i).sum(0))
            yield samp_abundance_freq_count

    def getFullRichnessEstimates(self):
        richness_estimates = []

        for num_obs, abundance_freqs in zip(self.getObservationCounts(),
                self.getAbundanceFrequencyCounts()):
            f1 = abundance_freqs[0]
            f2 = abundance_freqs[1]

            if f2 > 0:
                estimated_unobserved_count = f1**2 / (2 * f2)
            elif f2 == 0:
                estimated_unobserved_count = (f1 * (f1 - 1)) / (2 * (f2 + 1))
            else:
                raise ValueError("Encountered a negative f2 value (%d), which "
                                 "is invalid." % f2)

            richness_estimates.append(num_obs + estimated_unobserved_count)

        return richness_estimates
