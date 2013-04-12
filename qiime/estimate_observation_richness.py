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

from math import factorial
from numpy import ceil, sqrt

class EmptyTableError(Exception):
    pass

class EmptySampleError(Exception):
    pass

class ObservationRichnessEstimator(object):
    def __init__(self, biom_table, FullRichnessEstimator, PointEstimator):
        if biom_table.isEmpty():
            raise EmptyTableError("The input BIOM table cannot be empty.")
        self._biom_table = biom_table

        for c in self.getTotalIndividualCounts():
            if c < 1:
                raise EmptySampleError("Encountered a sample without any "
                                       "recorded observations.")

        self.FullRichnessEstimator = FullRichnessEstimator
        self.PointEstimator = PointEstimator

    def getSampleCount(self):
        return len(self._biom_table.SampleIds)

    def getTotalIndividualCounts(self):
        return [int(e) for e in self._biom_table.sum(axis='sample')]

    def getObservationCounts(self):
        return [(e > 0).sum(0) for e in self._biom_table.iterSampleData()]

    def getAbundanceFrequencyCounts(self):
        for samp_data, num_individuals in zip(
                self._biom_table.iterSampleData(),
                self.getTotalIndividualCounts()):
            samp_abundance_freq_count = []

            for i in range(1, num_individuals + 1):
                samp_abundance_freq_count.append((samp_data == i).sum(0))
            yield samp_abundance_freq_count

    def __call__(self, point_count=40):
        per_sample_results = []

        for samp_data, num_obs, n, abundance_freqs in zip(
                self._biom_table.iterSampleData(),
                self.getObservationCounts(),
                self.getTotalIndividualCounts(),
                self.getAbundanceFrequencyCounts()):
            samp_data = samp_data[samp_data > 0]

            if point_count > 1:
                step_size = int(ceil(n / (point_count - 1)))
                sizes = range(1, n, step_size)
            else:
                sizes = []
            sizes.append(n)

            # size <= n
            size_results = []
            for size in sizes:
                exp_obs_count = \
                    self.PointEstimator.estimateExpectedObservationCount(size,
                            n, abundance_freqs, num_obs,
                            self.FullRichnessEstimator)
                exp_obs_count_se = \
                    self.PointEstimator.estimateExpectedObservationCountStdErr(
                            size, n, abundance_freqs, num_obs, exp_obs_count,
                            self.FullRichnessEstimator)

                size_results.append((size, exp_obs_count, exp_obs_count_se))
            per_sample_results.append(size_results)

        return per_sample_results


class AbstractFullRichnessEstimator(object):
    def estimateFullRichness(self, abundance_frequency_counts,
                             observation_count):
        # S_est = S_obs + f_hat_0
        return observation_count + self.estimateUnobservedObservationCount(
                abundance_frequency_counts)

    def estimateUnobservedObservationCount(self, abundance_frequency_counts):
        raise NotImplementedError("Subclasses must implement "
                                  "estimateUnobservedObservationCount.")


class Chao1FullRichnessEstimator(AbstractFullRichnessEstimator):
    def estimateUnobservedObservationCount(self, abundance_frequency_counts):
        # Based on equation 15a and 15b of Colwell 2012.
        f1 = abundance_frequency_counts[0]
        f2 = abundance_frequency_counts[1]

        if f1 < 0 or f2 < 0:
            raise ValueError("Encountered a negative f1 or f2 value, which is "
                             "invalid.")

        if f2 > 0:
            estimated_unobserved_count = f1**2 / (2 * f2)
        else:
            estimated_unobserved_count = (f1 * (f1 - 1)) / (2 * (f2 + 1))

        return estimated_unobserved_count


class AbstractPointEstimator(object):
    def estimateExpectedObservationCount(self, m, n, fk, s_obs,
                                         full_richness_estimator):
        raise NotImplementedError("Subclasses must implement "
                                  "estimateExpectedObservationCount.")

    def estimateExpectedObservationCountStdErr(self, m, n, fk, s_obs, s_m,
                                               full_richness_estimator):
        raise NotImplementedError("Subclasses must implement "
                                  "estimateExpectedObservationCountStdErr.")


class MultinomialPointEstimator(AbstractPointEstimator):
    def estimateExpectedObservationCount(self, m, n, fk, s_obs,
                                         full_richness_estimator):
        # Equation 4 in Colwell 2012
        accumulation = 0

        for k in range(1, n + 1):
            alpha_km = self._calculate_alpha_km(n, k, m)

            # k is 1..n while fk idxs are 0..n-1.
            accumulation += alpha_km * fk[k - 1]

        return s_obs - accumulation

    def estimateExpectedObservationCountStdErr(self, m, n, fk, s_obs, s_m,
                                               full_richness_estimator):
        # Equation 5 in Colwell 2012 gives unconditional variance, but they
        # report the standard error (SE) (which is the same as the standard
        # deviation in this case) in their tables and use this to construct
        # confidence intervals. Thus, we compute SE as sqrt(variance).
        s_est = full_richness_estimator.estimateFullRichness(fk, s_obs)
        accumulation = 0

        for k in range(1, n + 1):
            alpha_km = self._calculate_alpha_km(n, k, m)

            # k is 1..n while fk idxs are 0..n-1.
            accumulation += (((1 - alpha_km)**2) * fk[k - 1])

        # Convert variance to standard error.
        return sqrt(accumulation - (s_m**2 / s_est))

    def _calculate_alpha_km(self, n, k, m):
        alpha_km = 0

        if k <= (n - m):
            alpha_km = ((factorial(n - k) * factorial(n - m)) /
                        (factorial(n) * factorial(n - k - m)))

        return alpha_km
