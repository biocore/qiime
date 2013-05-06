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

from bisect import insort
from csv import writer
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

    def __call__(self, start=1, stop=None, step_size=None):
        #results = RichnessEstimatesResults()
        results = []
        orig_indiv_counts = {}

        for samp_id, samp_data, num_obs, n, abundance_freqs in zip(
                self._biom_table.SampleIds,
                self._biom_table.iterSampleData(),
                self.getObservationCounts(),
                self.getTotalIndividualCounts(),
                self.getAbundanceFrequencyCounts()):
            # TODO samp_data not necessary?
            samp_data = samp_data[samp_data > 0]
            orig_indiv_counts[samp_id] = n

            # stop is inclusive. If the original individual count isn't
            # included in this range, add it in the correct spot.
            sizes = self._get_points_to_estimate(start, stop, step_size, n)

            for size in sizes:
                exp_obs_count = \
                    self.PointEstimator.estimateExpectedObservationCount(size,
                            n, abundance_freqs, num_obs,
                            self.FullRichnessEstimator)
                exp_obs_count_se = \
                    self.PointEstimator.estimateExpectedObservationCountStdErr(
                            size, n, abundance_freqs, num_obs, exp_obs_count,
                            self.FullRichnessEstimator)

                results.append((samp_id, size, exp_obs_count,
                                exp_obs_count_se))

        #return RichnessEstimatesResults(results, orig_indiv_counts)
        return results

    def _get_points_to_estimate(self, start, stop, step_size,
                                reference_individual_count):
        if start < 1 or step_size < 1:
            raise ValueError("The minimum observation count and step size "
                             "must both be greater than or equal to 1.")

        if start > stop:
            raise ValueError("The minimum observation count must be less than "
                             "or equal to the maximum observation count.")

        points = range(start, stop + 1, step_size)
        if reference_individual_count not in points:
            insort(points, reference_individual_count)
        return points


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

    #def evaluatePartialDerivative(self, abundance_frequency_counts):


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
        if m <= n:
            # Equation 4 in Colwell 2012.
            accumulation = 0

            for k in range(1, n + 1):
                alpha_km = self._calculate_alpha_km(n, k, m)

                # k is 1..n while fk idxs are 0..n-1.
                accumulation += alpha_km * fk[k - 1]

            estimate = s_obs - accumulation
        else:
            # Equation 9 in Colwell 2012.
            m_star = m - n
            f_hat = \
                full_richness_estimator.estimateUnobservedObservationCount(fk)
            f1 = fk[0]

            estimate = s_obs + f_hat * (1 - (1 - (f1 / (n * f_hat)))**m_star)

        return estimate

    def estimateExpectedObservationCountStdErr(self, m, n, fk, s_obs, s_m,
                                               full_richness_estimator):
        if m <= n:
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
            std_err_est = sqrt(accumulation - (s_m**2 / s_est))
        else:
            # Equation 10 in Colwell 2012.
            m_star = m - n
            f_hat = \
                full_richness_estimator.estimateUnobservedObservationCount(fk)
            a0 = self._calculate_a0(fk[0], f_hat, n)
            a = self._calculate_a(m_star, a0, n)
            b = self._calculate_b(m_star, a0, n)

            term1 = 0
            for i in range(1, n + 1):
                for j in range(1, n + 1):
                    term1 += self._calculate_covariance(fk[i - 1], fk[j - 1],
                                                        s_obs, f_hat, i==j)

            term2 = 0
            j = 1
            for i in range(1, n + 1):
                term2 += self._calculate_covariance(fk[i - 1], fk[j - 1],
                                                    s_obs, f_hat, i==j)
            term2 *= 2 * a

            term3 = 0
            j = 2
            for i in range(1, n + 1):
                term3 += self._calculate_covariance(fk[i - 1], fk[j - 1],
                                                    s_obs, f_hat, i==j)
            term3 *= 2 * b

            term4 = (a ** 2) * self._calculate_covariance(fk[0], fk[0], s_obs,
                                                          f_hat, True)

            term5 = (2 * a * b) * self._calculate_covariance(fk[1], fk[0],
                                                             s_obs, f_hat,
                                                             False)

            term6 = (b ** 2) * self._calculate_covariance(fk[1], fk[1], s_obs,
                                                          f_hat, True)

            variance_est = term1 + term2 - term3 + term4 - term5 + term6
            std_err_est = sqrt(variance_est)

        return std_err_est

    def _calculate_alpha_km(self, n, k, m):
        alpha_km = 0

        if k <= (n - m):
            alpha_km = ((factorial(n - k) * factorial(n - m)) /
                        (factorial(n) * factorial(n - k - m)))

        return alpha_km

    def _calculate_covariance(self, f_i, f_j, s_obs, f_hat, same_var):
        if same_var:
            cov = f_i * (1 - f_i / (s_obs + f_hat))
        else:
            cov = -(f_i * f_j) / (s_obs + f_hat)

        return cov

    def _calculate_a0(self, f_1, f_hat, n):
        return f_1 / (n * f_hat + f_1)

    def _calculate_a(self, m_star, a0, n):
        term1 = (2 * (1 - a0)) / n
        term2 = 1 / a0
        term3 = (1 - a0) ** m_star
        term4 = term2 + (m_star / 2)
        return term1 * (term2 - (term3 * term4))

    def _calculate_b(self, m_star, a0, n):
        term1 = (2 * ((1 - a0) ** 2)) / ((n ** 2) * a0)
        term2 = 1 / a0
        term3 = (1 - a0) ** m_star
        term4 = term2 + m_star
        return term1 * (term2 - (term3 * term4))


class RichnessEstimatesResults(object):
    _num_cols = 4
    _default_header = ['SampleID', 'Size', 'Estimate', 'Std Err']

    def __init__(self):
        self.RawEstimatesData = []
        self._samples = {}

    def getSampleCount(self):
        return len(self._reference_sample_counts)

    def addSample(self, sample_id, individual_count):
        if sample_id in self._samples:
            raise ValueError("Sample '%s' has already been added to the list "
                             "of results." % sample_id)
        else:
            self._samples[sample_id] = individual_count

    def addSampleEstimate(self, sample_id, size, estimate, estimate_std_err):
        if sample_id in self._samples:
            self.RawEstimatesData.append([sample_id, size, estimate,
                                          estimate_std_err])

    def saveRawData(self, out_fp, delimiter='\t', header=None):
        if header is None:
            header = self._default_header
        else:
            if len(header) != self._num_cols:
                raise ValueError("The supplied header must have exactly %d "
                                 "values." % self._num_cols)

        with open(out_fp, 'wb') as out_f:
            raw_data_writer = writer(out_f, delimiter=delimiter)

            raw_data_writer.writerow(header)
            raw_data_writer.writerows(self._estimates_data)
