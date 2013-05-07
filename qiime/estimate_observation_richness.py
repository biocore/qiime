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
from collections import defaultdict
from csv import writer
from itertools import izip

from numpy import ceil, sqrt

class EmptyTableError(Exception):
    pass

class EmptySampleError(Exception):
    pass

class ObservationRichnessEstimator(object):
    def __init__(self, PointEstimator, biom_table):
        if biom_table.isEmpty():
            raise EmptyTableError("The input BIOM table cannot be empty.")

        self._biom_table = biom_table

        for c in self.getTotalIndividualCounts():
            if c < 1:
                raise EmptySampleError("Encountered a sample without any "
                                       "recorded observations.")

        self.PointEstimator = PointEstimator

    def getSampleCount(self):
        return len(self._biom_table.SampleIds)

    def getTotalIndividualCounts(self):
        return [int(e) for e in self._biom_table.sum(axis='sample')]

    def getObservationCounts(self):
        return [(e > 0).sum(0) for e in self._biom_table.iterSampleData()]

    def getAbundanceFrequencyCounts(self):
        for samp_data, num_individuals in izip(
                self._biom_table.iterSampleData(),
                self.getTotalIndividualCounts()):
            samp_abundance_freq_count = defaultdict(int)

            for i in range(1, num_individuals + 1):
                fk_i = (samp_data == i).sum(0)
                if fk_i > 0:
                    samp_abundance_freq_count[i] = fk_i
            yield samp_abundance_freq_count

    def __call__(self, start=1, stop=None, step_size=None):
        results = RichnessEstimatesResults()

        for samp_id, num_obs, n, abundance_freqs in izip(
                self._biom_table.SampleIds,
                self.getObservationCounts(),
                self.getTotalIndividualCounts(),
                self.getAbundanceFrequencyCounts()):
            # stop is inclusive. If the original individual count isn't
            # included in this range, add it in the correct spot.
            sizes = self._get_points_to_estimate(start, stop, step_size, n)
            results.addSample(samp_id, n)

            for size in sizes:
                exp_obs_count, std_err = self.PointEstimator(size, n,
                                                             abundance_freqs,
                                                             num_obs)
                results.addSampleEstimate(samp_id, size, exp_obs_count,
                                          std_err)
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


class AbstractPointEstimator(object):
    def __call__(self, m, n, fk, s_obs):
        raise NotImplementedError("Subclasses must implement __call__.")


class Chao1MultinomialPointEstimator(AbstractPointEstimator):
    def __call__(self, m, n, fk, s_obs):
        if m <= n:
            # Equation 4 in Colwell 2012.
            accumulation = 0

            for k in range(1, n + 1):
                alpha_km = self._calculate_alpha_km(n, k, m)
                accumulation += alpha_km * fk[k]

            estimate = s_obs - accumulation

            # Equation 5 in Colwell 2012 gives unconditional variance, but they
            # report the standard error (SE) (which is the same as the standard
            # deviation in this case) in their tables and use this to construct
            # confidence intervals. Thus, we compute SE as sqrt(variance).
            s_est = self.estimateFullRichness(fk, s_obs)
            accumulation = 0

            for k in range(1, n + 1):
                alpha_km = self._calculate_alpha_km(n, k, m)
                accumulation += (((1 - alpha_km)**2) * fk[k])

            # Convert variance to standard error.
            std_err = sqrt(accumulation - (estimate**2 / s_est))
        else:
            # Equation 9 in Colwell 2012.
            m_star = m - n
            f_hat = self.estimateUnobservedObservationCount(fk)
            f1 = fk[1]

            estimate = s_obs + f_hat * (1 - (1 - (f1 / (n * f_hat)))**m_star)

            # Equation 10 in Colwell 2012.
            a0 = self._calculate_a0(f1, f_hat, n)
            a = self._calculate_a(m_star, a0, n)
            b = self._calculate_b(m_star, a0, n)

            term1 = 0
            for i in range(1, n + 1):
                for j in range(1, n + 1):
                    term1 += self._calculate_covariance(fk[i], fk[j], s_obs,
                                                        f_hat, i==j)

            term2 = 0
            j = 1
            for i in range(1, n + 1):
                term2 += self._calculate_covariance(fk[i], fk[j], s_obs, f_hat,
                                                    i==j)
            term2 *= 2 * a

            term3 = 0
            j = 2
            for i in range(1, n + 1):
                term3 += self._calculate_covariance(fk[i], fk[j], s_obs, f_hat,
                                                    i==j)
            term3 *= 2 * b

            term4 = (a ** 2) * self._calculate_covariance(fk[1], fk[1], s_obs,
                                                          f_hat, True)

            term5 = (2 * a * b) * self._calculate_covariance(fk[2], fk[1],
                                                             s_obs, f_hat,
                                                             False)

            term6 = (b ** 2) * self._calculate_covariance(fk[2], fk[2], s_obs,
                                                          f_hat, True)

            variance_est = term1 + term2 - term3 + term4 - term5 + term6
            std_err = sqrt(variance_est)

        return estimate, std_err

    def estimateFullRichness(self, abundance_frequency_counts,
                             observation_count):
        # S_est = S_obs + f_hat_0
        return observation_count + self.estimateUnobservedObservationCount(
                abundance_frequency_counts)

    def estimateUnobservedObservationCount(self, abundance_frequency_counts):
        # Based on equation 15a and 15b of Colwell 2012.
        f1 = abundance_frequency_counts[1]
        f2 = abundance_frequency_counts[2]

        if f1 < 0 or f2 < 0:
            raise ValueError("Encountered a negative f1 or f2 value, which is "
                             "invalid.")

        if f2 > 0:
            estimated_unobserved_count = f1**2 / (2 * f2)
        else:
            estimated_unobserved_count = (f1 * (f1 - 1)) / (2 * (f2 + 1))

        return estimated_unobserved_count

    def _calculate_alpha_km(self, n, k, m):
        alpha_km = 0
        diff = n - m

        if k <= diff:
            alpha_km = ((self._factorial(n - k) * self._factorial(diff)) /
                        (self._factorial(n) * self._factorial(diff - k)))

        return alpha_km

    def _factorial(self, n, cache=[1, 1]):
        # Taken and modified from
        # http://code.activestate.com/recipes/577241-faster-factorial
        if n < len(cache):
            return cache[n]
        else:
            last = len(cache) - 1
            total = cache[last]

            for i in range(last + 1, n + 1):
                total *= i
                cache.append(total)

            return total

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
    _default_header = ['SampleID', 'Size', 'Estimate', 'Std Err']
    _num_cols = len(_default_header)

    def __init__(self):
        # sample ID -> (ref individual count, {size -> (estimate, std err)})
        self._data = {}

    def getSampleCount(self):
        return len(self._data)

    def getReferenceIndividualCount(self, sample_id):
        if sample_id in self._data:
            return self._data[sample_id][0]
        else:
            raise ValueError("Unknown sample '%s'." % sample_id)

    def getEstimates(self, sample_id):
        if sample_id in self._data:
            results = []
            for size in sorted(self._data[sample_id][1]):
                estimate, std_err = self._data[sample_id][1][size]
                results.append((size, estimate, std_err))
            return results
        else:
            raise ValueError("Unknown sample '%s'." % sample_id)

    def addSample(self, sample_id, reference_individual_count):
        if sample_id in self._data:
            raise ValueError("Sample '%s' has already been added." % sample_id)
        else:
            self._data[sample_id] = (reference_individual_count, {})

    def addSampleEstimate(self, sample_id, size, estimate, std_err):
        if sample_id in self._data:
            estimates = self._data[sample_id][1]

            if size in estimates:
                raise ValueError("An estimate for sample '%s' already exists "
                                 "at size %d." % (sample_id, size))
            else:
                estimates[size] = (estimate, std_err)
        else:
            raise ValueError("An estimate of size %d was provided for an "
                             "unknown sample '%s'." % (size, sample_id))

    def toTable(self, out_f, delimiter='\t', header=None):
        if header is None:
            header = self._default_header
        else:
            if len(header) != self._num_cols:
                raise ValueError("The supplied header must have exactly %d "
                                 "values." % self._num_cols)

        table_writer = writer(out_f, delimiter=delimiter, lineterminator='\n')
        table_writer.writerow(header)

        for sample_id in sorted(self._data):
            estimates = self.getEstimates(sample_id)
            table_writer.writerows([[sample_id] + list(row)
                                     for row in estimates])
