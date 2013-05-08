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

from biom.util import compute_counts_per_sample_stats

from cogent.maths.stats.distribution import ndtri

from numpy import ceil, empty, sqrt

class EmptyTableError(Exception):
    pass


class EmptySampleError(Exception):
    pass


class ObservationRichnessEstimator(object):
    def __init__(self, biom_table, point_estimator_cls):
        if biom_table.isEmpty():
            raise EmptyTableError("The input BIOM table cannot be empty.")

        self._biom_table = biom_table
        self._point_estimator_cls = point_estimator_cls

    def getSampleCount(self):
        return len(self._biom_table.SampleIds)

    def __call__(self, start=1, stop=None, num_steps=10,
                 confidence_level=0.95):
        results = RichnessEstimatesResults()

        for samp_data, samp_id, _ in self._biom_table.iterSamples():
            point_estimator = self._point_estimator_cls(samp_data)
            ref_indiv_count = point_estimator.getTotalIndividualCount()

            # stop is inclusive. If the original individual count isn't
            # included in this range, add it in the correct spot.
            sizes = self._get_points_to_estimate(ref_indiv_count, start, stop,
                                                 num_steps)
            results.addSample(samp_id, ref_indiv_count)

            for size in sizes:
                exp_obs_count, std_err, ci_low, ci_high = point_estimator(size,
                        confidence_level=confidence_level)
                results.addSampleEstimate(samp_id, size, exp_obs_count,
                                          std_err, ci_low, ci_high)
        return results

    def _get_points_to_estimate(self, reference_individual_count, start=1,
                                stop=None, num_steps=10):
        if stop is None:
            # Compute base sample size as stopping point.
            min_size, max_size, _, _, _ = compute_counts_per_sample_stats(
                    self._biom_table)
            stop = int(max(2 * min_size, max_size))

        if start < 1 or num_steps < 1:
            raise ValueError("The minimum individual count and number of "
                             "steps must both be greater than or equal to 1.")

        if start > stop:
            raise ValueError("The minimum individual count must be less than "
                             "or equal to the maximum individual count.")

        step_size = max((stop - start) // num_steps, 1)

        points = range(start, stop + 1, step_size)
        if reference_individual_count not in points:
            insort(points, reference_individual_count)

        return points


class AbstractPointEstimator(object):
    def __init__(self, sample_data):
        n = self._calculate_total_individual_count(sample_data)
        if n < 1:
            raise EmptySampleError("Encountered a sample without any recorded "
                                   "observations.")
        else:
            self._n = n

        self._s_obs = self._calculate_observation_count(sample_data)
        self._fk = self._calculate_abundance_frequency_counts(sample_data, n)

    def getTotalIndividualCount(self):
        return self._n

    def getObservationCount(self):
        return self._s_obs

    def getAbundanceFrequencyCounts(self):
        return self._fk

    def __call__(self, size):
        raise NotImplementedError("Subclasses must implement __call__.")

    def _calculate_total_individual_count(self, sample_data):
        return int(sample_data.sum(0))

    def _calculate_observation_count(self, sample_data):
        return int((sample_data > 0).sum(0))

    def _calculate_abundance_frequency_counts(self, sample_data, n):
        fk = defaultdict(int)

        for i in range(1, n + 1):
            fk_i = (sample_data == i).sum(0)

            if fk_i > 0:
                fk[i] = int(fk_i)

        return fk


class Chao1MultinomialPointEstimator(AbstractPointEstimator):
    def __init__(self, sample_data):
        super(Chao1MultinomialPointEstimator, self).__init__(sample_data)
        self._f_hat = self._calculate_f_hat(self.getAbundanceFrequencyCounts())
        #self._cov_matrix = self._calculate_covariance_matrix(
        #        self.getAbundanceFrequencyCounts(),
        #        self.getTotalIndividualCount(),
        #        self.estimateFullRichness())

    def estimateUnobservedObservationCount(self):
        return self._f_hat

    def estimateFullRichness(self):
        # S_est = S_obs + f_hat
        return self.getObservationCount() + \
                self.estimateUnobservedObservationCount()

    def __call__(self, size, confidence_level=0.95):
        if confidence_level <= 0 or confidence_level >= 1:
            raise ValueError("Invalid confidence level: %.4f. Must be between "
                             "zero and one (exclusive)." % confidence_level)

        m = size
        fk = self.getAbundanceFrequencyCounts()
        n = self.getTotalIndividualCount()
        s_obs = self.getObservationCount()
        s_est = self.estimateFullRichness()

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
            s_est = self.estimateFullRichness()
            accumulation = 0

            for k in range(1, n + 1):
                alpha_km = self._calculate_alpha_km(n, k, m)
                accumulation += (((1 - alpha_km)**2) * fk[k])

            # Convert variance to standard error.
            std_err = sqrt(accumulation - (estimate**2 / s_est))
        else:
            # Equation 9 in Colwell 2012.
            m_star = m - n
            f1 = fk[1]
            f2 = fk[2]
            f_hat = self.estimateUnobservedObservationCount()
            pd_f1 = self._partial_derivative_f1(f1, f2, m_star, n)
            pd_f2 = self._partial_derivative_f2(f1, f2, m_star, n)

            estimate = s_obs + f_hat * (1 - (1 - (f1 / (n * f_hat))) ** m_star)

            # Equation 10 in Colwell 2012.
            accumulator = 0
            for i in range(1, n + 1):
                if i > 2:
                    pd_i = 1
                elif i == 1:
                    pd_i = pd_f1
                elif i == 2:
                    pd_i = pd_f2

                for j in range(1, n + 1):
                    if j > 2:
                        pd_j = 1
                    elif j == 1:
                        pd_j = pd_f1
                    elif j == 2:
                        pd_j = pd_f2

                    cov = self._calculate_covariance(fk[i], fk[j], s_est, i==j)
                    accumulator += (pd_i * pd_j * cov)

            std_err = sqrt(accumulator)

        
        z_crit = abs(ndtri((1 - confidence_level) / 2))
        ci_bound = z_crit * std_err

        return estimate, std_err, estimate - ci_bound, estimate + ci_bound

    def _calculate_f_hat(self, fk):
        # Based on equation 15a and 15b of Colwell 2012.
        f1 = fk[1]
        f2 = fk[2]

        if f1 < 0 or f2 < 0:
            raise ValueError("Encountered a negative f1 or f2 value, which is "
                             "invalid.")

        if f1 > 0 and f2 > 0:
            f_hat = f1 ** 2 / (2 * f2)
        else:
            f_hat = (f1 * (f1 - 1)) / (2 * (f2 + 1))

        return f_hat

    def _partial_derivative_f1(self, f1, f2, m_star, n):
        if f1 > 0 and f2 > 0:
            a_0 = self._calculate_a_0(f1, f2, n)
            term1 = (m_star * a_0 ** (m_star - 1)) / n
            term2 = (f1 * (1 - a_0 ** m_star)) / f2
            return 1 - term1 + term2
        else:
            a_1 = self._calculate_a_1(f1, f2, n)
            term1 = (m_star * f1) * a_1 ** (m_star - 1)
            term2 = n * (f1 - 1)
            term3 = (f1 - 1) * (1 - a_1 ** m_star)
            term4 = 2 * (f2 + 1)
            term5 = f1 * (1 - a_1 ** m_star)
            return 1 - (term1 / term2) + (term3 / term4) + (term5 / term4)

    def _partial_derivative_f2(self, f1, f2, m_star, n):
        if f1 > 0 and f2 > 0:
            a_0 = self._calculate_a_0(f1, f2, n)
            term1 = (f1 ** 2) * (1 - a_0 ** m_star)
            term2 = 2 * (f2 ** 2)
            term3 = (m_star * f1) * (a_0 ** (m_star - 1))
            term4 = n * f2
            return 1 - (term1 / term2) + (term3 / term4)
        else:
            a_1 = self._calculate_a_1(f1, f2, n)
            term1 = (m_star * f1) * a_1 ** (m_star - 1)
            term2 = n * (f2 + 1)
            term3 = (f1 * (f1 - 1)) * (1 - a_1 ** m_star)
            term4 = 2 * (f2 + 1) ** 2
            return 1 + (term1 / term2) - (term3 / term4)

    def _calculate_a_0(self, f1, f2, n):
        return 1 - ((2 * f2) / (n * f1))

    def _calculate_a_1(self, f1, f2, n):
        return 1 - ((2 * (f2 + 1)) / (n * (f1 - 1)))

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

    def _calculate_covariance(self, f_i, f_j, s_est, same_var):
        if same_var:
            cov = f_i * (1 - f_i / s_est)
        else:
            cov = -(f_i * f_j) / s_est

        return cov

    def _calculate_covariance_matrix(self, fk, n, s_est):
        result = empty((n, n))

        for i in range(0, n):
            f_i = fk[i + 1]

            for j in range(0, i + 1):
                if i != j:
                    f_j = fk[j + 1]
                    cov = -(f_i * f_j) / s_est
                    result[i, j] = cov
                    result[j, i] = cov
                else:
                    result[i, j] = f_i * (1 - f_i / s_est)

        return result


class RichnessEstimatesResults(object):
    _default_header = ['SampleID', 'Size', 'Estimate', 'Std Err', 'CI (lower)',
                       'CI (upper)']
    _num_cols = len(_default_header)

    def __init__(self):
        # sample ID -> (ref individual count,
        #               {size -> (estimate, std err, ci_low, ci_high)})
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
                estimate, std_err, ci_low, ci_high = \
                        self._data[sample_id][1][size]
                results.append((size, estimate, std_err, ci_low, ci_high))
            return results
        else:
            raise ValueError("Unknown sample '%s'." % sample_id)

    def addSample(self, sample_id, reference_individual_count):
        if sample_id in self._data:
            raise ValueError("Sample '%s' has already been added." % sample_id)
        else:
            self._data[sample_id] = (reference_individual_count, {})

    def addSampleEstimate(self, sample_id, size, estimate, std_err, ci_low,
                          ci_high):
        if sample_id in self._data:
            estimates = self._data[sample_id][1]

            if size in estimates:
                raise ValueError("An estimate for sample '%s' already exists "
                                 "at size %d." % (sample_id, size))
            else:
                estimates[size] = (estimate, std_err, ci_low, ci_high)
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
