#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2013, The QIIME Project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Contains functionality to estimate the observation richness of samples."""

from bisect import insort
from collections import defaultdict
from csv import writer

from biom.util import compute_counts_per_sample_stats
from cogent.maths.stats.distribution import ndtri
from numpy import empty, ones, sqrt, tensordot


class EmptyTableError(Exception):
    pass


class EmptySampleError(Exception):
    pass


class ObservationRichnessEstimator(object):

    """Class to estimate richness of samples in a table at varying depths.

    This is the main class that users of this module will interact with. It
    handles computing estimates at points within a user-specified range, and
    collates these point estimates into an object that provides various output
    formatting functionality.

    This class is agnostic to what type of point estimator is used. The type of
    point estimator is specified upon instantiation.
    """

    def __init__(self, biom_table, point_estimator_cls):
        """Construct a richness estimator instance.

        Arguments:
            biom_table - input BIOM table containing samples to estimate the
                richness of
            point_estimator_cls - the class of point estimator to use. Should
                be a subclass of AbstractPointEstimator, e.g.
                Chao1MultinomialPointEstimator
        """
        if biom_table.is_empty():
            raise EmptyTableError("The input BIOM table cannot be empty.")

        self._biom_table = biom_table
        self._point_estimator_cls = point_estimator_cls

    def getSampleCount(self):
        """Return the number of samples in the table."""
        return len(self._biom_table.ids())

    def __call__(self, start=1, stop=None, num_steps=10,
                 confidence_level=0.95):
        """Compute estimates using the provided range and collate results.

        For each depth, the estimate, its standard error, and confidence
        interval will be computed.

        If the reference sample size isn't included in this range, it will be
        added.

        Returns a RichnessEstimatesResults instance, which contains the
        collated results.

        Arguments:
            start - starting depth of the range
            stop - stopping depth of the range. May not necessarily be
                included; it is an upper bound. If None, the base sample size
                will be used (see script documentation for more details)
            num_steps - the number of steps to make between start and stop
            confidence_level - a float between 0 and 1 (exclusive) indicating
                the confidence level to use in the confidence intervals
        """
        results = RichnessEstimatesResults()

        for samp_data, samp_id, _ in self._biom_table.iter(axis='sample'):
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
        """Returns depths/sizes to estimate."""
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

    """Abstract class for a point estimator.

    A point estimator is concerned with a *single* sample's observation data
    (e.g. counts vector). It provides an estimate of richness given a
    depth/size.

    Subclasses must implement __call__.
    """

    def __init__(self, sample_data):
        """Construct an abstract point estimator.

        Arguments:
            sample_data - a 1-D numpy array containing observation counts in a
                sample
        """
        n = self._calculate_total_individual_count(sample_data)
        if n < 1:
            raise EmptySampleError("Encountered a sample without any recorded "
                                   "observations.")
        else:
            self._n = n

        self._s_obs = self._calculate_observation_count(sample_data)
        self._fk = self._calculate_abundance_frequency_counts(sample_data, n)

    def getTotalIndividualCount(self):
        """Return the reference sample size (total number of individuals).

        This is n in Colwell 2012.
        """
        return self._n

    def getObservationCount(self):
        """Return the number of (observed) observations in the sample.

        This is S_obs in Colwell 2012.
        """
        return self._s_obs

    def getAbundanceFrequencyCounts(self):
        """Return the abundance frequency counts in the sample.

        This is a defaultdict instance mapping abundance to frequency count.
        For example, if the sample contains two observations, each with an
        individual count of 3:

            {3: 2}

        This is f_k in Colwell 2012.
        """
        return self._fk

    def __call__(self, size, confidence_level=0.95):
        """Estimate the richness at the given size.

        Must return a tuple:

            (estimate, std_err, ci_low, ci_high)

        """
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

    """Point estimator using the multinomial model and Chao1.

    Supports both interpolation/rarefaction and extrapolation by implementing
    equations 4, 5, 9, 10, 15a, and 15b in Colwell et al. (2012), Journal of
    Plant Ecology.
    """

    def __init__(self, sample_data):
        super(Chao1MultinomialPointEstimator, self).__init__(sample_data)
        self._f_hat = self._calculate_f_hat(self.getAbundanceFrequencyCounts())

        n = self.getTotalIndividualCount()
        self._cov_matrix = self._calculate_covariance_matrix(
            self.getAbundanceFrequencyCounts(), n,
            self.estimateFullRichness())

        self._pd_matrix = ones((n, n))

    def estimateUnobservedObservationCount(self):
        """Return estimated number of observations not found in this sample.

        This is f_0_hat_Chao1 in Colwell 2012 (equations 15a and 15b).
        """
        return self._f_hat

    def estimateFullRichness(self):
        """Return estimated total number of observations in this sample.

        This is S_est_Chao1 in Colwell 2012, which is S_obs + f_0_hat_Chao1.
        """
        # S_est = S_obs + f_hat
        return self.getObservationCount() + \
            self.estimateUnobservedObservationCount()

    def __call__(self, size, confidence_level=0.95):
        if confidence_level <= 0 or confidence_level >= 1:
            raise ValueError("Invalid confidence level: %.4f. Must be between "
                             "zero and one (exclusive)." % confidence_level)

        # We'll use the variable names from Colwell 2012 for clarity and
        # brevity.
        m = size
        fk = self.getAbundanceFrequencyCounts()
        n = self.getTotalIndividualCount()
        s_obs = self.getObservationCount()
        s_est = self.estimateFullRichness()

        if m <= n:
            # Interpolation.

            # Equation 4 in Colwell 2012 for the estimate.
            estimate_acc = 0

            # Equation 5 in Colwell 2012 gives unconditional variance, but they
            # report the standard error (SE) (which is the same as the standard
            # deviation in this case) in their tables and use this to construct
            # confidence intervals. Thus, we compute SE as sqrt(variance).
            std_err_acc = 0

            for k in range(1, n + 1):
                alpha_km = self._calculate_alpha_km(n, k, m)
                estimate_acc += alpha_km * fk[k]
                std_err_acc += (((1 - alpha_km) ** 2) * fk[k])

            estimate = s_obs - estimate_acc

            # Convert variance to standard error.
            std_err = sqrt(std_err_acc - (estimate ** 2 / s_est))
        else:
            # Extrapolation.
            m_star = m - n
            f1 = fk[1]
            f2 = fk[2]
            f_hat = self.estimateUnobservedObservationCount()

            try:
                # Equation 9 in Colwell 2012.
                estimate = s_obs + f_hat * (1 -
                                            (1 - (f1 / (n * f_hat))) ** m_star)
            except ZeroDivisionError:
                # This can happen if we have exactly one singleton and no
                # doubletons, or no singletons and no doubletons.
                estimate = None
                std_err = None
            else:
                # Equation 10 in Colwell 2012. I used Wolfram Alpha to
                # calculate the analytic partial derivatives since they weren't
                # provided in the original paper. We have two partial
                # derivatives, wrt f1 and f2, that we really care about. All
                # other partial derivatives (e.g. wrt f3, f4, etc.) get a value
                # of 1.
                pd_f1 = self._partial_derivative_f1(f1, f2, m_star, n)
                pd_f2 = self._partial_derivative_f2(f1, f2, m_star, n)
                pd_f1f2 = pd_f1 * pd_f2

                # To do this efficiently, here's the algorithm:
                #
                # 1) Create nxn array filled with ones. Each element represents
                #    the multiplication of two partial derivatives.
                # 2) Fill in only what we need: the multiplication of partial
                #    derivatives wrt f1 and f2.
                # 3) Do an element-wise multiply between our partial derivative
                #    matrix and the covariance matrix. tensordot does this and
                #    also sums the result, which is exactly what we need. In
                #    the end, we've summed all n^2 elements, each of which are
                #    (pd_fi * pd_fj * cov_ij).
                self._pd_matrix[0, :] = pd_f1
                self._pd_matrix[1, :] = pd_f2
                self._pd_matrix[:, 0] = pd_f1
                self._pd_matrix[:, 1] = pd_f2

                self._pd_matrix[0, 0] = pd_f1 ** 2
                self._pd_matrix[0, 1] = pd_f1f2
                self._pd_matrix[1, 0] = pd_f1f2
                self._pd_matrix[1, 1] = pd_f2 ** 2

                std_err = sqrt(tensordot(self._pd_matrix, self._cov_matrix))

        # Compute CI based on std_err.
        ci_low = None
        ci_high = None
        if std_err is not None:
            # z_crit will be something like 1.96 for 95% CI.
            z_crit = abs(ndtri((1 - confidence_level) / 2))
            ci_bound = z_crit * std_err
            ci_low = estimate - ci_bound
            ci_high = estimate + ci_bound

        return estimate, std_err, ci_low, ci_high

    def _calculate_f_hat(self, fk):
        # Based on equations 15a and 15b in Colwell 2012.
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

    # I lost my sanity somewhere around this point... :P Sorry for anyone that
    # has to read this!

    def _partial_derivative_f1(self, f1, f2, m_star, n):
        """Derived from equation 9 using Wolfram Alpha, wrt f1."""
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
        """Derived from equation 9 using Wolfram Alpha, wrt f2."""
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
        # I made up the names a_0 and a_1 (they're not in the paper) to break
        # out common terms from the above partial derivatives.
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
        """Dynamic programming factorial!

        Using Python's built-in factorial was very slow. This implementation
        cut down a lot of time.

        Taken and modified from:
            http://code.activestate.com/recipes/577241-faster-factorial
        """
        if n < len(cache):
            return cache[n]
        else:
            last = len(cache) - 1
            total = cache[last]

            for i in range(last + 1, n + 1):
                total *= i
                cache.append(total)

            return total

    def _calculate_covariance_matrix(self, fk, n, s_est):
        # This is pretty expensive... need to find a way to speed this method
        # up.
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

    """Container to hold estimates results and provide output formatting.

    Currently only supports writing to a table, but may be able to create plots
    in the future.
    """

    _default_header = ['SampleID', 'Size', 'Estimate', 'Std Err', 'CI (lower)',
                       'CI (upper)']
    _num_cols = len(_default_header)

    def __init__(self):
        """Initialize an empty results container."""
        # sample ID -> (ref individual count,
        #               {size -> (estimate, std err, ci_low, ci_high)})
        self._data = {}

    def getSampleCount(self):
        """Return the number of samples that are currently stored."""
        return len(self._data)

    def getReferenceIndividualCount(self, sample_id):
        """Return the reference sample size for the given sample."""
        if sample_id in self._data:
            return self._data[sample_id][0]
        else:
            raise ValueError("Unknown sample '%s'." % sample_id)

    def getEstimates(self, sample_id):
        """Return a list of estimates (sorted by depth/size)."""
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
        """Add a new sample to the results container.

        Sample cannot already exist in the container.

        Arguments:
            sample_id - the sample ID to add
            reference_individual_count - the reference sample size (n)
        """
        if sample_id in self._data:
            raise ValueError("Sample '%s' has already been added." % sample_id)
        else:
            self._data[sample_id] = (reference_individual_count, {})

    def addSampleEstimate(self, sample_id, size, estimate, std_err, ci_low,
                          ci_high):
        """Add a richness estimate at the given size.

        Sample must already exist in the container.
        """
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

    def toTable(self, out_f, header=None, delimiter='\t', missing='N/A'):
        """Write results in tabular format to a file.

        Arguments:
            out_f - the output file(-like) object to write results to
            header - a list of strings (one per column) to be written as a
                header. If not provided, the class default will be used
            delimiter - string to delimit cells in the table
            missing - string used in place of any cells that are None
        """
        if header is None:
            header = self._default_header
        else:
            if len(header) != self._num_cols:
                raise ValueError("The supplied header must have exactly %d "
                                 "values." % self._num_cols)

        table_writer = writer(out_f, delimiter=delimiter, lineterminator='\n')
        table_writer.writerow(header)

        missing_f = lambda e: missing if e is None else e

        for sample_id in sorted(self._data):
            estimates = self.getEstimates(sample_id)
            table_writer.writerows([[sample_id] + map(missing_f, row)
                                    for row in estimates])
