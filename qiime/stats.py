#!/usr/bin/env python
from __future__ import division

__author__ = "Michael Dwan"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Michael Dwan", "Logan Knecht",
               "Damien Coy", "Levi McCracken", "Andrew Cochran",
               "Jose Carlos Clemente Litran", "Greg Caporaso",
               "Will Van Treuren"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""
This module provides functionality for the application of various statistical
methods to QIIME-formatted datasets.

The module provides an API that allows users to easily apply any number of
statistical analyses and just as easily retrieve the results. The module also
provides a hierarchy of statistical classes that can be inherited from to
create new statistical method implementations.
"""

from scipy.stats import (spearmanr, kruskal, mannwhitneyu, kendalltau,
                         power_divergence, ttest_1samp, ttest_ind)
from scipy.stats.distributions import (chi2, norm, f as fdist, t as tdist)

from scipy.special import ndtri

from collections import defaultdict
from os.path import join
from types import ListType
from copy import deepcopy
from itertools import combinations

from matplotlib import use
use('Agg', warn=False)

from matplotlib.pyplot import figure
from numpy import (argsort, array, ceil, empty, fill_diagonal, finfo,
                   log2, mean, ones, sqrt, tri, unique, zeros, ndarray, floor,
                   median, nan, min as np_min, max as np_max, absolute,
                   arctanh, asarray, e, hstack, isinf, isnan,
                   log, mean, nan, nonzero, sqrt, std, take, tanh,
                   transpose, seterr as np_seterr, var, arange, corrcoef,
                   trace, ravel, float as np_float, finfo, asarray, isnan,
                   isinf, abs)

from numpy.random import permutation, shuffle, randint
from biom.table import Table
from skbio.stats.distance import DistanceMatrix, mantel
from skbio.util import create_dir

from qiime.format import format_p_value_for_num_iters
from qiime.util import MetadataMap, write_biom_table

np_seterr(divide='warn')
MACHEP = finfo(np_float).eps

# Top-level stats functions.

tail_types = ['low', 'high', 'two-sided']
tail_type_desc = {
    'low': ('one-sided (low)', '<'),
    'high': ('one-sided (high)', '>'),
    'two-sided': ('two-sided', '!=')
}


def all_pairs_t_test(labels, dists, tail_type='two-sided',
                     num_permutations=999):
    """Perform two-sample t-test on all pairs of grouped distances.

    Performs Student's two-sample t-test on all pairs of distributions,
    optionally using Monte Carlo permutations to compute the nonparametric
    p-value in addition to the parametric p-value.

    Returns a formatted string (suitable for writing to a file) containing the
    results of the tests.

    This code is based on Jeremy Widmann's
    qiime.make_distance_histograms.monte_carlo_group_distances code from QIIME 1.8.0.

    Arguments:
        labels - list of labels corresponding to each of the distributions
        dists - list of lists, where each inner list is a distribution of
            numbers (observations)
        tail_type - type of hypothesis test to perform. One of 'two-sided',
            'high', or 'low'
        num_permutations - the number of Monte Carlo permutations to use. If
            zero, the nonparametric p-value will not be calculated and will be
            'N/A' in the returned string.
    """
    result = ''

    if len(labels) != len(dists):
        raise ValueError("The number of distribution labels must match the "
                         "number of distributions.")
    if tail_type not in tail_types:
        raise ValueError("Invalid tail type '%s'. Must be one of %r." %
                         (tail_type, tail_types))
    if num_permutations < 0:
        raise ValueError("Invalid number of permutations: %d. Must be greater "
                         "than or equal to zero." % num_permutations)

    result += '# The tests of significance were performed using a ' + \
              tail_type_desc[tail_type][0] + ' Student\'s two-sample t-test.\n'

    result += ('# Alternative hypothesis: Group 1 mean %s Group 2 mean\n' %
               tail_type_desc[tail_type][1])

    if num_permutations > 0:
        result += '# The nonparametric p-values were calculated using ' + \
                  '%d Monte Carlo permutations.\n' % num_permutations
        result += '# The nonparametric p-values contain the correct ' + \
                  'number of significant digits.\n'

    result += '# Entries marked with "N/A" could not be calculated because ' + \
              'at least one of the groups\n# of distances was empty, ' + \
              'both groups each contained only a single distance, or\n' + \
              '# the test could not be performed (e.g. no variance in ' + \
              'groups with the same mean).\nGroup 1\tGroup 2\t' + \
              't statistic\tParametric p-value\tParametric p-value ' + \
              '(Bonferroni-corrected)\tNonparametric p-value\t' + \
              'Nonparametric p-value (Bonferroni-corrected)\n'

    stats = _perform_pairwise_tests(labels, dists, tail_type, num_permutations)
    for stat in stats:
        stat = ['N/A' if e is nan else e for e in stat]
        result += '%s\t%s\t%s\t%s\t%s\t%s\t%s\n' % (stat[0], stat[1], stat[2],
                                                    stat[3], stat[4],
                                                    format_p_value_for_num_iters(stat[5], num_permutations) if
                                                    stat[5] != 'N/A' else 'N/A',
                                                    format_p_value_for_num_iters(stat[6], num_permutations) if
                                                    stat[6] != 'N/A' else 'N/A')
    return result


def _perform_pairwise_tests(labels, dists, tail_type, num_permutations):
    """Perform t-test for all pairs of distributions.

    Computes corrected p-values in addition to uncorrected.
    """
    result = []

    # Compare each pair of distributions, keeping track of the number of actual
    # tests that were successfully performed so that we can correct for
    # multiple comparisons.
    num_tests = 0
    for g1_idx, (g1_label, g1_dist) in enumerate(zip(labels[:-1], dists[:-1])):
        for g2_label, g2_dist in zip(
                labels[(g1_idx + 1):], dists[(g1_idx + 1):]):
            if ((len(g1_dist) == 1 and len(g2_dist) == 1) or
                    (len(g1_dist) < 1 or len(g2_dist) < 1)):
                # Not enough data to run the test.
                obs_t, param_p_val, nonparam_p_val = nan, nan, nan
            else:
                obs_t, param_p_val, _, nonparam_p_val = mc_t_two_sample(
                    g1_dist, g2_dist, tails=tail_type,
                    permutations=num_permutations)
            result.append([g1_label, g2_label, obs_t, param_p_val, None,
                           nonparam_p_val, None])
            if not isnan(obs_t):
                num_tests += 1

    # Correct the p-values for multiple comparisons, now that we know how many
    # tests succeeded.
    for stat in result:
        corr_param_p_val = stat[3]
        if corr_param_p_val is not None and not isnan(corr_param_p_val):
            corr_param_p_val = min(corr_param_p_val * num_tests, 1)
        stat[4] = corr_param_p_val

        corr_nonparam_p_val = stat[5]
        if corr_nonparam_p_val is not None and not isnan(corr_nonparam_p_val):
            corr_nonparam_p_val = min(corr_nonparam_p_val * num_tests, 1)
        stat[6] = corr_nonparam_p_val

    return result


def quantile(data, quantiles):
    """calculates quantiles of a dataset matching a given list of probabilities

    Input:
    data: 1-D list or numpy array with data to calculate the quantiles
    quantiles: list of probabilities, floating point values between 0 and 1

    Output:
    A list of elements drawn from 'data' that corresponding to the list of
    probabilities. This by default is using R. type 7 method for computation of
    the quantiles.
    """

    assert isinstance(data, list) or isinstance(data, ndarray), "Data must be either" +\
        " a Python list or a NumPy 1-D array"
    assert isinstance(quantiles, list) or isinstance(quantiles, ndarray), "Quantiles" +\
        " must be either a Python list or a NumPy 1-D array"
    assert all(map(lambda x: x >= 0 and x <= 1, quantiles)), "All the elements " +\
        "in the quantiles list must be greater than 0 and lower than one"

    # unless the user wanted, do not modify the data
    data = deepcopy(data)

    if not isinstance(data, ndarray):
        data = array(data)
    data.sort()

    output = []
    # if needed different quantile methods could be used
    for one_quantile in quantiles:
        output.append(_quantile(data, one_quantile))

    return output


def _quantile(data, quantile):
    """gets a single quantile value for a dataset using R. type 7 method

    Input:
    data: sorted 1-d numpy array with float or int elements
    quantile: floating point value between 0 and 1

    Output:
    quantile value of data

    This function is based on cogent.maths.stats.util.NumbersI
    """
    index = quantile * (len(data) - 1)
    bottom_index = int(floor(index))
    top_index = int(ceil(index))

    difference = index - bottom_index
    output = (1 - difference) * \
        data[bottom_index] + difference * data[top_index]

    return output


class DistanceMatrixStats(object):

    """Base class for distance matrix-based statistical methods.

    This class provides an interface to setting and accessing an arbitrary
    number of distance matrices. Users of this class can optionally specify the
    number of allowable distance matrices and their minimum allowable size (the
    default is no restrictions on either of these).

    It is the parent class of CorrelationStats.
    """

    def __init__(self, dms, num_dms=-1, min_dm_size=-1):
        """Default constructor.

        Initializes an instance with the provided list of distance matrices.

        Arguments:
            dms - a list of DistanceMatrix objects
            num_dms - the exact number of allowable distance matrices. If -1
                (the default), there is no restriction on how many distance
                matrices the user can set
            min_dm_size - the minimum size that all distance matrices must have
                that are stored by this instance. If -1, no size restriction
        """
        self._num_dms = num_dms
        self._min_dm_size = min_dm_size
        self.DistanceMatrices = dms

    @property
    def DistanceMatrices(self):
        """Returns the list of distance matrices."""
        return self._dms

    @DistanceMatrices.setter
    def DistanceMatrices(self, dms):
        """Sets the list of distance matrices to the supplied list.

        Arguments:
            dms - the new list of distance matrices being assigned
        """
        if not isinstance(dms, ListType):
            raise TypeError("The item passed in as the new list was not a "
                            "list data type.")
        if self._num_dms >= 0 and len(dms) != self._num_dms:
            raise ValueError("Cannot set %d distance matrices. Must provide "
                             "exactly %d distance matrices." % (len(dms),
                                                                self._num_dms))
        for dm in dms:
            if not isinstance(dm, DistanceMatrix):
                raise TypeError(
                    'Invalid type (%s); expected DistanceMatrix' %
                    dm.__class__.__name__)
            if self._min_dm_size >= 0 and dm.shape[0] < self._min_dm_size:
                raise ValueError("Distance matrix of size %dx%d is smaller "
                                 "than the minimum allowable distance matrix "
                                 "size of %dx%d for this analysis." %
                                 (dm.shape[0], dm.shape[0], self._min_dm_size,
                                  self._min_dm_size))
        self._dms = dms

    def __call__(self, num_perms=999):
        """Runs the statistical method and returns relevant results.

        The return value of this method is a python dictionary with arbitrary
        key/value pairs of results, since each statistical method returns
        different results.

        This method returns an empty result set (it is essentially not
        implemented) and should be implemented by subclasses to perform their
        specific statistical analysis. Subclasses should call the parent
        class' __call__ method first to obtain any results from the parent and
        then add more results to the dict that is obtained from the parent.

        Arguments:
            num_perms - the number of permutations to use in the statistical
                method. If the method is not permutation-based, simply ignore
                this argument
        """
        if num_perms < 0:
            raise ValueError("The number of permutations must be greater than "
                             "or equal to zero.")
        return {}


class CorrelationStats(DistanceMatrixStats):
    """Base class for distance matrix correlation statistical methods.

    It is subclassed by correlation methods such as partial Mantel and Mantel
    correlogram that compare two or more distance matrices.

    A valid instance of CorrelationStats must have at least one distance
    matrix, and all distance matrices must have matching dimensions and sample
    IDs (i.e. matching row/column labels). This check is in place to prevent
    the accidental comparison on two distance matrices that have sample IDs in
    different orders. Essentially, all of the distance matrices must be
    "compatible".

    Users of this class can optionally specify the number of allowable distance
    matrices and their minimum allowable size (the default is no restrictions
    on either of these).
    """

    @property
    def DistanceMatrices(self):
        # Must re-declare so we can override property setter below.
        return super(CorrelationStats, self).DistanceMatrices

    @DistanceMatrices.setter
    def DistanceMatrices(self, dms):
        """Sets the list of distance matrices to the supplied list.

        This method overrides the parent method and enforces more checks to
        ensure that at least one distance matrix is provided and that all of
        the distance matrices are compatible.

        Arguments:
            dms - the new list of distance matrices being assigned
        """
        # Must call superclass property setter this way (super doesn't work).
        DistanceMatrixStats.DistanceMatrices.fset(self, dms)

        if len(dms) < 1:
            raise ValueError("Must provide at least one distance matrix.")

        size = dms[0].shape[0]
        sample_ids = dms[0].ids
        for dm in dms:
            if dm.shape[0] != size:
                raise ValueError("All distance matrices must have the same "
                                 "number of rows and columns.")
            if dm.ids != sample_ids:
                raise ValueError("All distance matrices must have matching "
                                 "sample IDs.")


class MantelCorrelogram(CorrelationStats):

    """Class for the Mantel correlogram statistical method.

    This class provides the functionality to run a Mantel correlogram analysis
    on two distance matrices. In a nutshell, the distances are split into
    distance classes and a Mantel test is run over each distance class. A
    Mantel correlogram is created, which is basically a plot of distance
    classes versus Mantel statistics.

    Uses Sturge's rule to determine the number of distance classes, and
    Pearson's method to compute the correlation at each distance class. The
    corrected p-values are computed using Bonferroni correction.
    """

    def __init__(self, eco_dm, geo_dm, alpha=0.05,
                 variable_size_distance_classes=False):
        """Constructs a new MantelCorrelogram instance.

        Arguments:
            eco_dm - a DistanceMatrix object representing the ecological
                distances between samples (e.g. UniFrac distance matrix)
            geo_dm - a DistanceMatrix object representing some other distance
                measure between samples (most commonly geographical distances,
                but could also be distances in pH, temperature, etc.)
            alpha - the alpha value to use when marking the Mantel correlogram
                plot for significance
            variable_size_distance_classes - if True, distance classes (bins)
                will vary in size such that each distance class (bin) will have
                the same number of distances. If False, all distance classes
                will have the same size, though the number of distances in each
                class may not be equal. Having variable-sized distance classes
                can help maintain statistical power if there are large
                differences in the number of distances in each class
        """
        super(MantelCorrelogram, self).__init__([eco_dm, geo_dm], num_dms=2,
                                                min_dm_size=3)
        self.Alpha = alpha
        self.VariableSizeDistanceClasses = variable_size_distance_classes

    @property
    def Alpha(self):
        """Returns the alpha value."""
        return self._alpha

    @Alpha.setter
    def Alpha(self, alpha):
        """Sets the alpha value.

        Arguments:
            alpha - the value of alpha. Must be between 0 and 1, inclusive
        """
        if alpha >= 0 and alpha <= 1:
            self._alpha = alpha
        else:
            raise ValueError("Alpha must be between 0 and 1.")

    def __call__(self, num_perms=999):
        """Runs a Mantel correlogram test over the current distance matrices.

        Returns a dict containing the results. The following keys are set:
            method_name - name of the statistical method
            class_index - list of distance class indices (the center of each
                distance class)
            num_dist - list of the number of distances in each distance class
            mantel_r - list of the Mantel r statistics for each distance class
            mantel_p - list of the p-values for each distance class
            mantel_p_corr - list of the p-values for each distance class,
                corrected for multiple tests
            correlogram_plot - a matplotlib Figure object containing the
                correlogram

        Arguments:
            num_perms - the number of permutations to use when calculating the
                p-values

        Note: This code is heavily based on the implementation of
        mantel.correlog in R's vegan package.
        """
        results = super(MantelCorrelogram, self).__call__(num_perms)
        eco_dm = self.DistanceMatrices[0]
        geo_dm = self.DistanceMatrices[1]
        dm_size = eco_dm.shape[0]

        # Find the number of lower triangular elements (excluding the
        # diagonal).
        num_dists = dm_size * (dm_size - 1) // 2

        # Use Sturge's rule to determine the number of distance classes.
        num_classes = int(ceil(1 + log2(num_dists)))

        # Create the matrix of distance classes. Each element in the matrix
        # contains what distance class the original element is in. Also find
        # the distance class indices, which are the midpoints in each distance
        # class.
        dist_class_matrix, class_indices = self._find_distance_classes(
            geo_dm, num_classes)

        # Start assembling the results.
        results['method_name'] = 'Mantel Correlogram'
        results['class_index'] = []
        results['num_dist'] = []
        results['mantel_r'] = []
        results['mantel_p'] = []

        # Create a model matrix for each distance class, then compute a Mantel
        # test using it and the original eco distance matrix. A model matrix
        # contains ones for each element that is in the current distance class,
        # and zeros otherwise (zeros on the diagonal as well).
        for class_num in range(num_classes):
            results['class_index'].append(class_indices[class_num])
            model_matrix = zeros([dm_size, dm_size], dtype=int)
            for i in range(dm_size):
                for j in range(dm_size):
                    curr_ele = dist_class_matrix[i][j]
                    if curr_ele == class_num and i != j:
                        model_matrix[i][j] = 1
            model_matrix = DistanceMatrix(model_matrix, geo_dm.ids)

            # Count the number of distances in the current distance class.
            num_distances = int(model_matrix.data.sum())
            results['num_dist'].append(num_distances)
            if num_distances == 0:
                results['mantel_r'].append(None)
                results['mantel_p'].append(None)
            else:
                row_sums = model_matrix.data.sum(axis=1)
                row_sums = map(int, row_sums)
                has_zero_sum = 0 in row_sums

                # Only stop running Mantel tests if we've gone through half of
                # the distance classes and at least one row has a sum of zero
                # (i.e. the sample doesn't have any distances that fall in the
                # current class).
                if not (class_num > ((num_classes // 2) - 1) and has_zero_sum):
                    # Compute the correlation coefficient without performing
                    # permutation tests in order to check its sign below.
                    orig_stat, _, _ = mantel(
                        model_matrix, eco_dm, method='pearson',
                        permutations=0, strict=True)

                    # Negate the Mantel r statistic because we are using
                    # distance matrices, not similarity matrices (this is a
                    # necessary step, see Legendre's Numerical Ecology
                    # algorithm reference for more details).
                    results['mantel_r'].append(-orig_stat)

                    # Compute a one-tailed p-value in the direction of the
                    # sign.
                    if orig_stat < 0:
                        tail_type = 'less'
                    else:
                        tail_type = 'greater'

                    _, p_val, _ = mantel(
                        model_matrix, eco_dm, method='pearson',
                        permutations=num_perms, alternative=tail_type,
                        strict=True)

                    results['mantel_p'].append(p_val)
                else:
                    results['mantel_r'].append(None)
                    results['mantel_p'].append(None)

        # Correct p-values for multiple testing.
        results['mantel_p_corr'] = self._correct_p_values(results['mantel_p'])

        # Construct a correlogram of distance class versus mantel correlation
        # statistic and fill in each point that is statistically significant.
        results['correlogram_plot'] = self._generate_correlogram(
            results['class_index'], results['mantel_r'],
            results['mantel_p_corr'])
        return results

    def _find_distance_classes(self, dm, num_classes):
        """Computes a distance class matrix and distance class midpoints.

        Returns a matrix of the same dimensions as the input matrix but each
        element indicates which distance class (0..num_classes-1) the original
        element belongs to. The diagonal will always have a value of -1,
        indicating that it is not apart of any distance class. Also returns a
        list of distance class midpoints.

        Distance classes are determined by the minimum and maximum values in
        the input matrix and the number of specified classes. If
        self.VariableSizeDistanceClasses is True, distance classes will each
        contain the same number of distances (but may vary in size). If False,
        distance classes will be of equal size (but possibly with unequal
        numbers of distances).

        Arguments:
            dm - the input DistanceMatrix object to compute distance classes on
            num_classes - the number of desired distance classes
        """
        if num_classes < 1:
            raise ValueError("Cannot have fewer than one distance class.")

        dm_lower_flat = dm.condensed_form()
        size = dm.shape[0]

        if self.VariableSizeDistanceClasses:
            class_size = int(ceil(len(dm_lower_flat) / num_classes))
            order = argsort(array(dm_lower_flat))

            # Create the matrix of distance classes. Every element in the
            # matrix tells what distance class the original element belongs to.
            # Each element in the original matrix is traversed in sorted
            # (min -> max) order, and the current distance class is incremented
            # once it is "filled" with class_size distances.
            dist_class_matrix = empty([size, size], dtype=int)
            class_indices = []
            curr_class = 0
            class_start = dm_lower_flat[order[0]]
            for i, sorted_idx in enumerate(order):
                row_idx, col_idx = self._find_row_col_indices(sorted_idx)
                class_end = dm_lower_flat[sorted_idx]

                # Matrix is symmetric.
                dist_class_matrix[row_idx][col_idx] = curr_class
                dist_class_matrix[col_idx][row_idx] = curr_class

                # Check if we've filled up our current class or are at the last
                # iteration (the final distance class may not completely fill
                # up).
                if (i + 1) % class_size == 0 or i == len(order) - 1:
                    curr_class += 1
                    class_indices.append(class_start +
                                         (class_end - class_start) / 2)
                    class_start = class_end

            if curr_class < num_classes:
                # Our last class was empty, so record the last distance seen
                # (which will be the max) as the class index.
                class_indices.append(class_end)

            # Fill diagonal with -1, as it does not belong to any distance
            # class.
            fill_diagonal(dist_class_matrix, -1)
        else:
            # Compute the breakpoints of the distance classes based on the
            # number of specified classes and the ranges of values in the lower
            # triangular portion of the distance matrix (excluding the
            # diagonal).
            break_points = self._find_break_points(np_min(dm_lower_flat),
                                                   np_max(dm_lower_flat),
                                                   num_classes)

            # Find the class indices (the midpoints between breakpoints).
            class_indices = []
            for bp_index, break_point in \
                    enumerate(break_points[0:num_classes]):
                next_bp = break_points[bp_index + 1]
                class_indices.append(break_point +
                                     (0.5 * (next_bp - break_point)))

            # Create the matrix of distance classes. Every element in the
            # matrix tells what distance class the original element belongs to.
            dist_class_matrix = empty([size, size], dtype=int)
            for i in range(size):
                for j in range(size):
                    if i != j:
                        curr_ele = dm[i][j]
                        bps = [(k - 1) for k, bp in enumerate(break_points)
                               if bp >= curr_ele]
                        min_bp = min(bps)

                        # If we somehow got a negative breakpoint (possible
                        # sometimes due to rounding error), put it in the first
                        # distance class.
                        dist_class_matrix[i][j] = min_bp if min_bp >= 0 else 0
                    else:
                        dist_class_matrix[i][j] = -1

        return dist_class_matrix, class_indices

    def _find_row_col_indices(self, idx):
        """Returns row, col for idx into flattened lower triangular matrix.

        It is assumed that the index points to a matrix that was flattened,
        containing only the lower triangular elements (excluding the diagonal)
        in left-to-right, top-to-bottom order (such as that given by
        DistanceMatrix.condensed_form()).
        """
        if idx < 0:
            raise IndexError("The index %d must be greater than or equal to "
                             "zero." % idx)

        # First find the row we're at. The number of elements at each row
        # increases by one each time.
        curr_idx = 0
        delta = 1

        while curr_idx <= idx:
            curr_idx += delta
            delta += 1

        # We subtract one because delta gives us one row past our target.
        row = delta - 1

        # Now that we know the row index, we subtract the number of elements
        # below the row (given by (n*n-n)/2) to find the column that idx is at.
        col = int(idx - ((row * row - row) / 2))

        return row, col

    def _find_break_points(self, start, end, num_classes):
        """Finds the points to break a range into equal width classes.

        Returns a list of floats indicating breakpoints in the range.

        Arguments:
            start - the minimum value in the range
            end - the maximum value in the range
            num_classes - the number of classes to break the range into
        """
        if start >= end:
            raise ValueError("Cannot find breakpoints because the starting "
                             "point is greater than or equal to the ending "
                             "point.")
        if num_classes < 1:
            raise ValueError("Cannot have fewer than one distance class.")

        width = (end - start) / num_classes
        break_points = [start + width * class_num
                        for class_num in range(num_classes)]
        break_points.append(float(end))

        # Move the first breakpoint a little bit to the left. Machine epsilon
        # is taken from:
        # http://en.wikipedia.org/wiki/Machine_epsilon#
        #     Approximation_using_Python
        epsilon = finfo(float).eps
        break_points[0] = break_points[0] - epsilon

        return break_points

    def _correct_p_values(self, p_vals):
        """Corrects p-values for multiple testing using Bonferroni correction.

        This method of correction is non-progressive. If any of the p-values
        are None or NaN, they are not counted towards the number of tests used
        in the correction.

        Returns a list of Bonferroni-corrected p-values for those that are not
        None/NaN. Those that are None/NaN are simply returned. The ordering of
        p-values is maintained.

        Arguments:
            p_vals - list of p-values (of type float or None)
        """
        num_tests = len([p_val for p_val in p_vals
                         if p_val is not None and not isnan(p_val)])

        corrected_p_vals = []
        for p_val in p_vals:
            if p_val is not None and not isnan(p_val):
                corrected_p_vals.append(min(p_val * num_tests, 1))
            else:
                corrected_p_vals.append(p_val)
        return corrected_p_vals

    def _generate_correlogram(self, class_indices, mantel_stats,
                              corrected_p_vals):
        """Generates a matplotlib plot of the Mantel correlogram.

        Returns a matplotlib Figure instance, which can then be manipulated
        further or saved to a file as necessary.

        Arguments:
            class_indices - list of distance class indices (for the x-axis)
            mantel_stats - list of Mantel r stats (for the y-axis)
            corrected_p_vals - list of corrected p-values (for filling in
                points to indicate significance)
        """
        # Plot distance class index versus mantel correlation statistic.
        fig = figure()
        ax = fig.add_subplot(111)
        ax.plot(class_indices, mantel_stats, 'ks-', mfc='white', mew=1)

        # Fill in each point that is significant (based on alpha).
        signif_classes = []
        signif_stats = []
        for idx, p_val in enumerate(corrected_p_vals):
            if p_val is not None and not isnan(p_val) and p_val <= self.Alpha:
                signif_classes.append(class_indices[idx])
                signif_stats.append(mantel_stats[idx])
        ax.plot(signif_classes, signif_stats, 'ks', mfc='k')

        ax.set_title("Mantel Correlogram")
        ax.set_xlabel("Distance class index")
        ax.set_ylabel("Mantel correlation statistic")
        return fig


class PartialMantel(CorrelationStats):

    """Class for the partial Mantel matrix correlation statistical method.

    This class provides the functionality to run a partial Mantel analysis on
    three distance matrices. A partial Mantel test essentially computes the
    Pearson correlation between two distance matrices after first controlling
    for the effects of a third distance matrix (the control matrix).
    """

    def __init__(self, dm1, dm2, cdm):
        """Constructs a new PartialMantel instance.

        Arguments:
            dm1 - first DistanceMatrix object to be compared
            dm2 - second DistanceMatrix object to be compared
            cdm - the control DistanceMatrix object
        """
        super(PartialMantel, self).__init__([dm1, dm2, cdm], num_dms=3,
                                            min_dm_size=3)

    def __call__(self, num_perms=999):
        """Runs a partial Mantel test on the current distance matrices.

        Returns a dict containing the results. The following keys are set:
            method_name - name of the statistical method
            mantel_p - the p-value computed by the test
            mantel_r - the Mantel r statistic computed by the test

        Arguments:
            num_perms - the number of times to permute the distance matrix
                while calculating the p-value

        Credit: The code herein is based loosely on the implementation found in
        R's vegan package.
        """
        res = super(PartialMantel, self).__call__(num_perms)

        # Calculate the correlation statistic.
        corr = lambda rxy, rxz, ryz: (rxy - rxz * ryz) / (sqrt(1 -
                                                               rxz ** 2) * sqrt(1 - ryz ** 2))
        # Load initial/placeholder values in the results dictionary.
        res['method_name'] = 'Partial Mantel'
        res['mantel_r'] = None
        res['mantel_p'] = None

        dm1, dm2, cdm = self.DistanceMatrices
        dm1_flat = dm1.condensed_form()
        dm2_flat = dm2.condensed_form()
        cdm_flat = cdm.condensed_form()

        # Get the initial r-values before permuting.
        rval1 = pearson(dm1_flat, dm2_flat)
        rval2 = pearson(dm1_flat, cdm_flat)
        rval3 = pearson(dm2_flat, cdm_flat)

        # Calculate the original test statistic (r-value).
        orig_stat = corr(rval1, rval2, rval3)

        # Calculate permuted r-values and p-values, storing them for use in the
        # calculation of the final statistic.
        perm_stats = []
        numerator = 0
        for i in range(0, num_perms):
            # Permute the first distance matrix and calculate new r and
            # p-values.
            p1 = permute_2d(dm1, permutation(dm1.shape[0]))
            dm1_perm = DistanceMatrix(p1, dm1.ids)
            dm1_perm_flat = dm1_perm.condensed_form()
            rval1 = pearson(dm1_perm_flat, dm2_flat)
            rval2 = pearson(dm1_perm_flat, cdm_flat)
            perm_stats.append(corr(rval1, rval2, rval3))

            if perm_stats[-1] >= orig_stat:
                numerator += 1
        # Load the final statistics into the result dictionary.
        res['mantel_r'] = orig_stat
        res['mantel_p'] = (numerator + 1) / (num_perms + 1)
        return res


def paired_difference_analyses(personal_ids_to_state_values,
                               analysis_categories,
                               state_values,
                               output_dir,
                               line_color="black",
                               ymin=None,
                               ymax=None):
    """run paired difference analysis one sample t-tests and generate plots

       Apply one-sample t-tests and generate plots to test for changes in
       certain values with a state change. A state change here refers to a
       pre/post-type experimental design, such as pre-treatment to
       post-treatment, and the values that are being tested for change can
       be things like alpha diversity, abundance of specific taxa, a principal
       coordinate value (e.g., PC1 value before and after treatment), and so
       on.

       The one-sample t-test is applied on each pair of differences. So, if
       experiment was based on looking for changes in proteobacteria abundance
       with treatment, you would have pre- and post-treatment proteobacteria
       abundances for a number of individuals. The difference would be computed
       between those, and the null hypothesis is that the mean of those differences
       is equal to zero (i.e., no change with treatment).

       Line plots are also generated to show the change on a per-individual basis.

     personal_ids_to_state_values: a 2d dictionary mapping personal ids to potential
      analysis categories, which each contain a pre/post value. this might look like
      the following:
       {'firmicutes-abundance':
            {'subject1':[0.45,0.55],
             'subject2':[0.11,0.52]},
           'bacteroidetes-abundace':
             {'subject1':[0.28,0.21],
              'subject2':[0.11,0.01]}
        }
       examples of functions that can be useful for generating these data are
        qiime.parse.extract_per_individual_state_metadata_from_sample_metadata and
        qiime.parse.extract_per_individual_state_metadata_from_sample_metadata_and_biom

     analysis_categories: a list of categories to include in analyses (e.g,
       ['firmicutes-abundance', 'bacteroidetes-abundace'])

     state_values: an ordered list describing each of the states being compared (these
       are the x labels in the resulting plots)

     output_dir: directory where output should be written (will be created if
       it doesn't exist)

     ymin: minimum y-value in plots (if it should be consistent across
       plots - by default will be chosen on a per-plot basis)

     ymax: maximum y-value in plots (if it should be consistent across
       plots - by default will be chosen on a per-plot basis)
    """

    if len(state_values) != 2:
        raise ValueError("Only two state values can be provided. "
                         "Support currently exists only for pre/post experimental design.")

    # create the output directory if it doesn't already exist
    create_dir(output_dir)

    num_analysis_categories = len(analysis_categories)
    x_values = range(len(state_values))

    paired_difference_output_fp = \
        join(output_dir, 'paired_difference_comparisons.txt')
    paired_difference_output_f = open(paired_difference_output_fp, 'w')
    # write header line to output file
    paired_difference_output_f.write(
        "#Metadata category\tNum differences (i.e., n)\tMean difference\t"
        "Median difference\tt one sample\tt one sample parametric p-value\t"
        "t one sample parametric p-value (Bonferroni-corrected)\n")

    paired_difference_t_test_results = {}

    biom_table_fp = join(output_dir, 'differences.biom')
    biom_sids_fp = join(output_dir, 'differences_sids.txt')
    biom_observation_ids = []
    biom_data = []
    # need a list of personal_ids to build the biom table -
    # ugly, but get it working first
    personal_ids = []
    for c in personal_ids_to_state_values.values():
        personal_ids.extend(c.keys())
    personal_ids = list(set(personal_ids))

    # initiate list of output file paths to return
    output_fps = [paired_difference_output_fp,
                  biom_table_fp,
                  biom_sids_fp]

    num_successful_tests = 0
    included_personal_ids = defaultdict(list)
    for category_number, analysis_category in enumerate(analysis_categories):
        personal_ids_to_state_metadatum = personal_ids_to_state_values[
            analysis_category]
        analysis_category_fn_label = analysis_category.replace(' ', '-')
        plot_output_fp = join(
            output_dir,
            '%s.pdf' %
            analysis_category_fn_label)
        fig = figure()
        axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])

        # initialize a list to store the distribution of changes
        # with state change
        differences = []
        pre_values = []
        post_values = []
        store_biom_datum = True

        for personal_id in personal_ids:
            data = personal_ids_to_state_metadatum[personal_id]
            if None in data:
                # if any of the data points are missing, don't store
                # a difference for this individual, and don't store
                # the category in the BIOM table
                store_biom_datum = False
                raise ValueError("Some data points are missing, "
                                 "cannot create biom file.")
            else:
                # otherwise compute the difference between the ending
                # and starting state
                pre_value = data[0]
                post_value = data[1]
                included_personal_ids[personal_id].append(pre_value)
                included_personal_ids[personal_id].append(post_value)
                pre_values.append(pre_value)
                post_values.append(post_value)
                difference = post_value - pre_value
                differences.append(difference)
                # and plot the start and stop values as a line
                axes.plot(x_values, data, line_color, linewidth=0.5)

        if store_biom_datum:
            biom_observation_ids.append(analysis_category)
            biom_data.append(differences)

        # run stats for current analysis category
        n = len(differences)
        mean_differences = mean(differences)
        median_differences = median(differences)
        t_one_sample_results = t_one_sample(differences)
        t = t_one_sample_results[0]
        p_value = t_one_sample_results[1]
        if p_value is not None:
            num_successful_tests += 1
        # analysis_category gets stored as the key and the first entry
        # in the value to faciliate sorting the values and writing to
        # file
        paired_difference_t_test_results[analysis_category] = \
            [analysis_category,
             n,
             mean_differences,
             median_differences,
             t,
             p_value]

        # Finalize plot for current analysis category
        axes.plot(x_values,
                  [median(pre_values), median(post_values)],
                  line_color,
                  linewidth=3,
                  ls='--')
        axes.set_ylabel(analysis_category)
        axes.set_xticks(range(len(state_values)))
        axes.set_xticklabels(state_values)
        axes.set_ylim(ymin=ymin, ymax=ymax)
        fig.savefig(plot_output_fp, transparent=True)
        output_fps.append(plot_output_fp)

    # write a biom table based on differences and
    # a list of the sample ids that could be converted
    # to a mapping file for working with this biom table

    biom_table = Table(biom_data,
                               biom_observation_ids,
                               personal_ids,
                               input_is_dense=True)
    write_biom_table(biom_table, biom_table_fp)
    biom_sids_f = open(biom_sids_fp, 'w')
    sid_headers = ['#SampleID']
    for e in analysis_categories:
        sid_headers.append('Pre-%s' % e)
        sid_headers.append('Post-%s' % e)
    biom_sids_f.write('%s\n' % ('\t'.join(sid_headers)))
    for sid, data in included_personal_ids.iteritems():
        data_str = '\t'.join(map(str,data))
        biom_sids_f.write('%s\t%s\n' % (sid, data_str))
    biom_sids_f.close()

    # sort stats output by uncorrected p-value, compute corrected p-value,
    # and write results to file
    paired_difference_t_test_lines = \
        paired_difference_t_test_results.values()
    paired_difference_t_test_lines.sort(key=lambda x: x[5])
    for r in paired_difference_t_test_lines:
        p_value = r[5]
        if p_value is None:
            bonferroni_p_value = None
        else:
            bonferroni_p_value = min([p_value * num_successful_tests, 1.0])
        paired_difference_output_f.write(
            '\t'.join(map(str, r + [bonferroni_p_value])))
        paired_difference_output_f.write('\n')
    paired_difference_output_f.close()

    return output_fps, paired_difference_t_test_results


class ZeroExpectedError(ValueError):

    """Class for handling tests where an expected value was zero."""
    pass


def G_2_by_2(a, b, c, d, williams=1, directional=1):
    """G test for independence in a 2 x 2 table.

    Usage: G, prob = G_2_by_2(a, b, c, d, willliams, directional)

    Cells are in the order:

        a b
        c d

    a, b, c, and d can be int, float, or long.
    williams is a boolean stating whether to do the Williams correction.
    directional is a boolean stating whether the test is 1-tailed.

    Briefly, computes sum(f ln f) for cells - sum(f ln f) for
    rows and columns + f ln f for the table.

    Always has 1 degree of freedom

    To generalize the test to r x c, use the same protocol:
    2*(cells - rows/cols + table), then with (r-1)(c-1) df.

    Note that G is always positive: to get a directional test,
    the appropriate ratio (e.g. a/b > c/d) must be tested
    as a separate procedure. Find the probability for the
    observed G, and then either halve or halve and subtract from
    one depending on whether the directional prediction was
    upheld.

    The default test is now one-tailed (Rob Knight 4/21/03).

    See Sokal & Rohlf (1995), ch. 17. Specifically, see box 17.6 (p731).
    """
    cells = [a, b, c, d]
    n = sum(cells)
    # return 0 if table was empty
    if not n:
        return (0, 1)
    # raise error if any counts were negative
    if min(cells) < 0:
        raise ValueError(
            "G_2_by_2 got negative cell counts(s): must all be >= 0.")

    G = 0
    # Add x ln x for items, adding zero for items whose counts are zero
    for i in filter(None, cells):
        G += i * log(i)
    # Find totals for rows and cols
    ab = a + b
    cd = c + d
    ac = a + c
    bd = b + d
    rows_cols = [ab, cd, ac, bd]
    # exit if we are missing a row or column entirely: result counts as
    # never significant
    if min(rows_cols) == 0:
        return (0, 1)
    # Subtract x ln x for rows and cols
    for i in filter(None, rows_cols):
        G -= i * log(i)
    # Add x ln x for table
    G += n * log(n)
    # Result needs to be multiplied by 2
    G *= 2

    # apply Williams correction
    if williams:
        q = 1 + \
            ((((n / ab) + (n / cd)) - 1) * (((n / ac) + (n / bd)) - 1)) / \
            (6 * n)
        G /= q

    p = chi2prob(G, 1, direction='high')

    # find which tail we were in if the test was directional
    if directional:
        is_high = ((b == 0) or (d != 0 and (a / b > c / d)))
        p = tail(p, is_high)
        if not is_high:
            G = -1 * G
    return G, p


def safe_sum_p_log_p(a, base=None):
    """Calculates p * log(p) safely for an array that may contain zeros."""
    flat = ravel(a)
    nz = take(flat, nonzero(flat)[0])
    logs = log(nz)
    if base:
        logs /= log(base)
    return sum(nz * logs, 0)


def g_fit(data, williams=True):
    """Calculate the G statistic aka log-likelihood ratio test.

    Parameters
    ----------
    data : iterable of 1-D array_like
        Each element of the iterable is 1D with any length and represents the
        observed frequencies of a given OTU in one of the sample classes.
    williams : boolean
        Whether or not to apply the Williams correction before comparing to the
        chi-squared distribution.

    Returns
    -------
    G : float
        The G statistic that is additive between all groups.
    pval : float
        The pvalue associated with the given G statistic.

    Notes
    -----
    For discussion read [1]_. This function compares the calculated G statistic
    (with Williams correction) to the chi-squared distribution with the
    appropriate number of degrees of freedom. If the data do not pass sanity
    checks for basic assumptions then this function will return nans.

    This function wraps the scipy function scipy.stats.power_divergence. When
    comparing this to the original implementation. Unfortunately, scipy does
    not have the Williams correction, so this script adds this as an option.
    For discussion read [1]_ pg. 695-699.

    The G testis normally applied to data when you have only one observation of
    any given sample class (e.g. you observe 90 wildtype and 30 mutants). In
    microbial ecology it is normal to have multiple samples which contain a
    given feature where those samples share a metadata class (eg. you observe
    OTUX at certain frequencies in 12 samples, 6 of which are treatment
    samples, 6 of which are control samples). To reconcile these approaches
    this function averages the frequency of the given feature (OTU) across all
    samples in the metadata class (e.g. in the 6 treatment samples, the value
    for OTUX is averaged, and this forms the average frequency which represents
    all treatment samples in aggregate). This means that this version of the G
    stat cannot detect sample heterogeneity as a replicated goodness of fit
    test would be able to. In addition, this function assumes the extrinsic
    hypothesis is that the mean frequency in all the samples groups is the
    same.

    References
    ----------
    .. [1] Sokal and Rohlf. "Biometry: The Principles and Practices of
       Statistics in Biological Research". ISBN: 978-0716724117
    """
    r_data = [array(i).mean() for i in data]
    G, p = power_divergence(r_data, lambda_="log-likelihood")
    if williams:
        G_corr = williams_correction(sum(r_data), len(r_data), G)
        return G_corr, chi2prob(G_corr, len(r_data) - 1, direction='high')
    else:
        return G, p


def williams_correction(n, a, G):
    """Return the Williams corrected G statistic for G goodness of fit test.

    For discussion read [1]_ pg 698-699.

    Parameters
    ----------
    n : int
        Sum of observed frequencies.
    a : int
        Number of groups that are being compared.
    G : float
        Uncorrected G statistic

    Notes
    -----
    The equation given in this function is simplified from [1]_
    q = 1. + (a**2 - 1)/(6.*n*a - 6.*n) == 1. + (a+1.)/(6.*n)

    References
    ----------
    .. [1] Sokal and Rohlf. "Biometry: The Principles and Practices of
       Statistics in Biological Research". ISBN: 978-0716724117
    """
    q = 1. + (a + 1.) / (6. * n)
    return G / q


def t_paired(a, b, tails='two-sided', exp_diff=0):
    """Returns t and prob for TWO RELATED samples of scores a and b.

    From Sokal and Rohlf (1995), p. 354.
    Calculates the vector of differences and compares it to exp_diff
    using the 1-sample t test.

    Usage:   t, prob = t_paired(a, b, tails, exp_diff)

    t is a float; prob is a probability.
    a and b should be equal-length lists of paired observations (numbers).
    tails should be None (default), 'high', or 'low'.
    exp_diff should be the expected difference in means (a-b); 0 by default.
    """
    if len(a) != len(b):
        raise ValueError('Unequal length lists in ttest_paired.')
    return t_one_sample(array(a) - array(b), popmean=exp_diff, tails=tails)


def t_one_sample(a, popmean=0, tails='two-sided'):
    '''Peform a one sample t-test against a given population mean.

    Parameters
    ----------
    a : array-like
        A vector of observations.
    popmean : float
        The population mean to test against.
    tails : str
        The hypothesis to test, one of 'low', 'high', 'two-sided'.

    Returns
    -------
    t : float
        t statstic.
    p : float
        p-value assocaited with the t-statistic given the tails.
    '''
    t, _ = ttest_1samp(a, popmean)  # returns array([t]), p
    if isnan(t) or isinf(t):
        return nan, nan
    p = tprob(t, len(a) - 1, tails)
    return float(t), p  # cast t to a float rather than 0-d array


def t_two_sample(a, b, tails='two-sided', exp_diff=0):
    '''scipy t_two_sample.'''
    if len(a) == 1 or len(b) == 1: #need to use t_one_observation
        if len(a) <= len(b):
            t, p = t_one_observation(a, b, tails, exp_diff)
        else:
            t, p = t_one_observation(b, a, tails, exp_diff)
        return t, p
    t, _ = ttest_ind(asarray(a) - exp_diff, asarray(b), axis=0, equal_var=True)
    if isnan(t) or isinf(t):
        return nan, nan
    p = tprob(t, len(a) + len(b) - 2., tails)
    return float(t), p


def mc_t_two_sample(x_items, y_items, tails='two-sided', permutations=999,
                    exp_diff=0):
    """Performs a two-sample t-test with Monte Carlo permutations.

    x_items and y_items must be INDEPENDENT observations (sequences of
    numbers). They do not need to be of equal length.

    Returns the observed t statistic, the parametric p-value, a list of t
    statistics obtained through Monte Carlo permutations, and the nonparametric
    p-value obtained from the Monte Carlo permutations test.

    This code is partially based on Jeremy Widmann's
    qiime.make_distance_histograms.monte_carlo_group_distances code.

    Arguments:
        x_items - the first list of observations
        y_items - the second list of observations
        tails - if None (the default), a two-sided test is performed. 'high'
            or 'low' for one-tailed tests
        permutations - the number of permutations to use in calculating the
            nonparametric p-value. Must be a number greater than or equal to 0.
            If 0, the nonparametric test will not be performed. In this case,
            the list of t statistics obtained from permutations will be empty,
            and the nonparametric p-value will be NaN
        exp_diff - the expected difference in means (x_items - y_items)
    """
    if permutations < 0:
        raise ValueError("Invalid number of permutations: %d. Must be greater "
                         "than or equal to zero." % permutations)

    if (len(x_items) == 1 and len(y_items) == 1) or \
       (len(x_items) < 1 or len(y_items) < 1):
        raise ValueError("At least one of the sequences of observations is "
                         "empty, or the sequences each contain only a single "
                         "observation. Cannot perform the t-test.")

    # Perform t-test using original observations.
    obs_t, param_p_val = t_two_sample(x_items, y_items, tails=tails,
                                      exp_diff=exp_diff)

    # Only perform the Monte Carlo test if we got a sane answer back from the
    # initial t-test and we have been specified permutations.
    nonparam_p_val = nan
    perm_t_stats = []
    if permutations > 0 and not isnan(obs_t) and not isnan(param_p_val):
        perm_t_stats = zeros(permutations, dtype=float)
        px, py = _permute_observations(x_items, y_items, permutations)
        for i in range(permutations):
            perm_t_stats[i] = t_two_sample(px[i], py[i], tails=tails,
                                           exp_diff=exp_diff)[0]

        # Compute nonparametric p-value based on the permuted t-test results.
        if tails == 'two-sided':
            better = (abs(perm_t_stats) >= abs(obs_t)).sum()
        elif tails == 'low':
            better = ((perm_t_stats) <= obs_t).sum()
        elif tails == 'high':
            better = ((perm_t_stats) >= obs_t).sum()
        nonparam_p_val = (better + 1) / (permutations + 1)

    return obs_t, param_p_val, perm_t_stats, nonparam_p_val


def _permute_observations(x, y, num_perms):
    """Return num_perms pairs of permuted vectors x,y.

    Parameters
    ----------
    x : 1-D array-like
        Lists or arrays of values to be permuted.
    y : 1-D array-like
        Lists or arrays of values to be permuted.

    Returns
    -------
    xs : list of arrays
        Permuted vectors x
    ys : list of arrays
        Permuted vectors y
    """
    vals = hstack([array(x), array(y)])
    lenx = len(x)
    # sorting step is unnecessary for this code, but it ensure that test code
    # which relies on seeding the prng works (if we dont do this then different
    # observation orders in x and y for eg. the mc_t_two_sample test will fail
    # to produce the same results)
    vals.sort()
    inds = arange(vals.size)
    xs, ys = [], []
    for i in range(num_perms):
        shuffle(inds)
        xs.append(vals[inds[:lenx]])
        ys.append(vals[inds[lenx:]])
    return xs, ys


def t_one_observation(x, sample, tails='two-sided', exp_diff=0):
    """Returns t-test for significance of single observation versus a sample.

    Parameters
    ----------
    x : float
        The single observation to test against the sample.
    sample : array-like
        Vector of observations for the sample to test against x.
    tails : str
        The hypothesis to test, one of 'low', 'high', 'two-sided'.
    exp_diff : float
        The expected difference between the sample mean and the observation.

    Returns
    -------
    t : float
        t statstic.
    p : float
        p-value assocaited with the t-statistic given the tails.

    Notes
    -----
    Equation for 1-observation t [1]_ p 228:
    t = obs - mean - exp_diff / (var * sqrt((n+1)/n))
    df = n - 1

    References
    ----------
    .. [1] Sokal and Rohlf. "Biometry: The Principles and Practices of
       Statistics in Biological Research". ISBN: 978-0716724117

    """
    try:
        sample_mean = mean(sample)
        sample_std = std(sample, ddof=1)

        if sample_std == 0:  # no variance means can't compute t, p
            return (nan, nan)

        else:  # The list varies.
            n = len(sample)
            t = ((x - sample_mean - exp_diff) / sample_std / sqrt((n + 1) / n))
            prob = tprob(t, n - 1, tails)
            return (float(t), prob)

    except (ZeroDivisionError, ValueError, AttributeError, TypeError,
            FloatingPointError):
        return (nan, nan)


def pearson(v1, v2):
    '''Pearson correlation using numpy.corrcoef.

    Parameters
    ----------
    v1 : array-like
        List or array of ints or floats to be correlated.
    v2 : array-like
        List or array of ints or floats to be correlated.

    Returns
    -------
    corrcoef : float
        Pearson correlation between the vectors.

    Raises
    ------
    ValueError
        If the vectors are not equally sized or if they are only a single
        element a ValueError will be returned.

    Examples
    --------
    >>> from qiime.stats import pearson
    >>> v1 = [.1, .2, .5, .3, .4]
    >>> v2 = [.9, .01, .5, .6, .7]
    >>> pearson(v1, v2)
    -0.052364331421504685
    '''
    v1, v2 = array(v1), array(v2)
    if not (v1.size == v2.size > 1):
        raise ValueError('One or more vectors isn\'t long enough to correlate '
                         ' or they have unequal lengths.')
    return corrcoef(v1, v2)[0][1]  # 2x2 symmetric unit matrix


def spearman(v1, v2):
    """Returns Spearman's rho.

    Parameters
    ----------
    v1 : array-like
        List or array of ints or floats to be correlated.
    v2 : array-like
        List or array of ints or floats to be correlated.

    Returns
    -------
    rho : float
        Spearman correlation between the vectors.

    Raises
    ------
    ValueError
        If the vectors are not equally sized or if they are only a single
        element a ValueError will be returned.

    See Also
    --------
    scipy.stats.spearmanr

    Notes
    -----
    This will always be a value between -1.0 and +1.0. v1 and v2 must
    be the same length, and cannot have fewer than 2 elements each. If one or
    both of the input vectors do not have any variation, the return value will
    be nan.
    """
    v1, v2 = array(v1), array(v2)
    if not (v1.size == v2.size > 1):
        raise ValueError('One or more vectors isn\'t long enough to correlate '
                         ' or they have unequal lengths.')
    return spearmanr(v1, v2)[0]  # return only the rho-value


def kendall(v1, v2):
    """Compute Kendall's Tau between v1 and v2 using scipy.stats.kendalltau

    Parameters
    ----------
    v1 : array-like
        List or array of ints or floats to be correlated.
    v2 : array-like
        List or array of ints or floats to be correlated.

    Returns
    -------
    rho : float
        Spearman correlation between the vectors.

    Raises
    ------
    ValueError
        If the vectors are not equally sized or if they are only a single
        element a ValueError will be returned.
    """
    v1, v2 = array(v1), array(v2)
    if not (v1.size == v2.size > 1):
        raise ValueError('One or more vectors isn\'t long enough to correlate '
                         ' or they have unequal lengths.')
    return kendalltau(v1, v2)[0]  # return only the tau correlation coeff


def kendall_pval(tau, n):
    '''Calculate the p-value for the passed tau and vector length n.'''
    test_stat = tau / ((2 * (2 * n + 5)) / float(9 * n * (n - 1))) ** .5
    return normprob(test_stat, direction='two-sided')


def assign_correlation_pval(corr, n, method, permutations=None,
                            perm_test_fn=None, v1=None, v2=None):
    """Assign pval to a correlation score with given method.

    This function will assign significance to the correlation score passed
    given the method that is passed. Some of the methods are appropriate only
    for certain types of data and there is no way for this test to determine
    the appropriateness, thus you must use this function only when with the
    proper prior knowledge. The 'parametric_t_distribution' method is described
    in [1]_ pg. 576, the 'fisher_z_transform' method
    is described on pg 576 and 577. The 'bootstrap' method calculates the given
    correlation permutations number of times using perm_test_fn.
    Also note, this does *not* take the place of FDR correction.

    Paramters
    ---------
    corr : float
        Correlation score from Kendall's Tau, Spearman's Rho, or Pearson.
    n : int
        Length of the vectors that were correlated.
    method : str
        One of ['parametric_t_distribution', 'fisher_z_transform',
        'bootstrapped', 'kendall'].
    permutations : int
        Number of permutations to use if bootstrapped selected.
    perm_test_fn : function
        Used to use to calculate correlation if permuation test desired.
    v1 : array-like or None
        List or array of ints or floats to be correlated. Passed if
        method='bootstrapped'.
    v2 : array-like or None
        List or array of ints or floats to be correlated. Passed if
        method='bootstrapped'

    References
    ----------
    .. [1] Sokal and Rohlf. "Biometry: The Principles and Practices of
       Statistics in Biological Research". ISBN: 978-0716724117
    """
    if method == 'parametric_t_distribution':
        df = n - 2
        if df <= 1:
            raise ValueError("Must have more than 1 degree of freedom. "
                             "Can't Continue.")
        try:
            ts = corr * ((df / (1. - corr ** 2)) ** .5)
            return tprob(ts, df, tails='two-sided')
            # two tailed test because H0 is corr=0
        except (ValueError, FloatingPointError, ZeroDivisionError):
            # something unpleasant happened, most likely r or rho where +- 1
            # which means the parametric p val should be 1 or 0 or nan
            return nan
    elif method == 'fisher_z_transform':
        # Sokal and Rohlf indicate that for n<50, the Fisher Z transform for
        # assigning correlation probabilities is not accurate. Currently no
        # check is in place
        z = fisher_z_transform(corr)
        # the z transform pval compares against a t distribution with inf
        # degrees of freedom which is equal to a z distribution.
        return z_transform_pval(z, n)
    elif method == 'bootstrapped':
        if any([i is None for i in [v1, v2, permutations, perm_test_fn]]):
            raise ValueError('You must specify vectors, permutation '
                             'function, and number of permutations to calc '
                             'bootstrapped pvalues. Cant continue.')
        if any([isnan(corr), isinf(corr)]):
            return nan
        else:
            r = empty(permutations)
            for i in range(permutations):
                r[i] = perm_test_fn(v1, permutation(v2))
            return (abs(r) >= abs(corr)).sum() / float(permutations)
    elif method == 'kendall':
        return kendall_pval(corr, n)
    else:
        raise ValueError("'%s' method is unknown." % method)


def correlation_t(x_items, y_items, method='pearson', tails='two-sided',
                  permutations=999, confidence_level=0.95):
    """Computes the correlation between two vectors and its significance.

    Computes a parametric p-value by using Student's t-distribution with df=n-2
    to perform the test of significance, as well as a nonparametric p-value
    obtained by permuting one of the input vectors the specified number of
    times given by the permutations parameter. A confidence interval is also
    computed using Fisher's Z transform if the number of observations is
    greater than 3. Please see Sokal and Rohlf pp. 575-580 and pg. 598-601 for
    more details regarding these techniques.

    Warning: the parametric p-value is unreliable when the method is spearman
    and there are less than 11 observations in each vector.

    Returns the correlation coefficient (r or rho), the parametric p-value, a
    list of the r or rho values obtained from permuting the input, the
    nonparametric p-value, and a tuple for the confidence interval, with the
    first element being the lower bound of the confidence interval and the
    second element being the upper bound for the confidence interval. The
    confidence interval will be (None, None) if the number of observations is
    not greater than 3.

    x_items and y_items must be the same length, and cannot have fewer than 2
    elements each. If one or both of the input vectors do not have any
    variation, r or rho will be 0.0.

    Note: the parametric portion of this function is based on the correlation
    function in this module.

    Arguments:
        x_items - the first list of observations
        y_items - the second list of observations
        method - 'pearson' or 'spearman'
        tails - if None (the default), a two-sided test is performed. 'high'
            for a one-tailed test for positive association, or 'low' for a
            one-tailed test for negative association. This parameter affects
            both the parametric and nonparametric tests, but the confidence
            interval will always be two-sided
        permutations - the number of permutations to use in the nonparametric
            test. Must be a number greater than or equal to 0. If 0, the
            nonparametric test will not be performed. In this case, the list of
            correlation coefficients obtained from permutations will be empty,
            and the nonparametric p-value will be None
        confidence_level - the confidence level to use when constructing the
            confidence interval. Must be between 0 and 1 (exclusive)
    """
    # Perform some initial error checking.
    if method == 'pearson':
        corr_fn = pearson
    elif method == 'spearman':
        corr_fn = spearman
    else:
        raise ValueError("Invalid method '%s'. Must be either 'pearson' or "
                         "'spearman'." % method)
    if permutations < 0:
        raise ValueError("Invalid number of permutations: %d. Must be greater "
                         "than or equal to zero." % permutations)
    if confidence_level <= 0 or confidence_level >= 1:
        raise ValueError("Invalid confidence level: %.4f. Must be between "
                         "zero and one." % confidence_level)

    # Calculate the correlation coefficient.
    corr_coeff = corr_fn(x_items, y_items)

    # Perform the parametric test first.
    x_items, y_items = array(x_items), array(y_items)
    n = len(x_items)
    df = n - 2
    if n < 3:
        parametric_p_val = 1
    else:
        try:
            t = corr_coeff / sqrt((1 - (corr_coeff * corr_coeff)) / df)
            parametric_p_val = tprob(t, df, tails)
        except (ZeroDivisionError, FloatingPointError):
            # r/rho was presumably 1.
            parametric_p_val = 0

    # Perform the nonparametric test.
    perm_ccs = zeros(permutations, dtype=float)
    nonparametric_p_val = None
    better = 0
    for i in range(permutations):
        perm_ccs[i] = corr_fn(x_items, y_items[permutation(n)])

    if tails == 'two-sided':
        better = (abs(perm_ccs.round(15)) >= abs(round(corr_coeff, 15))).sum()
    elif tails == 'high':
        better = (perm_ccs.round(15) >= round(corr_coeff, 15)).sum()
    elif tails == 'low':
        better = (perm_ccs.round(15) <= round(corr_coeff, 15)).sum()
    else:
        # Not strictly necessary since this was checked above, but included
        # for safety in case the above check gets removed or messed up. We
        # don't want to return a p-value of 0 if someone passes in a bogus
        # tail type somehow.
        raise ValueError("Invalid tail type '%s'. Must be either None, "
                         "'high', or 'low'." % tails)
    if permutations > 0:
        nonparametric_p_val = (better + 1) / (permutations + 1)

    # Compute the confidence interval for corr_coeff using Fisher's Z
    # transform.
    z_crit = abs(ndtri((1 - confidence_level) / 2))
    ci_low, ci_high = None, None

    if n > 3:
        try:
            ci_low = tanh(arctanh(corr_coeff) - (z_crit /
                                                 sqrt(n - 3)))
            ci_high = tanh(arctanh(corr_coeff) + (z_crit /
                                                  sqrt(n - 3)))
        except (ZeroDivisionError, FloatingPointError):
            # r/rho was presumably 1 or -1. Match what R does in this case.
            ci_low, ci_high = corr_coeff, corr_coeff

    return (corr_coeff, parametric_p_val, perm_ccs,
            nonparametric_p_val, (ci_low, ci_high))


def fisher(probs):
    """Uses Fisher's method to combine multiple tests of a hypothesis.

    -2 * SUM(ln(P)) gives chi-squared distribution with 2n degrees of freedom.
    """
    try:
        return chi2prob(-2 * sum(log(probs)), 2 * len(probs), direction='high')
    except OverflowError:
        return 0.0


def ANOVA_one_way(a):
    """Performs a one way analysis of variance

    a is a list of lists of observed values. Each list is the values
    within a category. The analysis must include 2 or more categories(lists).
    Each category of the list, and overall list, is converted to a numpy array.

    An F value is first calculated as the variance of the group means
    divided by the mean of the within-group variances.
    """
    group_means = []
    group_variances = []
    num_cases = 0  # total observations in all groups
    all_vals = []
    for i in a:
        num_cases += len(i)
        group_means.append(mean(i))
        group_variances.append(i.var(ddof=1) * (len(i) - 1))
        all_vals.extend(i)

    # Get within Group variances (denominator)
    dfd = num_cases - len(group_means)
    # need to add a check -- if the sum of the group variances is zero it will
    # error, but only if the between_Groups value is not zero
    within_Groups = sum(group_variances) / dfd
    if within_Groups == 0.:
        return nan, nan
    # Get between Group variances (numerator)
    all_vals = array(all_vals)
    grand_mean = all_vals.mean()
    between_Groups = 0
    for i in a:
        diff = i.mean() - grand_mean
        diff_sq = diff * diff
        x = diff_sq * len(i)
        between_Groups += x

    dfn = len(group_means) - 1
    between_Groups = between_Groups / dfn
    F = between_Groups / within_Groups
    return F, fprob(F, dfn, dfd, direction='high')


def _average_rank(start_rank, end_rank):
    ave_rank = sum(range(start_rank, end_rank + 1)) / \
        (1 + end_rank - start_rank)
    return ave_rank


def _get_bootstrap_sample(x, y, num_reps):
    """yields num_reps random samples drawn with replacement from x and y"""
    combined = hstack([x, y])
    total_obs = len(combined)
    num_x = len(x)
    for i in range(num_reps):
        # sampling with replacement
        indices = randint(0, total_obs, total_obs)
        sampled = combined.take(indices)
        # split into the two populations
        sampled_x = sampled[:num_x]
        sampled_y = sampled[num_x:]
        yield sampled_x, sampled_y


def mw_t(x, y, continuity=True, two_sided=True):
    '''Compute the Mann Whitney U statistic using scipy.stats.mannwhitneyu

    This wrapper controls whether the continuity correction will be applied
    and whether or not a two sided hypothesis is specified.

    Parameters
    ----------
    x : array-like
        List or array of numeric values to be tested.
    y : array-like
        List or array of numeric values to be tested.
    continuity : boolean
        Whether or not to use the continuity correction.
    two_sided: boolean
        Whether or not to use a two sided test. See Notes.

    Returns
    -------
    U stat : float
        The MWU U statistic.
    p-value : float
        The pvalue associated with the given U statistic assuming a normal
        probability distribution.

    See Also
    --------
    scipy.stats.mannwhitneyu

    Notes
    -----
    Two tails is appropriate because we do not know which of our groups has a
    higher mean, thus our alternate hypothesis is that the distributions from
    which the two samples come are not the same (FA!=FB) and we must account
    for E[FA] > E[FB] and E[FB] < E[FA]. See [1]_ pgs 427-431.

    References
    ----------
    .. [1] Sokal and Rohlf. "Biometry: The Principles and Practices of
       Statistics in Biological Research". ISBN: 978-0716724117
    '''
    u, pval = mannwhitneyu(x, y, continuity)
    if two_sided:
        return u, 2. * pval
    else:
        return u, pval


def mw_boot(x, y, num_reps=999):
    """Bootstrapped version of Mann-Whitney-U test

    Parameters
    ----------
    x : array-like
        List or array of numeric values to be tested.
    y : array-like
        List or array of numeric values to be tested.
    num_reps : int
        Number of permutations tests to do.

    Returns
    -------
    observed_stat : float
        Value of the U statistic for the comparison of x and y.
    pval : float
        Number of times a U statistic as small or smaller than the observed U
        statistic was found.

    Notes
    -----
    The u statistic must be smaller than the observed u statistic to count as
    more extreme according to [1]_. Only a two tailed test is allowed through
    this function.

    Examples
    --------
    >>> from qiime.stats import mw_boot
    >>> x = [1.5, 4.6, 7.8, 10.2, 23.4]
    >>> y = [3.4, 10.1, 100.3, 45.6, 45.6, 78.9]
    >>> mw_boot(x, y, num_reps = 999)
    (6.0, 0.079)

    References
    ----------
    .. [1] http://docs.scipy.org/doc/scipy-0.13.0/reference/generated/scipy.sta
       ts.mannwhitneyu.html
    """
    tol = MACHEP * 100
    observed_stat, obs_p = mw_t(x, y)
    u_stats_as_or_more_extreme = 0
    for sampled_x, sampled_y in _get_bootstrap_sample(x, y, num_reps):
        try:
            sample_stat, sample_p = mw_t(sampled_x, sampled_y)
            if sample_stat <= (observed_stat - tol):
                # the u statistic must be smaller than the observed u statistic
                # to count as more extreme. see [1]
                u_stats_as_or_more_extreme += 1
        except ValueError:  # the mwu test got identical x,y items
            pass  # we don't add to the u stats, this was not more extreme
    return observed_stat, (u_stats_as_or_more_extreme + 1) / (num_reps + 1)


def kruskal_wallis(data):
    '''Calculate Kruskal Wallis U stat and pval using scipy.stats.kruskal

    Parameters
    ----------
    data : list of array-likes
        data is a nested list whose elements are arrays of float data. The
        different lists correpsond to the groups being tested.

    Returns
    -------
    U stat : float
        The Kruskal Wallis U statistic.
    pval : float
        The pvalue associated with the given U statistic assuming a chi-squared
        probability distribution.

    Examples
    --------
    >>> from qiime.stats import kruskal_wallis
    >>> data = [[1, 4.5, 67, 100, 2], [145, 100, 3, 14.5, -19], [2, 1.1, 5.5,
    ...         3.3, 16.7, 18, 100.3]]
    >>> rho, pval = kruskal_wallis(data)
    >>> print rho == 0.16848016848016789
    True
    >>> print pval == 0.91921054163678728
    True
    '''
    return kruskal(*data)


def permute_2d(m, p):
    """Performs 2D permutation of matrix m according to p."""
    return m[p][:, p]


def is_symmetric_and_hollow(matrix):
    """Return True if matrix is symmetric and hollow, otherwise False."""
    return (matrix.T == matrix).all() and (trace(matrix) == 0)


def tail(prob, test):
    """If test is true, returns prob/2. Otherwise returns 1-(prob/2).
    """
    prob /= 2
    if test:
        return prob
    else:
        return 1 - prob


def bonferroni_correction(pvals):
    """Adjust pvalues for multiple tests using the Bonferroni method.

    In short: multiply all pvals by the number of comparisons.

    Parameters
    ----------
    pvals : list or array
        List or array of floats.

    Returns
    -------
    list of pvals
        Returns the list of pvals multiplied by their length. Pvals are
        still unsorted (i.e. order has not changed).

    See Also
    --------
    benjamini_hochberg_step_down

    Examples
    --------
    >>> from qiime.stats import bonferroni_correction
    >>> bonferroni_correction([0.1, 0.21, 0.5, 0.2, 0.6])
    array([ 0.5 ,  1.05,  2.5 ,  1.  ,  3.  ])
    """
    return array(pvals, dtype=float) * len(pvals)  # float conv: Nones->nans


def fdr_correction(pvals):
    """Adjust pvalues for multiple tests using the false discovery rate method.

    Parameters
    ----------
    pvals : list or array
        List or array of floats.

    Returns
    -------
    list of pvals
        Returns the list of pvals properly adjusted based on the FDR. Pvals are
        still unsorted (i.e. order has not changed).

    See Also
    --------
    benjamini_hochberg_step_down

    Notes
    -----
    In short: ranks the p-values in ascending order and multiplies each p-value
    by the number of comparisons divided by the rank of the p-value in the
    sorted list. Input is list of floats.  Does *not* assume pvals is sorted.

    Examples
    --------
    >>> from qiime.stats import fdr_correction
    >>> fdr_correction([.01, .2, .5, .1, .3])
    array([ 0.05      ,  0.33333333,  0.5       ,  0.25      ,  0.375     ])
    """
    tmp = array(pvals).astype(float)  # this converts Nones to nans
    return tmp * tmp.size / (1. + argsort(argsort(tmp)).astype(float))


def benjamini_hochberg_step_down(pvals):
    """Perform Benjamini and Hochberg's 1995 FDR step down procedure.

    Parameters
    ----------
    pvals : list or array
        List or array of floats.

    Returns
    -------
    list of pvals
        Returns the list of pvals properly adjusted based on the and then
        adjusted according to the BH rules.

    See Also
    --------
    fdr_correction

    Notes
    -----
    In short, computes the fdr adjusted pvals (ap_i's), and working from
    the largest to smallest, compare ap_i to ap_i-1. If ap_i < ap_i-1 set
    ap_i-1 equal to ap_i. Does *not* assume pvals is sorted. Described in [1]_.

    Examples
    --------
    >>> from qiime.stats import fdr_correction
    >>> benjamini_hochberg_step_down([0.1, 0.21, 0.5, 0.2, 0.6])
    array([ 0.35,  0.35,  0.6 ,  0.35,  0.6 ])

    References
    ----------
    .. [1] Controlling the False Discovery Rate: A Practical and Powerful
       Approach to Multiple Testing' Yoav Benjamini and Yosef Hochberg. Journal
       of the Royal Statistical Society. Series B (Methodological), Vol. 57,
       No. 1 (1995) 289-300.
    """
    tmp = fdr_correction(pvals)
    corrected_vals = empty(len(pvals))
    max_pval = 1.
    for i in argsort(pvals)[::-1]:
        if tmp[i] < max_pval:
            corrected_vals[i] = tmp[i]
            max_pval = tmp[i]
        else:
            corrected_vals[i] = max_pval
    return corrected_vals


def fisher_z_transform(r):
    """Calculate the Fisher Z transform of a correlation coefficient.

    Relies on formulation in [1_] pg 575.

    Parameters
    ----------
    r : float
        Correlation coefficient to transform.

    Returns
    -------
    z value of r

    References
    ----------
    .. [1] Sokal and Rohlf. "Biometry: The Principles and Practices of
       Statistics in Biological Research". ISBN: 978-0716724117
    """
    if abs(r) >= 1:  # fisher z transform is undefined, have to return nan
        return nan
    return .5 * log((1. + r) / (1. - r))


def inverse_fisher_z_transform(z):
    """Calculate the inverse of the Fisher Z transform on a z value.

    Relies on formulation in [1_] pg 576.

    Parameters
    ----------
    z : float
        z value of a correlation coefficient that has undergone transformation.

    Returns
    -------
    r : float
        Rho or correlation coefficient that would produce given z score.

    References
    ----------
    .. [1] Sokal and Rohlf. "Biometry: The Principles and Practices of
       Statistics in Biological Research". ISBN: 978-0716724117
    """
    return ((e ** (2 * z)) - 1.) / ((e ** (2 * z)) + 1.)


def z_transform_pval(z, n):
    '''Calculate two tailed probability of value as or more extreme than z.

    Relies on formulation in [1_] pg. 576.

    Parameters
    ----------
    z : float
        z-score
    n : int or float
        Number of samples that were used to generate the z-score.

    Returns
    -------
    zprob : float
        Probability of getting a zscore as or more extreme than the passed z
        given the total number of samples that generated it (n).

    References
    ----------
    .. [1] Sokal and Rohlf. "Biometry: The Principles and Practices of
       Statistics in Biological Research". ISBN: 978-0716724117
    '''
    if n <= 3:  # sample size must be greater than 3 otherwise this transform
        # isn't supported.
        return nan
    return normprob(z * ((n - 3) ** .5), direction='two-sided')


def normprob(z, direction='two-sided', mean=0, std=1):
    '''Calculate probability from normal distribution

    Paramaters
    ----------
    z : float
        Value of z statistic
    direction : str
        One of 'low', 'high', or 'two-sided'. Determines the bounds of the
        integration of the PDF. 'high' calculates the probability that a
        random variable Z will take a value as great or greater than z. 'low'
        will calculate the probability that Z will take a value less than or
        equal to z. 'two-sided' will calculate the probability that Z will take
        a value more extreme than z (i.e. abs(Z) >= z).
    mean : float
        Mean of the distirbution.
    std : float
        Standard deviation of the distribution.

    Returns
    -------
    p-value

    Notes
    -----
    scipy.stats.norm calculates the 'lower tail' of the distribution, i.e. the
    probability of a random variable Z taking a value smaller than or equal to
    the given z value.
    '''
    if direction == 'two-sided':
        if z >= 0:
            return 2 * (1. - norm.cdf(z, mean, std))
        else:
            return 2 * norm.cdf(z, mean, std)
    elif direction == 'high':
        return 1 - norm.cdf(z, mean, std)
    elif direction == 'low':
        return norm.cdf(z, mean, std)
    else:
        raise ValueError('Unknown direction.')


def chi2prob(x, df, direction='high'):
    '''Return the chi-squared statistic.

    Paramaters
    ----------
    x : float
        Value of x statistic.
    direction : str
        One of 'low' or 'high'. Determines the bounds of the
        integration of the PDF. 'high' calculates the probability that a
        random variable X will take a value as great or greater than x. 'low'
        will calculate the probability that X will take a value less than or
        equal to x.

    Returns
    -------
    p-value

    Notes
    -----
    scipy's chi2.cdf returns the 'lower tail' of the chi-squared distribution,
    that is p(X <= x). This necessitates adjustment of 1 - p for most qiime
    applications. However, with negative x a value of 0.0 is returned. Negative
    x are outside the domain of the CDF of chi-squared (and we should return a
    pval of nan in this case).
    '''
    if x <= 0:
        return nan
    elif direction == 'high':
        return 1. - chi2.cdf(x, df)
    elif direction == 'low':
        return chi2.cdf(x, df)
    else:
        raise ValueError('Unknown direction.')


def tprob(t, df, tails='high'):
    '''Calculate probability from t distribution

    Paramaters
    ----------
    t : float
        Value of t statistic
    tails : str
        One of 'low', 'high', or 'two-sided'. Determines the bounds of the
        integration of the PDF. 'high' calculates the probability that a
        random variable T will take a value as great or greater than t. 'low'
        will calculate the probability that T will take a value less than or
        equal to t. 'two-sided' will calculate the probability that T will take
        a value more extreme than t (i.e. abs(T) >= t).

    Returns
    -------
    p-value

    Notes
    -----
    scipy.stats.t calculates the 'lower tail' of the distribution, i.e. the
    probability of a random variable T taking a value smaller than or equal to
    the given t value.
    '''
    if tails == 'two-sided':
        if t >= 0:
            return 2 * (1. - tdist.cdf(t, df))
        else:
            return 2 * tdist.cdf(t, df)
    elif tails == 'high':
        return 1 - tdist.cdf(t, df)
    elif tails == 'low':
        return tdist.cdf(t, df)
    else:
        raise ValueError('Unknown direction.')


def fprob(f, dfn, dfd, direction='high'):
    '''Calculate probability from F distribution

    Paramaters
    ----------
    f : float
        Value of f statistic
    dfn : float
        Degrees of freedom for ???
    dfd : float
        Degrees of freedom for ???
    direction : str
        One of 'low' or 'high'. Determines the bounds of the
        integration of the PDF. 'high' calculates the probability that a
        random variable F will take a value as great or greater than f. 'low'
        will calculate the probability that F will take a value less than or
        equal to f.

    Returns
    -------
    p-value

    Notes
    -----
    scipy.stats.f calculates the 'lower tail' of the F distribution, ie the
    probability of a random variable F taking a value smaller than or equal to
    the given f value.
    '''
    if f < 0.:
        return nan
    elif direction == 'high':
        return 1. - fdist.cdf(f, dfn, dfd)
    elif direction == 'low':
        return fdist.cdf(f, dfn, dfd)
    else:
        raise ValueError('Unknown direction.')


def fisher_population_correlation(corrcoefs, sample_sizes):
    """Calculate population rho, homogeneity from corrcoefs using Z transform.

    Parameters
    ----------
    corrcoefs : array-like
        A list or array of floats.
    sample_sizes : array-like
        A list or array of ints.

    Returns
    -------
    rho : float
        Combined rho for the population.
    h_val : float
        probablity that the rhos of the different samples were homogenous.

    Notes
    -----
    This function calculates the correlation of a population that is relying
    on multiple different studies that have calculated correlation coefficients
    of their own. The procedure is detailed in [1]_ pgs
    576 - 578. Pvals that are nans will be excluded by this function.

    References
    ----------
    .. [1] Sokal and Rohlf. "Biometry: The Principles and Practices of
       Statistics in Biological Research". ISBN: 978-0716724117

    Examples
    --------
    >>> from qiime.stats import fisher_population_correlation
    >>> cc = [.4, .6, .7, .9]
    >>> ss = [10, 25, 14, 50]
    >>> fisher_population_correlation(cc, ss)
    (0.80559851846605035, 0.0029974748507499201)
    """
    tmp_rs = array(corrcoefs)
    tmp_ns = array(sample_sizes)
    # make checks for nans and exclude them as they will cause things to break
    rs = tmp_rs[~isnan(tmp_rs)]
    ns = tmp_ns[~isnan(tmp_rs)]
    if not (ns > 3).all():
        # not all samples have size > 3 which causes 0 varaince estimation.
        # thus we must return nan for pval and h_val
        return nan, nan
    if not len(ns) > 1:
        # only one sample, because of reduced degrees of freedom must have at
        # least two samples to calculate the homogeneity.
        return nan, nan
    if (rs >= 1.0).any():
        # a failure will occur in chi_high calculation where an non-terminating
        # loop will be initiated.
        raise ValueError('A correlation coefficient >= 1 was passed. This is '
                         'a non real valured correlation coefficient and it\'s'
                         ' Fisher Z transform cannot be computed.')
    # calculate zs
    zs = array([fisher_z_transform(float(i)) for i in rs])
    # calculate variance weighted z average = z_bar
    z_bar = (zs * (ns - 3)).sum() / float((ns - 3).sum())
    rho = inverse_fisher_z_transform(z_bar)
    # calculate homogeneity
    x_2 = ((ns - 3) * (zs - z_bar) ** 2).sum()
    h_val = chi2prob(x_2, len(ns) - 1, direction='high')
    return rho, h_val


def cscore(v1, v2):
    '''Calculate C-score between v1 and v2 according to Stone and Roberts 1990.

    Parameters
    ----------
    v1 : array-like
        List or array of numeric values to be tested.
    v2 : array-like
        List or array of numeric values to be tested.

    Returns
    -------
    cscore : float
        C-score between v1 and v2

    Notes
    -----
    This function calculates the C-score between equal length vectors v1 and v2
    according to the formulation given in [1]_.

    References
    ----------
    .. [1] Stone and Roberts. 1990, Oecologia 85:74-79
    '''
    v1_b = v1.astype(bool)
    v2_b = v2.astype(bool)
    sij = (v1_b * v2_b).sum()
    return (v1_b.sum() - sij) * (v2_b.sum() - sij)

def correlate(v1, v2, method):
    '''Correlate vectors using method.

    Parameters
    ----------
    v1 : array-like
        List or array of ints or floats to be correlated.
    v2 : array-like
        List or array of ints or floats to be correlated.
    method : str
        One of 'spearman', 'pearson', 'kendall', 'cscore'.

    Returns
    -------
    rho : float
        Correlation between the vectors.
    '''
    if method == 'pearson':
        corr_fn = pearson
    elif method == 'spearman':
        corr_fn = spearman
    elif method == 'kendall':
        corr_fn = kendall
    elif method == 'cscore':
        corr_fn = cscore
    else:
        raise ValueError('Correlation function not recognized.')
    return corr_fn(v1, v2)
