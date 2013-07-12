#!/usr/bin/env python
from __future__ import division

__author__ = "Michael Dwan"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Michael Dwan", "Logan Knecht",
               "Damien Coy", "Levi McCracken", "Andrew Cochran",
               "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"
__status__ = "Development"

"""
This module provides functionality for the application of various statistical
methods to QIIME-formatted datasets.

The module provides an API that allows users to easily apply any number of
statistical analyses and just as easily retrieve the results. The module also
provides a hierarchy of statistical classes that can be inherited from to
create new statistical method implementations.
"""

from os.path import join
from types import ListType
from copy import deepcopy
from matplotlib import use
use('Agg', warn=False)
from matplotlib.pyplot import figure
from numpy import (argsort, array, asarray, ceil, empty, fill_diagonal, finfo,
        log2, mean, ones, sqrt, tri, unique, zeros, ndarray, floor, median)
from numpy import argsort, min as np_min, max as np_max
from numpy.random import permutation
from cogent.util.misc import combinate, create_dir
from cogent.maths.stats.test import t_one_sample

from qiime.pycogent_backports.test import (mantel_test, mc_t_two_sample,
                                           pearson, permute_2d, spearman)
from qiime.format import format_p_value_for_num_iters
from qiime.util import DistanceMatrix, MetadataMap

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
    qiime.make_distance_histograms.monte_carlo_group_distances code.

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
        stat = ['N/A' if e is None else e for e in stat]
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

    # Convert our notion of tail type into the format expected by
    # PyCogent.
    if tail_type == 'two-sided':
        tail_type = None

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
                obs_t, param_p_val, nonparam_p_val = None, None, None
            else:
                obs_t, param_p_val, _, nonparam_p_val = mc_t_two_sample(
                        g1_dist, g2_dist, tails=tail_type,
                        permutations=num_permutations)
            result.append([g1_label, g2_label, obs_t, param_p_val, None,
                           nonparam_p_val, None])
            if obs_t is not None:
                num_tests += 1

    # Correct the p-values for multiple comparisons, now that we know how many
    # tests succeeded.
    for stat in result:
        stat[4] = stat[3] if stat[3] is None else min(stat[3] * num_tests, 1)
        stat[6] = stat[5] if stat[5] is None else min(stat[5] * num_tests, 1)
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

    assert type(data) == list or type(data) == ndarray, "Data must be either"+\
        " a Python list or a NumPy 1-D array"
    assert type(quantiles) == list or type(quantiles) == ndarray, "Quantiles"+\
        " must be either a Python list or a NumPy 1-D array"
    assert all(map(lambda x: x>=0 and x<=1, quantiles)), "All the elements "+\
        "in the quantiles list must be greater than 0 and lower than one"

    # unless the user wanted, do not modify the data
    data = deepcopy(data)

    if type(data) != ndarray:
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
    index = quantile*(len(data)-1)
    bottom_index = int(floor(index))
    top_index = int(ceil(index))

    difference = index-bottom_index
    output = (1-difference)*data[bottom_index]+difference*data[top_index]

    return output

class DistanceMatrixStats(object):
    """Base class for distance matrix-based statistical methods.

    This class provides an interface to setting and accessing an arbitrary
    number of distance matrices. Users of this class can optionally specify the
    number of allowable distance matrices and their minimum allowable size (the
    default is no restrictions on either of these).

    It is the parent class of CorrelationStats and CategoryStats.
    """

    def __init__(self, dms, num_dms=-1, min_dm_size=-1,
                 suppress_symmetry_and_hollowness_check=False):
        """Default constructor.

        Initializes an instance with the provided list of distance matrices.

        Arguments:
            dms - a list of DistanceMatrix objects
            num_dms - the exact number of allowable distance matrices. If -1
                (the default), there is no restriction on how many distance
                matrices the user can set
            min_dm_size - the minimum size that all distance matrices must have
                that are stored by this instance. If -1, no size restriction
            suppress_symmetry_and_hollowness_check - by default, the input
                distance matrices will be checked for symmetry and hollowness.
                It is recommended to leave this check in place for safety, as
                the check is fairly fast. However, if you *know* you have
                symmetric and hollow distance matrices, you can disable this
                check for small performance gains on extremely large distance
                matrices. Alternatively, if the statistical method works on
                asymmetric and/or non-hollow distance matrices, you can disable
                this check to allow for these types of distance matrices
        """
        self._num_dms = num_dms
        self._min_dm_size = min_dm_size
        self._suppress_symmetry_and_hollowness_check = \
                suppress_symmetry_and_hollowness_check
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
                raise TypeError('Invalid type: %s; expected DistanceMatrix' %
                                dm.__class__.__name__)
            if self._min_dm_size >= 0 and dm.Size < self._min_dm_size:
                raise ValueError("Distance matrix of size %dx%d is smaller "
                                 "than the minimum allowable distance matrix "
                                 "size of %dx%d for this analysis." %
                                 (dm.Size, dm.Size, self._min_dm_size,
                                  self._min_dm_size))
            if not self._suppress_symmetry_and_hollowness_check:
                if not dm.is_symmetric_and_hollow():
                    raise ValueError("The distance matrix must be symmetric "
                                     "and hollow.")
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
    that compare two or more distance matrices.

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

    def __init__(self, dms, num_dms=-1, min_dm_size=-1,
                 suppress_symmetry_and_hollowness_check=False):
        """Default constructor.

        Creates a new instance with the provided list of distance matrices.

        Arguments:
            dms - a list of DistanceMatrix objects
            num_dms - the exact number of allowable distance matrices. If -1
                (the default), there is no restriction on how many distance
                matrices the user can set
            min_dm_size - the minimum size that all distance matrices must have
                that are stored by this instance. If -1, no size restriction
            suppress_symmetry_and_hollowness_check - by default, the input
                distance matrices will be checked for symmetry and hollowness.
                It is recommended to leave this check in place for safety, as
                the check is fairly fast. However, if you *know* you have
                symmetric and hollow distance matrices, you can disable this
                check for small performance gains on extremely large distance
                matrices. Alternatively, if the statistical method works on
                asymmetric and/or non-hollow distance matrices, you can disable
                this check to allow for these types of distance matrices
        """
        super(CorrelationStats, self).__init__(dms, num_dms, min_dm_size,
                suppress_symmetry_and_hollowness_check)

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
        size = dms[0].Size
        sample_ids = dms[0].SampleIds
        for dm in dms:
            if dm.Size != size:
                raise ValueError("All distance matrices must have the same "
                                 "number of rows and columns.")
            if dm.SampleIds != sample_ids:
                raise ValueError("All distance matrices must have matching "
                                 "sample IDs.")


class CategoryStats(DistanceMatrixStats):
    """Base class for categorical statistical analyses.

    It is subclassed by categorical statistical methods such as DB-RDA or BEST.
    Categorical statistical methods usually have some categorical grouping of
    samples, and the significance of this grouping is usually what is tested.
    For example, are treatment samples significantly different from control
    samples? This is not always the case (e.g. DB-RDA is an ordination
    technique), but most of the categorical methods follow this general design.

    A valid instance of CategoryStats must have at least one distance matrix
    and a single metadata map containing the sample IDs of the distance matrix
    or matrices.
    """

    def __init__(self, mdmap, dms, cats, num_dms=-1, min_dm_size=-1,
                 random_fn=permutation,
                 suppress_symmetry_and_hollowness_check=False,
                 suppress_category_uniqueness_check=False,
                 suppress_numeric_category_check=True,
                 suppress_single_category_value_check=False):
        """Default constructor.

        Creates a new instance with the provided distance matrices,
        metadata map, and list of categories.

        Arguments:
            mdmap - a MetadataMap instance
            dms - a list of DistanceMatrix objects
            cats - a list of strings denoting categories in the metadata map
                that will be used by this analysis (i.e. the grouping
                variable(s))
            num_dms - the exact number of allowable distance matrices. If -1
                (the default), there is no restriction on how many distance
                matrices the user can set
            min_dm_size - the minimum size that all distance matrices must have
                that are stored by this instance. If -1, no size restriction
            random_fn - the function to use when randomizing the grouping
                of samples in a category during calculation of the p-value. It
                must return a value and must be callable
            suppress_symmetry_and_hollowness_check - by default, the input
                distance matrices will be checked for symmetry and hollowness.
                It is recommended to leave this check in place for safety, as
                the check is fairly fast. However, if you *know* you have
                symmetric and hollow distance matrices, you can disable this
                check for small performance gains on extremely large distance
                matrices. Alternatively, if the statistical method works on
                asymmetric and/or non-hollow distance matrices, you can disable
                this check to allow for these types of distance matrices
            suppress_category_uniqueness_check - by default, each input
                category will be checked to ensure that not all values are
                unique (i.e. some duplicated values exist). In other words,
                this check makes sure that the category will group samples such
                that there exists at least one group of samples that is
                composed of two or more samples. Many categorical statistical
                methods (such as ANOSIM and PERMANOVA) need this requirement to
                be met to avoid erroneous math (e.g. division by zero) due to a
                lack of 'within' distances. An example of a unique category
                would be SampleID, where each category value is unique
            suppress_numeric_category_check - if False, each category's values
                will be checked to ensure they can be converted to a float.
                Useful for methods that only accept numeric categories
            suppress_single_category_value_check - if False, each category's
                values will be checked to ensure they are not all the same.
                Many of the methods will not work if every value is the same
                (i.e. there is only one single group of samples)
        """
        super(CategoryStats, self).__init__(dms, num_dms, min_dm_size,
                suppress_symmetry_and_hollowness_check)
        self._suppress_category_uniqueness_check = \
                suppress_category_uniqueness_check
        self._suppress_numeric_category_check = suppress_numeric_category_check
        self._suppress_single_category_value_check = \
                suppress_single_category_value_check
        self.MetadataMap = mdmap
        self.Categories = cats
        self.RandomFunction = random_fn
        self._validate_compatibility()

    @property
    def MetadataMap(self):
        """Returns the instance's metadata map.

        The metadata map is returned as a MetadataMap class instance.
        """
        return self._metadata_map

    @MetadataMap.setter
    def MetadataMap(self, new_mdmap):
        """Sets the instance's metadata map.

        Arguments:
            new_mdmap - A MetadataMap object instance
        """
        if not isinstance(new_mdmap, MetadataMap):
            raise TypeError('Invalid type: %s; not MetadataMap' %
                            new_mdmap.__class__.__name__)
        self._metadata_map = new_mdmap

    @property
    def Categories(self):
        """Gets the instance's categories.

        Returns a list of mapping file category name strings.
        """
        return self._categories

    @Categories.setter
    def Categories(self, new_categories):
        """Sets the instance's list of categories.

        Arguments:
            new_categories - A list of category name strings. These must be
                present in the current metadata map
        """
        if not isinstance(new_categories, ListType):
            raise TypeError("The supplied categories must be a list of "
                            "strings.")
        for new_cat in new_categories:
            if not isinstance(new_cat, str):
                raise TypeError("Invalid category: not of type 'string'")
            elif new_cat not in self._metadata_map.CategoryNames:
                raise ValueError("The category '%s' is not in the mapping "
                                 "file." % new_cat)

            if not self._suppress_numeric_category_check:
                if not self._metadata_map.isNumericCategory(new_cat):
                    raise TypeError("The category '%s' is not numeric. Not "
                                    "all values could be converted to numbers."
                                    % new_cat)

            if not self._suppress_category_uniqueness_check:
                if self._metadata_map.hasUniqueCategoryValues(new_cat):
                    raise ValueError("All values in category '%s' are unique. "
                                     "This statistical method cannot operate "
                                     "on a category with unique values (e.g. "
                                     "there are no 'within' distances because "
                                     "each group of samples contains only a "
                                     "single sample)." % new_cat)

            if not self._suppress_single_category_value_check:
                if self._metadata_map.hasSingleCategoryValue(new_cat):
                    raise ValueError("All values in category '%s' are the "
                                     "same. This statistical method cannot "
                                     "operate on a category that creates only "
                                     "a single group of samples (e.g. there "
                                     "are no 'between' distances because "
                                     "there is only a single group)."
                                     % new_cat)

        self._categories = new_categories

    @property
    def RandomFunction(self):
        """Returns the randomization function used in p-value calculations."""
        return self._random_fn

    @RandomFunction.setter
    def RandomFunction(self, random_fn):
        """Setter for the randomization function used in p-value calcs.

        Arguments:
            random_fn - the function to use when randomizing the grouping
                during calculation of the p-value. It must return a value and
                must be callable
        """
        if hasattr(random_fn, '__call__'):
            self._random_fn = random_fn
        else:
            raise TypeError("The supplied function reference is not callable.")

    def _validate_compatibility(self):
        """Checks that the current dms, map, and categories are compatible.

        This method will raise an error if any of the sample IDs in any of the
        distance matrices are not found in the metadata map. Ordering of
        sample IDs is not taken into account. An error will also be raised if
        any categories cannot be found in the mapping file.

        This method exists because we do not have a method to set distance
        matrices, metadata map, and categories at the same time.
        """
        for dm in self.DistanceMatrices:
            for samp_id in dm.SampleIds:
                if samp_id not in self.MetadataMap.SampleIds:
                    raise ValueError("The sample ID '%s' was not found in the "
                                     "metadata map." % samp_id)
        for cat in self.Categories:
            if cat not in self.MetadataMap.CategoryNames:
                raise ValueError("The category '%s' was not found in the "
                                 "metadata map." % cat)

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
        # Make sure the *current* distance matrices and metadata map are
        # compatible before continuing.
        self._validate_compatibility()
        return super(CategoryStats, self).__call__(num_perms)


class Anosim(CategoryStats):
    """Class for the ANOSIM categorical statistical analysis.

    Briefly, ANOSIM tests whether two or more groups of samples are
    significantly different. The user of the class specifies a category in the
    metadata map to group samples by.

    This code is heavily based on Andrew Cochran's original procedural version.
    """

    def __init__(self, mdmap, dm, cat, random_fn=permutation,
                 suppress_symmetry_and_hollowness_check=False):
        """Initializes an instance with the specified analysis parameters.

        WARNING: Only symmetric, hollow distance matrices may be used as input.
        Asymmetric distance matrices, such as those obtained by the UniFrac
        Gain metric (i.e. beta_diversity.py -m unifrac_g), should not be used
        as input.

        Arguments:
            mdmap - the MetadataMap instance to obtain grouping info from
            dm - the DistanceMatrix instance to obtain distances from
            cat - the category string to group samples by (must be in the
                metadata map)
            random_fn - the function to use when randomizing the grouping
                during calculation of the p-value. It must return a value and
                must be callable
            suppress_symmetry_and_hollowness_check - by default, the input
                distance matrix will be checked for symmetry and hollowness.
                It is recommended to leave this check in place for safety, as
                the check is fairly fast. However, if you *know* you have
                a symmetric and hollow distance matrix, you can disable this
                check for small performance gains on extremely large distance
                matrices
        """
        super(Anosim, self).__init__(mdmap, [dm], [cat], num_dms=1,
                random_fn=random_fn, suppress_symmetry_and_hollowness_check=\
                suppress_symmetry_and_hollowness_check)

    def __call__(self, num_perms=999):
        """Runs ANOSIM on the current distance matrix and sample grouping.

        Returns a dict containing the results. The following keys are set:
            method_name - name of the statistical method
            r_value - the ANOSIM R statistic computed by the test
            p_value - the p-value computed by the test, or 'NA' if the number
                of permutations was zero
            num_perms - the number of permutations used when calculating the
                p-value

        Arguments:
            num_perms - the number of permutations to use when calculating the
                p-value
        """
        results = super(Anosim, self).__call__(num_perms)
        category = self.Categories[0]
        samples = self.DistanceMatrices[0].SampleIds

        # Create the group map, which maps sample ID to category value (e.g.
        # sample 1 to 'control' and sample 2 to 'fast').
        group_map = {}
        for samp_id in samples:
            group_map[samp_id] = self.MetadataMap.getCategoryValue(
                    samp_id, category)

        # Calculate the R statistic with the grouping found in the current
        # metadata map.
        r_stat = self._anosim(group_map)

        if num_perms > 0:
            # Calculate the p-value based on the number of permutations.
            perm_stats = []
            for i in range(num_perms):
                # Randomize grouping. We don't use values() in order to
                # preserve ordering in case the user's random function doesn't
                # change the order of the items in the list.
                grouping_random = [group_map[sample] for sample in samples]
                grouping_random = self.RandomFunction(grouping_random)
                for j, sample in enumerate(samples):
                    group_map[sample] = grouping_random[j]
                perm_stats.append(self._anosim(group_map))
            # Calculate the p-value.
            p_value = (sum(perm_stats >= r_stat) + 1) / (num_perms + 1)
        else:
            p_value = 1.0

        results['method_name'] = 'ANOSIM'
        results['r_value'] = r_stat
        results['p_value'] = p_value
        results['num_perms'] = num_perms

        return results

    def _anosim(self, group_map):
        """Computes ANOSIM on the supplied grouping, returning the R value.

        The R value is between -1 and 1 and indicates the strength of the
        grouping.

        Arguments:
            group_map - a python dict mapping sample ID to category value (e.g.
                sample 1 to 'control' and sample 2 to 'fast'). This map must
                contain a key for each sample ID in the current distance
                matrix
        """
        dm = self.DistanceMatrices[0]
        dm_size = dm.Size

        # Create grouping matrix, where a one means that the two samples are in
        # the same group (e.g. control) and a zero means that they aren't.
        within_between = zeros((dm_size, dm_size))
        for i, i_sample in enumerate(dm.SampleIds):
            for j, j_sample in enumerate(dm.SampleIds):
                if group_map[i_sample] == group_map[j_sample]:
                    within_between[i][j] = 1

        # Extract upper triangle from the distance and grouping matrices.
        distances = dm.DataMatrix[tri(dm_size) == 0]
        grouping = within_between[tri(dm_size) == 0]

        # Sort extracted data.
        sorted_distances = []
        sorted_grouping = []
        for idx in argsort(distances):
             sorted_distances.append(distances[idx])
             sorted_grouping.append(grouping[idx])

        # Account for rank ties, then compute R statistic.
        rank_list = range(1, len(sorted_distances) + 1)
        adjusted_rank_list = self._remove_ties(sorted_distances, rank_list)
        return self._compute_r_value(adjusted_rank_list, sorted_grouping,
                                     dm_size)

    def _remove_ties(self, sorted_dists, ranks):
        """Replaces repeat values with the average of them.

        Returns a list containing the adjusted ranks.

        Arguments:
            sorted_dists: list of the sorted distances
            ranks: list containing the ranks of each of the differences
        """
        result = []
        ties = []
        tie_count = 0
        tie_flag = 0

        for i in range(len(sorted_dists) - 1):
            # Store state information.
            curr_dist = sorted_dists[i]
            next_dist = sorted_dists[i+1]
            rank_val = ranks[i]

            # A tie has not occured yet.
            if tie_flag == 0:
                if curr_dist == next_dist:
                    # We have a tie, so add the current rank to the tie list.
                    tie_count = tie_count + 1
                    ties.append(rank_val)
                    first_tie_index = i
                    tie_flag = 1
                else:
                    # If no tie, fill in the list with the current rank.
                    result.append(rank_val)
            else:
                # A tie has already occured.
                if curr_dist == next_dist:
                    # If another tie occurs, add the current rank to the tie
                    # list.
                    tie_count = tie_count + 1
                    ties.append(rank_val)
                else:
                    # No more ties, average their values and attach to adjusted
                    # rank list.
                    ties.append(rank_val)
                    last_tie_index = i
                    result.extend(self._get_adjusted_vals(ties,
                            first_tie_index, last_tie_index))
                    tie_flag = 0
                    tie_count = 0
                    ties = []
        # If there is a tie that extends to the final position, we must process
        # it here to avoid out of list bounds errors.
        if tie_flag == 1:
            ties.append(ranks[i+1])
            last_tie_index = i + 1
            result.extend(self._get_adjusted_vals(ties, first_tie_index,
                                                  last_tie_index))
        else:
            result.append(ranks[i+1])
        return result

    def _get_adjusted_vals(self, ties, first_tie_idx, last_tie_idx):
        """Helper function to _remove_ties. Consolidates repeated code."""
        adjusted_val = sum(ties) / len(ties)
        return [adjusted_val] * ((last_tie_idx - first_tie_idx) + 1)

    def _compute_r_value(self, adjusted_ranks, sorted_groups, num_samps):
        """Code that performs the actual math involved in solving ANOSIM.

        Returns the ANOSIM R value (between -1 and 1).

        Arguments:
            adjusted_ranks - list of the ranks, adjusted for ties
            sorted_groups - list associating distances to groups
            num_samps: how many total samples
        """
        adjusted_ranks = array(adjusted_ranks)
        sorted_groups = array(sorted_groups)

        # Compute r_W and r_B.
        r_W = mean(adjusted_ranks[sorted_groups==1])
        r_B = mean(adjusted_ranks[sorted_groups==0])
        divisor = num_samps * ((num_samps - 1) / 4)
        return (r_B - r_W) / divisor


class Permanova(CategoryStats):
    """Class for the PERMANOVA statistical method.

    This is a non-parametric, permutation-based method to determine the
    significance of sample grouping.

    This code is heavily based on Andrew Cochran's original procedural version.
    """

    def __init__(self, mdmap, dm, cat, random_fn=permutation,
                 suppress_symmetry_and_hollowness_check=False):
        """Initializes an instance with the specified analysis parameters.

        WARNING: Only symmetric, hollow distance matrices may be used as input.
        Asymmetric distance matrices, such as those obtained by the UniFrac
        Gain metric (i.e. beta_diversity.py -m unifrac_g), should not be used
        as input.

        Arguments:
            mdmap - the MetadataMap instance to obtain grouping info from
            dm - the DistanceMatrix instance to obtain distances from
            cat - the category string to group samples by (must be in the
                metadata map)
            num_perms - the number of permutations to use when calculating the
                p-value. If zero, the p-value will not be calculated. Must be
                greater than or equal to zero
            random_fn - the function to use when randomizing the grouping
                during calculation of the p-value. It must return a value and
                must be callable
            suppress_symmetry_and_hollowness_check - by default, the input
                distance matrix will be checked for symmetry and hollowness.
                It is recommended to leave this check in place for safety, as
                the check is fairly fast. However, if you *know* you have
                a symmetric and hollow distance matrix, you can disable this
                check for small performance gains on extremely large distance
                matrices
        """
        super(Permanova, self).__init__(mdmap, [dm], [cat], num_dms=1,
                random_fn=random_fn, suppress_symmetry_and_hollowness_check=\
                suppress_symmetry_and_hollowness_check)

    def __call__(self, num_perms=999):
        """Runs PERMANOVA on the current distance matrix and sample grouping.

        Returns a dict containing the results. The following keys are set:
            method_name - name of the statistical method
            f_value - the PERMANOVA F statistic computed by the test
            p_value - the p-value computed by the test, or 'NA' if the number
                of permutations was zero
            num_perms - the number of permutations used when calculating the
                p-value

        Arguments:
            num_perms - the number of permutations to use when calculating the
                p-value
        """
        results = super(Permanova, self).__call__(num_perms)
        category = self.Categories[0]
        samples = self.DistanceMatrices[0].SampleIds

        # Create the group map, which maps sample ID to category value (e.g.
        # sample 1 to 'control' and sample 2 to 'fast').
        group_map = {}
        for samp_id in samples:
            group_map[samp_id] = self.MetadataMap.getCategoryValue(
                    samp_id, category)

        # Calculate the F statistic with the grouping found in the current
        # metadata map.
        f_stat = self._permanova(group_map)

        if num_perms > 0:
            # Calculate the p-value based on the number of permutations.
            perm_stats = []
            for i in range(num_perms):
                # Randomize grouping. We don't use values() in order to
                # preserve ordering in case the user's random function doesn't
                # change the order of the items in the list.
                grouping_random = [group_map[sample] for sample in samples]
                grouping_random = self.RandomFunction(grouping_random)
                for j, sample in enumerate(samples):
                    group_map[sample] = grouping_random[j]
                perm_stats.append(self._permanova(group_map))
            # Calculate the p-value.
            p_value = (sum(perm_stats >= f_stat) + 1) / (num_perms + 1)
        else:
            p_value = 1.0

        results['method_name'] = 'PERMANOVA'
        results['f_value'] = f_stat
        results['p_value'] = p_value
        results['num_perms'] = num_perms

        return results

    def _permanova(self, grouping):
        """Computes PERMANOVA pseudo-F-statistic.

        Arguments:
            grouping - a python dict mapping sample ID to category value (e.g.
                sample 1 to 'control' and sample 2 to 'fast'). This map must
                contain a key for each sample ID in the current distance
                matrix
        """
        samples = self.DistanceMatrices[0].SampleIds
        dm = self.DistanceMatrices[0]
        # Number of samples in each group.
        unique_n = []
        group_map = {}

        # Extract the unique list of group labels.
        gl_unique = unique(array(grouping.values()))

        # Calculate number of groups and unique 'n's.
        number_groups = len(gl_unique)
        for i, i_string in enumerate(gl_unique):
            group_map[i_string] = i
            unique_n.append(grouping.values().count(i_string))

        # Create grouping matrix.
        grouping_matrix = -1 * ones((dm.Size, dm.Size))
        for i, i_sample in enumerate(samples):
            grouping_i = grouping[i_sample]
            for j, j_sample in enumerate(samples):
                if grouping_i == grouping[j_sample]:
                    grouping_matrix[i][j] = group_map[grouping[i_sample]]

        # Extract upper triangle.
        distances = dm[tri(dm.Size) == 0]
        groups = grouping_matrix[tri(len(grouping_matrix)) == 0]

        # Compute F value.
        return self._compute_f_value(distances, groups, dm.Size,
                                     number_groups, unique_n)

    def _compute_f_value(self, distances, groupings, number_samples,
                         number_groups, unique_n):
        """Performs the calculations for the F value.

        Arguments:
            distances - a list of the distances
            groupings - a list associating the distances to their groups
            number_samples - how many samples there are
            number_groups - how many groups there are
            unique_n - list containing how many samples are in each within
                group
        """
        a = number_groups
        N = number_samples

        # Calculate s_T.
        s_T = sum(distances * distances) / N

        # Calculate s_W for each group, this accounts for different group
        # sizes.
        s_W = 0
        for i in range(number_groups):
            group_ix = groupings==i
            diffs = distances[group_ix]
            s_W = s_W + sum(diffs**2) / unique_n[i]

        # Execute the formula.
        s_A = s_T - s_W
        return (s_A / (a-1)) / (s_W / (N-a))


class Best(CategoryStats):
    """Class for the BEST/BioEnv statistical analysis.
    
    Based on vegan::bioenv function, which is an implementation of the BEST
    statistical method.
    """

    def __init__(self, dm, metadata_map, cats,
                 suppress_symmetry_and_hollowness_check=False):
        """Default constructor.

        WARNING: Only symmetric, hollow distance matrices may be used as input.
        Asymmetric distance matrices, such as those obtained by the UniFrac
        Gain metric (i.e. beta_diversity.py -m unifrac_g), should not be used
        as input.

        Arguments:
            dm - the DistanceMatrix instance to run the analysis on
            metadata_map - the MetadataMap instance to obtain category
                information from
            cats - list of category strings in the metadata map that will be
                used to determine the combination that "best" explains the
                variability in the data. Each category must be numeric
            suppress_symmetry_and_hollowness_check - by default, the input
                distance matrix will be checked for symmetry and hollowness.
                It is recommended to leave this check in place for safety, as
                the check is fairly fast. However, if you *know* you have
                a symmetric and hollow distance matrix, you can disable this
                check for small performance gains on extremely large distance
                matrices
        """
        # BEST doesn't require non-unique categories or non-single value
        # categories, but *does* require only numeric categories.
        super(Best, self).__init__(metadata_map, [dm], cats, num_dms=1,
              suppress_symmetry_and_hollowness_check=\
              suppress_symmetry_and_hollowness_check,
              suppress_category_uniqueness_check=True,
              suppress_numeric_category_check=False,
              suppress_single_category_value_check=True)

    def __call__(self, num_perms=999):
        """Runs the BEST/BioEnv analysis on a distance matrix using specified
        metadata map categories.

        num_perms is ignored, but maintained for a consistent interface with
        the other statistical method classes.

        Returns a dictionary which contains the resulting data. Keys:
            method_name - name of the statistical method
            num_vars - number of categories (variables)
            vars - mapping of category names to indices
            rho_vals - spearman correlation statistics, one for each
                combination of vars
        """
        res = super(Best, self).__call__()
        cats = self.Categories
        dm = self.DistanceMatrices[0]
        dm_flat = dm.flatten()

        row_count = dm.Size
        col_count = len(cats)
        sum = 0
        stats = [(-777777777, '') for c in range(col_count+1)]
        for i in range(1, col_count+1):
            combo = list(combinate([j for j in range(1,col_count+1)], i))

            for c in range(len(combo)):
                cat_mat = self._make_cat_mat(cats, combo[c])
                cat_dm = self._derive_euclidean_dm(cat_mat, row_count)
                cat_dm_flat = cat_dm.flatten()
                r = spearman(dm_flat, cat_dm_flat)
                if r > stats[i-1][0]:
                    stats[i-1] = (r, ','.join(str(s) for s in combo[c]))

        res['method_name'] = 'BEST'
        res['num_vars'] = col_count
        res['vars'] = ['%s = %d' % (name,val+1)
                       for val,name in enumerate(cats)]
        res['rho_vals'] = stats[:-1]

        return res

    def _derive_euclidean_dm(self, cat_mat, dim):
        """Returns an n x n, euclidean distance matrix, where n = len(cats)."""
        dm_labels = self.DistanceMatrices[0].SampleIds
        res_mat = []
        for i in range(dim):
            res_mat.append([0 for k in range(dim)])
            for j in range(i):
                res_mat[i][j] = self._vector_dist(cat_mat[i], cat_mat[j])
                res_mat[j][i] = res_mat[i][j]

        return DistanceMatrix(asarray(res_mat), dm_labels, dm_labels)

    def _vector_dist(self, vec1, vec2):
        """Calculates the Euclidean distance between two vectors."""
        return sqrt(sum([(float(v1) - float(v2))**2 for v1,v2 in
                            zip(vec1,vec2)]))

    def _make_cat_mat(self, cats, combo):
        """Returns a matrix with columns pulled from category values.

        Returns a matrix with len(sample_ids) rows of columns pulled from
        category values, the number of columns for each category is
        determined by the current combination (combo).
        """
        dm = self.DistanceMatrices[0]
        md_map = self.MetadataMap
        res = []
        for i in combo:
            res.append(md_map.getCategoryValues(dm.SampleIds, cats[i-1]))
        return zip(*res)


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
                 suppress_symmetry_and_hollowness_check=False,
                 variable_size_distance_classes=False):
        """Constructs a new MantelCorrelogram instance.

        WARNING: Only symmetric, hollow distance matrices may be used as input.
        Asymmetric distance matrices, such as those obtained by the UniFrac
        Gain metric (i.e. beta_diversity.py -m unifrac_g), should not be used
        as input.

        Arguments:
            eco_dm - a DistanceMatrix object representing the ecological
                distances between samples (e.g. UniFrac distance matrix)
            geo_dm - a DistanceMatrix object representing some other distance
                measure between samples (most commonly geographical distances,
                but could also be distances in pH, temperature, etc.)
            alpha - the alpha value to use when marking the Mantel
                correlogram plot for significance
            suppress_symmetry_and_hollowness_check - by default, the input
                distance matrices will be checked for symmetry and hollowness.
                It is recommended to leave this check in place for safety, as
                the check is fairly fast. However, if you *know* you have
                symmetric and hollow distance matrices, you can disable this
                check for small performance gains on extremely large distance
                matrices
            variable_size_distance_classes - if True, distance classes (bins)
                will vary in size such that each distance class (bin) will have
                the same number of distances. If False, all distance classes
                will have the same size, though the number of distances in each
                class may not be equal. Having variable-sized distance classes
                can help maintain statistical power if there are large
                differences in the number of distances in each class
        """
        super(MantelCorrelogram, self).__init__([eco_dm, geo_dm], num_dms=2,
                min_dm_size=3, suppress_symmetry_and_hollowness_check=\
                suppress_symmetry_and_hollowness_check)
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
        dm_size = eco_dm.Size

        # Find the number of lower triangular elements (excluding the
        # diagonal).
        num_dists = dm_size * (dm_size - 1) // 2

        # Use Sturge's rule to determine the number of distance classes.
        num_classes = int(ceil(1 + log2(num_dists)))

        # Create the matrix of distance classes. Each element in the matrix
        # contains what distance class the original element is in. Also find
        # the distance class indices, which are the midpoints in each distance
        # class.
        dist_class_matrix, class_indices = self._find_distance_classes(geo_dm,
            num_classes)

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
            model_matrix = DistanceMatrix(model_matrix, geo_dm.SampleIds,
                                          geo_dm.SampleIds)

            # Count the number of distances in the current distance class.
            num_distances = int(model_matrix.sum())
            results['num_dist'].append(num_distances)
            if num_distances == 0:
                results['mantel_r'].append(None)
                results['mantel_p'].append(None)
            else:
                row_sums = model_matrix.sum(axis='observation')
                row_sums = map(int, row_sums)
                has_zero_sum = 0 in row_sums

                # Only stop running Mantel tests if we've gone through half of
                # the distance classes and at least one row has a sum of zero
                # (i.e. the sample doesn't have any distances that fall in the
                # current class).
                if not (class_num > ((num_classes // 2) - 1) and has_zero_sum):
                    mantel_test = Mantel(model_matrix, eco_dm,
                                         tail_type='greater')
                    mantel_test_results = mantel_test(num_perms)
                    p_val, orig_stat, perm_stats = (
                            mantel_test_results['p_value'],
                            mantel_test_results['r_value'],
                            mantel_test_results['perm_stats'])

                    # Negate the Mantel r statistic because we are using
                    # distance matrices, not similarity matrices (this is a
                    # necessary step, see Legendre's Numerical Ecology
                    # algorithm reference for more details).
                    results['mantel_r'].append(-orig_stat)

                    # The mantel function produces a one-tailed p-value
                    # (H1: r>0). Here, compute a one-tailed p-value in the
                    # direction of the sign.
                    if orig_stat < 0:
                        perm_sum = sum([1 for ps in perm_stats \
                            if ps <= orig_stat]) + 1
                        p_val = perm_sum / (num_perms + 1)
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

        dm_lower_flat = dm.flatten()
        size = dm.Size

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
        DistanceMatrix.flatten(lower=True).
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
        are None, they are not counted towards the number of tests used in the
        correction.

        Returns a list of Bonferroni-corrected p-values for those that are not
        None. Those that are None are simply returned. The ordering of p-values
        is maintained.

        Arguments:
            p_vals - list of p-values (of type float or None)
        """
        num_tests = len([p_val for p_val in p_vals if p_val is not None])
        corrected_p_vals = []
        for p_val in p_vals:
            if p_val is not None:
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
            if p_val <= self.Alpha:
                signif_classes.append(class_indices[idx])
                signif_stats.append(mantel_stats[idx])
        ax.plot(signif_classes, signif_stats, 'ks', mfc='k')

        ax.set_title("Mantel Correlogram")
        ax.set_xlabel("Distance class index")
        ax.set_ylabel("Mantel correlation statistic")
        return fig


class Mantel(CorrelationStats):
    """Class for the Mantel matrix correlation statistical method.

    This class provides the functionality to run a Mantel analysis on two
    distance matrices. A Mantel test essentially computes the Pearson
    correlation between the two distance matrices.
    """

    def __init__(self, dm1, dm2, tail_type='two sided',
                 suppress_symmetry_and_hollowness_check=False):
        """Constructs a new Mantel instance.

        WARNING: Only symmetric, hollow distance matrices may be used as input.
        Asymmetric distance matrices, such as those obtained by the UniFrac
        Gain metric (i.e. beta_diversity.py -m unifrac_g), should not be used
        as input.

        Arguments:
            dm1 - first DistanceMatrix object to be compared
            dm2 - second DistanceMatrix object to be compared
            tail_type - the type of Mantel test to perform (i.e. hypothesis
                test). Can be "two sided", "less", or "greater"
            suppress_symmetry_and_hollowness_check - by default, the input
                distance matrices will be checked for symmetry and hollowness.
                It is recommended to leave this check in place for safety, as
                the check is fairly fast. However, if you *know* you have
                symmetric and hollow distance matrices, you can disable this
                check for small performance gains on extremely large distance
                matrices
        """
        super(Mantel, self).__init__([dm1, dm2], num_dms=2, min_dm_size=3,
                suppress_symmetry_and_hollowness_check=\
                suppress_symmetry_and_hollowness_check)
        self.TailType = tail_type

    @property
    def TailType(self):
        """Returns the tail type being used for the Mantel test."""
        return self._tail_type

    @TailType.setter
    def TailType(self, tail_type):
        """Sets the tail type that will be used for the Mantel test.

        Arguments:
            tail_type - the tail type to use when calculating the p-value.
                Valid types are 'two sided', 'less', or 'greater'.
        """
        if tail_type not in ("two sided", "greater", "less"):
            raise ValueError("Unrecognized alternative hypothesis (tail "
                             "type). Must be either 'two sided', 'greater', "
                             "or 'less'.")
        self._tail_type = tail_type

    def __call__(self, num_perms=999):
        """Runs a Mantel test over the current distance matrices.

        Returns a dict containing the results. The following keys are set:
            method_name - name of the statistical method
            dm1 - the first DistanceMatrix instance that was used
            dm2 - the second DistanceMatrix instance that was used
            num_perms - the number of permutations used to compute the p-value
            p_value - the p-value computed by the test
            r_value - the Mantel r statistic computed by the test
            perm_stats - a list of Mantel r statistics, one for each
                permutation
            tail_type - the type of Mantel test performed

        Arguments:
            num_perms - the number of times to permute the distance matrix
                while calculating the p-value

        Note: R's mantel function will always perform a one-sided test (type
        'greater'), so the p-values may differ from R unless you explicitly
        specify the tail type of 'greater'.
        """
        m1, m2 = self.DistanceMatrices
        alt = self.TailType

        # We suppress the symmetric and hollowness check since that is
        # guaranteed to have already happened when the distance matrices were
        # set (we don't need to do it a second time here).
        results = mantel_test(m1.DataMatrix, m2.DataMatrix, num_perms, alt=alt,
                              suppress_symmetry_and_hollowness_check=True)

        resultsDict = super(Mantel, self).__call__(num_perms)
        resultsDict['method_name'] = "Mantel"
        resultsDict['dm1'] = self.DistanceMatrices[0]
        resultsDict['dm2'] = self.DistanceMatrices[1]
        resultsDict['num_perms'] = num_perms
        resultsDict['p_value'] = results[0]
        resultsDict['r_value'] = results[1]
        resultsDict['perm_stats'] = results[2]
        resultsDict['tail_type'] = self.TailType

        return resultsDict


class PartialMantel(CorrelationStats):
    """Class for the partial Mantel matrix correlation statistical method.

    This class provides the functionality to run a partial Mantel analysis on
    three distance matrices. A partial Mantel test essentially computes the
    Pearson correlation between two distance matrices after first controlling
    for the effects of a third distance matrix (the control matrix).
    """

    def __init__(self, dm1, dm2, cdm,
                 suppress_symmetry_and_hollowness_check=False):
        """Constructs a new PartialMantel instance.

        WARNING: Only symmetric, hollow distance matrices may be used as input.
        Asymmetric distance matrices, such as those obtained by the UniFrac
        Gain metric (i.e. beta_diversity.py -m unifrac_g), should not be used
        as input.

        Arguments:
            dm1 - first DistanceMatrix object to be compared
            dm2 - second DistanceMatrix object to be compared
            cdm - the control DistanceMatrix object
            suppress_symmetry_and_hollowness_check - by default, the input
                distance matrices will be checked for symmetry and hollowness.
                It is recommended to leave this check in place for safety, as
                the check is fairly fast. However, if you *know* you have
                symmetric and hollow distance matrices, you can disable this
                check for small performance gains on extremely large distance
                matrices
        """
        super(PartialMantel, self).__init__([dm1, dm2, cdm], num_dms=3,
                min_dm_size=3, suppress_symmetry_and_hollowness_check=\
                suppress_symmetry_and_hollowness_check)

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
        corr = lambda rxy, rxz, ryz: (rxy - rxz*ryz)/(sqrt(1 -
                                      rxz**2)*sqrt(1 - ryz**2))
        # Load initial/placeholder values in the results dictionary.
        res['method_name'] = 'Partial Mantel'
        res['mantel_r'] = None
        res['mantel_p'] = None

        dm1, dm2, cdm = self.DistanceMatrices
        dm1_flat = dm1.flatten()
        dm2_flat = dm2.flatten()
        cdm_flat = cdm.flatten()

        # Get the initial r-values before permuting.
        rval1 = pearson(dm1_flat, dm2_flat)
        rval2 = pearson(dm1_flat, cdm_flat)
        rval3 = pearson(dm2_flat, cdm_flat)

        # Calculate the original test statistic (r-value).
        orig_stat = corr(rval1, rval2, rval3)

        # Calculate permuted r-values and p-values, storing
        # them for use in the calculation of the final statistic.
        perm_stats = []
        numerator = 0
        for i in range(0, num_perms):
            # Permute the first distance matrix and calculate new
            # r and p-values.
            p1 = permute_2d(dm1, permutation(dm1.Size))
            dm1_perm = DistanceMatrix(p1, dm1.SampleIds, dm1.SampleIds)
            dm1_perm_flat = dm1_perm.flatten()
            rval1 = pearson(dm1_perm_flat, dm2_flat)
            rval2 = pearson(dm1_perm_flat, cdm_flat)
            perm_stats.append(corr(rval1, rval2, rval3))

            # Sum the permuted statistics for calculation of the final
            # statistic.
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
       {'subject1':{'firmicutes-abundance':[0.45,0.55],
                   'bacteroidetes-abundace':[0.22,0.11]},
        'subject2':{'firmicutes-abundance':[0.11,0.52],
                   'bacteroidetes-abundace':[0.28,0.21]},
         ...
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
        raise ValueError, ("Only two state values can be provided. "
        "Support currently exists only for pre/post experimental design.")
    
    # create the output directory if it doesn't already exist
    create_dir(output_dir)
    
    num_analysis_categories = len(analysis_categories)
    x_values = range(len(state_values))
    
    paired_difference_output_fp = \
     join(output_dir,'paired_difference_comparisons.txt')
    paired_difference_output_f = open(paired_difference_output_fp,'w')
    # write header line to output file
    paired_difference_output_f.write(
     "#Metadata category\tNum differences (i.e., n)\tMean difference\t"
     "Median difference\tt one sample\tt one sample parametric p-value\t"
     "t one sample parametric p-value (Bonferroni-corrected)\n")
    
    paired_difference_t_test_results = []
    # initiate list of output file paths to return 
    output_fps = [paired_difference_output_fp]

    for category_number, analysis_category in enumerate(analysis_categories):
        personal_ids_to_state_metadatum = personal_ids_to_state_values[analysis_category]
        plot_output_fp = join(output_dir,'%s.pdf' % analysis_category.replace(' ','-'))
        fig = figure()
        axes = fig.add_axes([0.1, 0.1, 0.8, 0.8])
        
        # initialize a list to store the distribution of changes 
        # with state change
        differences = []
        
        for pid, data in personal_ids_to_state_metadatum.items():
            if None in data:
                # if any of the data points are missing, skip this 
                # individual
                continue
            else:
                # otherwise compute the difference between the ending
                # and starting state
                differences.append(data[1] - data[0])
                # and plot the start and stop values as a line
                axes.plot(x_values,data,"black",linewidth=0.5)
        
        # run stats for current analysis category
        t_one_sample_results = t_one_sample(differences)
        t = t_one_sample_results[0]
        p_value = t_one_sample_results[1]
        bonferroni_p_value = min([p_value * num_analysis_categories,1.0])
        paired_difference_t_test_results.append([analysis_category,
                                        len(differences),
                                        mean(differences),
                                        median(differences),
                                        t,
                                        p_value,
                                        bonferroni_p_value])
        
        # Finalize plot for current analysis category
        axes.set_ylabel(analysis_category)
        axes.set_xticks(range(len(state_values)))
        axes.set_xticklabels(state_values)
        axes.set_ylim(ymin=ymin,ymax=ymax)
        fig.savefig(plot_output_fp)
        output_fps.append(plot_output_fp)
    
    # sort output by uncorrected p-value and write results
    # to file
    paired_difference_t_test_results.sort(key=lambda x: x[5])
    for r in paired_difference_t_test_results:
        paired_difference_output_f.write('\t'.join(map(str,r)))
        paired_difference_output_f.write('\n')
    paired_difference_output_f.close()
    
    return output_fps
