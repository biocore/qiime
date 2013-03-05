#!/usr/bin/env python
"""Provides standard statistical tests. Tests produce statistic and P-value.
"""
from __future__ import division
import warnings
from cogent.maths.stats.distribution import chi_high, z_low, z_high, zprob, \
    t_high, t_low, tprob, f_high, f_low, fprob, binomial_high, binomial_low, \
    ndtri
from cogent.maths.stats.special import lgam, log_one_minus, one_minus_exp,\
    MACHEP
from cogent.maths.stats.ks import psmirnov2x, pkstwo
from cogent.maths.stats.kendall import pkendall, kendalls_tau
from cogent.maths.stats.special import Gamma

from numpy import absolute, arctanh, array, asarray, concatenate, transpose, \
        ravel, take, nonzero, log, sum, mean, cov, corrcoef, fabs, any, \
        reshape, tanh, clip, nan, isnan, isinf, sqrt, trace, exp, \
        median as _median, zeros, ones
        #, std - currently incorrect
from numpy.random import permutation, randint
from cogent.maths.stats.util import Numbers
from operator import add
from random import choice

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Gavin Huttley", "Rob Knight", "Catherine Lozupone",
               "Sandra Smit", "Micah Hamady", "Daniel McDonald",
               "Greg Caporaso", "Jai Ram Rideout", "Michael Dwan"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"

class IndexOrValueError(IndexError, ValueError): pass

var = cov   #cov will calculate variance if called on a vector

def std_(x, axis=None):
    """Returns standard deviations by axis (similiar to numpy.std)

    The result is unbiased, matching the result from MLab.std
    """
    x = asarray(x)
    
    if axis is None:
        d = x - mean(x)
        return sqrt(sum(d**2)/(len(x)-1))
    elif axis == 0:
        result = []
        for col in range(x.shape[1]):
            vals = x[:,col]
            d = vals - mean(vals)
            result.append(sqrt(sum(d**2)/(len(x)-1)))
        return result
    elif axis == 1:
        result = []
        for row in range(x.shape[0]):
            vals = x[row,:]
            d = vals - mean(vals)
            result.append(sqrt(sum(d**2)/(len(x)-1)))
        return result
    else:
        raise ValueError, "axis out of bounds"
    
# tested only by std
def var(x, axis=None):
    """Returns unbiased standard deviations over given axis.

    Similar with numpy.std, except that it is unbiased. (var = SS/n-1)

    x: a float ndarray or asarray(x) is a float ndarray.
    axis=None: computed for the flattened array by default, or compute along an
     integer axis.

    Implementation Notes:
    Change the SS calculation from:
        SumSq(x-x_bar) to SumSq(x) - SqSum(x)/n
        See p. 37 of Zar (1999) Biostatistical Analysis.
    """
    x = asarray(x)
    #figure out sample size along the axis
    if axis is None:
        n = x.size
    else:
        n = x.shape[axis]
    #compute the sum of squares from the mean(s)
    sample_SS = sum(x**2, axis) - sum(x, axis)**2/n
    return sample_SS / (n-1)

def std(x, axis=None):
    """computed unbiased standard deviations along given axis or flat array.

    Similar with numpy.std, except that it is unbiased. (var = SS/n-1)

    x: a float ndarray or asarray(x) is a float ndarray.
    axis=None: computed for the flattened array by default, or compute along an
      given integer axis.
    """
    try:
        sample_variance = var(x, axis=axis)
    except IndexError, e: #just to avoid breaking the old test code
        raise IndexOrValueError(e)
    return sqrt(sample_variance)

def median(m, axis=None):
    """Returns medians by axis (similiar to numpy.mean)

    numpy.median does not except an axis parameter. Is safe for substition for 
    numpy.median
    """
    median_vals = []
    rows, cols = m.shape

    if axis is None:
        return _median(ravel(m))
    elif axis == 0:
        for col in range(cols):
            median_vals.append(_median(m[:,col]))
    elif axis == 1 or axis == -1:
        for row in range(rows):
            median_vals.append(_median(m[row,:]))
    else:
        raise ValueError, "axis(=%s) out of bounds" % axis

    return array(median_vals)

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
    #return 0 if table was empty
    if not n:
        return (0, 1)
    #raise error if any counts were negative
    if min(cells) < 0:
        raise ValueError, \
        "G_2_by_2 got negative cell counts(s): must all be >= 0."
    
    G = 0
    #Add x ln x for items, adding zero for items whose counts are zero
    for i in filter(None, cells):
        G += i * log(i)
    #Find totals for rows and cols
    ab = a + b
    cd = c + d
    ac = a + c
    bd = b + d
    rows_cols = [ab, cd, ac, bd]
    #exit if we are missing a row or column entirely: result counts as
    #never significant
    if min(rows_cols) == 0:
        return (0, 1)
    #Subtract x ln x for rows and cols
    for i in filter(None, rows_cols):
        G -= i * log(i)
    #Add x ln x for table
    G += n * log(n)
    #Result needs to be multiplied by 2 
    G *= 2

    #apply Williams correction
    if williams:
        q = 1 + ((  ( (n/ab) + (n/cd) ) -1 ) * ( ( (n/ac) + (n/bd) ) -1))/(6*n)
        G /= q

    p = chi_high(max(G,0), 1)
    
    #find which tail we were in if the test was directional
    if directional:
        is_high =  ((b == 0) or (d != 0 and (a/b > c/d)))
        p = tail(p, is_high)
        if not is_high:
            G = -G
    return G, p

def safe_sum_p_log_p(a, base=None):
    """Calculates p * log(p) safely for an array that may contain zeros."""
    flat = ravel(a)
    nz = take(flat, nonzero(flat)[0])
    logs = log(nz)
    if base:
        logs /= log(base)
    return sum(nz * logs, 0)

def G_ind(m, williams=False):
    """Returns G test for independence in an r x c table.
    
    Requires input data as a numpy array. From Sokal and Rohlf p 738.
    """
    f_ln_f_elements = safe_sum_p_log_p(m)
    f_ln_f_rows = safe_sum_p_log_p(sum(m,0))
    f_ln_f_cols = safe_sum_p_log_p(sum(m,1))
    tot = sum(ravel(m))
    f_ln_f_table = tot * log(tot)

    df = (len(m)-1) * (len(m[0])-1)
    G = 2*(f_ln_f_elements-f_ln_f_rows-f_ln_f_cols+f_ln_f_table)
    if williams:
        q = 1+((tot*sum(1.0/sum(m,1))-1)*(tot*sum(1.0/sum(m,0))-1)/ \
            (6*tot*df))
        G = G/q
    return G, chi_high(max(G,0), df)


def calc_contingency_expected(matrix):
        """Calculates expected frequencies from a table of observed frequencies

        The input matrix is a dict2D object and represents a frequency table
        with different variables in the rows and columns. (observed
        frequencies as values)
                                
        The expected value is calculated with the following equation:
            Expected = row_total x column_total / overall_total
                                            
        The returned matrix (dict2D) has lists of the observed and the 
        expected frequency as values
        """
        #transpose matrix for calculating column totals
        t_matrix = matrix.copy()
        t_matrix.transpose()

        overall_total = sum(list(matrix.Items))
        #make new matrix for storing results
        result = matrix.copy()

        #populate result with expected values
        for row in matrix:
            row_sum = sum(matrix[row].values())
            for item in matrix[row]:
                column_sum = sum(t_matrix[item].values())
                #calculate expected frequency
                Expected = (row_sum * column_sum)/overall_total
                result[row][item] = [result[row][item]]
                result[row][item].append(Expected)
        return result

def G_fit(obs, exp, williams=1):
    """G test for fit between two lists of counts.

    Usage: test, prob = G_fit(obs, exp, williams)
    
    obs and exp are two lists of numbers.
    williams is a boolean stating whether to do the Williams correction.
    
    SUM(2 f(obs)ln (f(obs)/f(exp)))
    
    See Sokal and Rohlf chapter 17.
    """
    k = len(obs)
    if k != len(exp):
        raise ValueError, "G_fit requires two lists of equal length."
    G = 0
    n = 0
    
    for o, e in zip(obs, exp):
        if o < 0:
            raise ValueError, \
            "G_fit requires all observed values to be positive."
        if e <= 0:
            raise ZeroExpectedError, \
            "G_fit requires all expected values to be positive."
        if o:   #if o is zero, o * log(o/e) must be zero as well.
            G += o * log(o/e)
            n += o
    
    G *= 2
    if williams:
        q = 1 + (k + 1)/(6*n)
        G /= q

    return G, chi_high(G, k - 1)

def G_fit_from_Dict2D(data):
    """G test for fit on a Dict2D

    data is a dict2D. Values are a list containing the observed
    and expected frequencies (can be created with calc_contingency_expected)
    """
    obs_counts = []
    exp_counts = []
    for item in data.Items:
        if len(item) == 2:
            obs_counts.append(item[0])
            exp_counts.append(item[1])
    g_val, prob = G_fit(obs_counts, exp_counts)
    return g_val, prob

def chi_square_from_Dict2D(data):
    """Chi Square test on a Dict2D

    data is a Dict2D. The values are a list of the observed (O)
    and expected (E) frequencies,(can be created with calc_contingency_expected)

    The chi-square value (test) is the sum of (O-E)^2/E over the items in data

    degrees of freedom are calculated from data as:
    (r-1)*(c-1) if cols and rows are both > 1
    otherwise is just 1 - the # of rows or columns
    (whichever is greater than 1)
    
    """
    test =  sum([((item[0] - item[1]) * (item[0] - item[1]))/item[1] \
                   for item in data.Items])
    num_rows = len(data)
    num_cols = len([col for col in data.Cols])
    if num_rows == 1:
        df = num_cols - 1
    elif num_cols == 1:
        df = num_rows - 1
    elif num_rows == 0 or num_cols == 0:
        raise ValueError, "data matrix must have data"
    else:
        df = (len(data) - 1) * (len([col for col in data.Cols]) - 1)
    
    return test, chi_high(test, df)
    

def likelihoods(d_given_h, priors):
    """Calculate likelihoods through marginalization, given Pr(D|H) and priors.

    Usage: scores = likelihoods(d_given_h, priors)

    d_given_h and priors are equal-length lists of probabilities. Returns
    a list of the same length of numbers (not probabilities).
    """
    #check that the lists of Pr(D|H_i) and priors are equal
    length = len(d_given_h)
    if length != len(priors):
        raise ValueError, "Lists not equal lengths."
    #find weighted sum of Pr(H_i) * Pr(D|H_i)
    wt_sum = 0
    for d, p in zip(d_given_h, priors):
        wt_sum += d * p
    #divide each Pr(D|H_i) by the weighted sum and multiply by its prior
    #to get its likelihood
    return [d/wt_sum for d in d_given_h]

def posteriors(likelihoods, priors):
    """Calculate posterior probabilities given priors and likelihoods.

    Usage: probabilities = posteriors(likelihoods, priors)

    likelihoods is a list of numbers. priors is a list of probabilities.
    Returns a list of probabilities (0-1).
    """
    #Check that there is a prior for each likelihood
    if len(likelihoods) != len(priors):
        raise ValueError, "Lists not equal lengths."
    #Posterior probability is defined as prior * likelihood
    return [l * p for l, p in zip(likelihoods, priors)]

def bayes_updates(ds_given_h, priors = None):
    """Successively apply lists of Pr(D|H) to get Pr(H|D) by marginalization.

    Usage: final_probs = bayes_updates(ds_given_h, [priors])

    ds_given_h is a list (for each form of evidence) of lists of probabilities.
    priors is optionally a list of the prior probabilities.
    Returns a list of posterior probabilities.
    """
    try:
        first_list = ds_given_h[0]
        length = len(first_list)
        #calculate flat prior if none was passed
        if not priors:
            priors = [1/length] * length
        #apply each form of data to the priors to get posterior probabilities
        for index, d in enumerate(ds_given_h):
            #first, ignore the form of data if all the d's are the same
            all_the_same = True
            first_element = d[0]
            for i in d:
                if i != first_element:
                    all_the_same = False
                    break
            if not all_the_same:    #probabilities won't change
                if len(d) != length:
                    raise ValueError, "bayes_updates requires equal-length lists."
                liks = likelihoods(d, priors)
                pr = posteriors(liks, priors)
                priors = pr
        return priors #posteriors after last calculation are 'priors' for next
    #return column of zeroes if anything went wrong, e.g. if the sum of one of
    #the ds_given_h is zero.
    except (ZeroDivisionError, FloatingPointError):
        return [0] * length

def t_paired (a,b, tails=None, exp_diff=0):
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
    n = len(a)
    if n != len(b):
        raise ValueError, 'Unequal length lists in ttest_paired.'
    try:
        diffs = array(a) - array(b)
        return t_one_sample(diffs, popmean=exp_diff, tails=tails)
    except (ZeroDivisionError, ValueError, AttributeError, TypeError, \
        FloatingPointError):
        return (None, None)
    

def t_one_sample(a,popmean=0, tails=None):
    """Returns t for ONE group of scores a, given a population mean.

    Usage:   t, prob = t_one_sample(a, popmean, tails)

    t is a float; prob is a probability.
    a should support Mean, StandardDeviation, and Count.
    popmean should be the expected mean; 0 by default.
    tails should be None (default), 'high', or 'low'.
"""
    try:
        n = len(a)
        t = (mean(a) - popmean)/(std(a)/sqrt(n))
    except (ZeroDivisionError, ValueError, AttributeError, TypeError, \
        FloatingPointError):
        return None, None
    if isnan(t) or isinf(t):
        return None, None

    prob = t_tailed_prob(t, n-1, tails)
    return t, prob


def t_two_sample(a, b, tails=None, exp_diff=0, none_on_zero_variance=True):
    """Returns t, prob for two INDEPENDENT samples of scores a, and b.  
    
    From Sokal and Rohlf, p 223.  

    Usage:   t, prob = t_two_sample(a,b, tails, exp_diff)

    t is a float; prob is a probability.
    a and b should be sequences of observations (numbers). Need not be equal
        lengths.
    tails should be None (default), 'high', or 'low'.
    exp_diff should be the expected difference in means (a-b); 0 by default.
    none_on_zero_variance: if True, will return (None,None) if both a and b
        have zero variance (e.g. a=[1,1,1] and b=[2,2,2]). If False, the
        following values will be returned:

            Two-tailed test (tails=None):
                a < b: (-inf,0.0)
                a > b: (+inf,0.0)

            One-tailed 'high':
                a < b: (-inf,1.0)
                a > b: (+inf,0.0)

            One-tailed 'low':
                a < b: (-inf,0.0)
                a > b: (+inf,1.0)

        If a and b both have no variance and have the same single value (e.g.
        a=[1,1,1] and b=[1,1,1]), (None,None) will always be returned.
    """
    if tails is not None and tails != 'high' and tails != 'low':
        raise ValueError("Invalid tail type '%s'. Must be either None, "
                         "'high', or 'low'." % tails)

    try:
        #see if we need to back off to the single-observation for single-item
        #groups
        n1 = len(a)
        if n1 < 2:
            return t_one_observation(sum(a), b, tails, exp_diff,
                    none_on_zero_variance=none_on_zero_variance)

        n2 = len(b)
        if n2 < 2:
            t, prob = t_one_observation(sum(b), a, reverse_tails(tails),
                    exp_diff, none_on_zero_variance=none_on_zero_variance)

            # Negate the t-statistic because we swapped the order of the inputs
            # in the t_one_observation call, as well as tails.
            if t != 0:
                t = -t

            return (t, prob)

        #otherwise, calculate things properly
        x1 = mean(a)
        x2 = mean(b)
        var1 = var(a)
        var2 = var(b)

        if var1 == 0 and var2 == 0:
            # Both lists do not vary.
            if x1 == x2 or none_on_zero_variance:
                result = (None, None)
            else:
                result = _t_test_no_variance(x1, x2, tails)
        else:
            # At least one list varies.
            df = n1+n2-2
            svar = ((n1-1)*var1 + (n2-1)*var2)/df
            t = (x1-x2-exp_diff)/sqrt(svar*(1/n1 + 1/n2))

            if isnan(t) or isinf(t):
                result = (None, None)
            else:
                prob = t_tailed_prob(t, df, tails)
                result = (t, prob)
    except (ZeroDivisionError, ValueError, AttributeError, TypeError,
            FloatingPointError), e:
        #invalidate if the sample sizes are wrong, the values aren't numeric or
        #aren't present, etc.
        result = (None, None)

    return result

def _t_test_no_variance(mean1, mean2, tails):
    """Handles case where two distributions have no variance."""
    if tails is not None and tails != 'high' and tails != 'low':
        raise ValueError("Invalid tail type '%s'. Must be either None, "
                         "'high', or 'low'." % tails)

    if tails is None:
        if mean1 < mean2:
            result = (float('-inf'), 0.0)
        else:
            result = (float('inf'), 0.0)
    elif tails == 'high':
        if mean1 < mean2:
            result = (float('-inf'), 1.0)
        else:
            result = (float('inf'), 0.0)
    else:
        if mean1 < mean2:
            result = (float('-inf'), 0.0)
        else:
            result = (float('inf'), 1.0)

    return result

def mc_t_two_sample(x_items, y_items, tails=None, permutations=999,
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
            and the nonparametric p-value will be None
        exp_diff - the expected difference in means (x_items - y_items)
    """
    if tails is not None and tails != 'high' and tails != 'low':
        raise ValueError("Invalid tail type '%s'. Must be either None, "
                         "'high', or 'low'." % tails)
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
                                      exp_diff=exp_diff,
                                      none_on_zero_variance=False)

    # Only perform the Monte Carlo test if we got a sane answer back from the
    # initial t-test and we have been specified permutations.
    nonparam_p_val = None
    perm_t_stats = []
    if permutations > 0 and obs_t is not None and param_p_val is not None:
        # Permute observations between x_items and y_items the specified number
        # of times.
        perm_x_items, perm_y_items = _permute_observations(x_items, y_items,
                                                           permutations)
        perm_t_stats = [t_two_sample(perm_x_items[n], perm_y_items[n],
                                     tails=tails, exp_diff=exp_diff,
                                     none_on_zero_variance=False)[0]
                        for n in range(permutations)]

        # Compute nonparametric p-value based on the permuted t-test results.
        if tails is None:
            better = (absolute(array(perm_t_stats)) >= absolute(obs_t)).sum()
        elif tails == 'low':
            better = (array(perm_t_stats) <= obs_t).sum()
        elif tails == 'high':
            better = (array(perm_t_stats) >= obs_t).sum()
        nonparam_p_val = (better + 1) / (permutations + 1)

    return obs_t, param_p_val, perm_t_stats, nonparam_p_val

def _permute_observations(x_items, y_items, permutations,
                          permute_f=permutation):
    """Returns permuted versions of the sequences of observations.

    Values are permuted between x_items and y_items (i.e. shuffled between the
    two input sequences of observations).

    This code is based on Jeremy Widmann's
    qiime.make_distance_histograms.permute_between_groups code.
    """
    num_x = len(x_items)
    num_y = len(y_items)
    num_total_obs = num_x + num_y
    combined_obs = concatenate((x_items, y_items))

    # Generate a list of all permutations.
    perms = [permute_f(num_total_obs) for i in range(permutations)]

    # Use random permutations to split into groups.
    rand_xs = [combined_obs[perm[:num_x]] for perm in perms]
    rand_ys = [combined_obs[perm[num_x:num_total_obs]] for perm in perms]
    return rand_xs, rand_ys

def t_one_observation(x, sample, tails=None, exp_diff=0,
                      none_on_zero_variance=True):
    """Returns t-test for significance of single observation versus a sample.
    
    Equation for 1-observation t (Sokal and Rohlf 1995 p 228):
    t = obs - mean - exp_diff / (var * sqrt((n+1)/n)) 
    df = n - 1

    none_on_zero_variance: see t_two_sample for details. If sample has no
        variance and its single value is the same as x (e.g. x=1 and
        sample=[1,1,1]), (None,None) will always be returned
    """
    try:
        sample_mean = mean(sample)
        sample_std = std(sample)

        if sample_std == 0:
            # The list does not vary.
            if sample_mean == x or none_on_zero_variance:
                result = (None, None)
            else:
                result = _t_test_no_variance(x, sample_mean, tails)
        else:
            # The list varies.
            n = len(sample)
            t = (x - sample_mean - exp_diff)/sample_std/sqrt((n+1)/n)
            prob = t_tailed_prob(t, n-1, tails)
            result = (t, prob)
    except (ZeroDivisionError, ValueError, AttributeError, TypeError,
            FloatingPointError):
        result = (None, None)

    return result

def pearson(x_items, y_items):
    """Returns Pearson's product moment correlation coefficient.
    
    This will always be a value between -1.0 and +1.0. x_items and y_items must
    be the same length, and cannot have fewer than 2 elements each. If one or
    both of the input vectors do not have any variation, the return value will
    be 0.0.

    Arguments:
        x_items - the first list of observations
        y_items - the second list of observations
    """
    x_items, y_items = array(x_items), array(y_items)

    if len(x_items) != len(y_items):
        raise ValueError("The length of the two vectors must be the same in "
                         "order to calculate the Pearson correlation "
                         "coefficient.")
    if len(x_items) < 2:
        raise ValueError("The two vectors must both contain at least 2 "
                "elements. The vectors are of length %d." % len(x_items))

    sum_x = sum(x_items)
    sum_y = sum(y_items)
    sum_x_sq = sum(x_items*x_items)
    sum_y_sq = sum(y_items*y_items)
    sum_xy = sum(x_items*y_items)
    n = len(x_items)

    try:
        r = 1.0 * ((n * sum_xy) - (sum_x * sum_y)) / \
           (sqrt((n * sum_x_sq)-(sum_x*sum_x))*sqrt((n*sum_y_sq)-(sum_y*sum_y)))
    except (ZeroDivisionError, ValueError, FloatingPointError): #no variation
        r = 0.0
    #check we didn't get a naughty value for r due to rounding error
    if r > 1.0:
        r = 1.0
    elif r < -1.0:
        r = -1.0
    return r

def spearman(x_items, y_items):
    """Returns Spearman's rho.

    This will always be a value between -1.0 and +1.0. x_items and y_items must
    be the same length, and cannot have fewer than 2 elements each. If one or
    both of the input vectors do not have any variation, the return value will
    be 0.0.

    Arguments:
        x_items - the first list of observations
        y_items - the second list of observations
    """
    x_items, y_items = array(x_items), array(y_items)

    if len(x_items) != len(y_items):
        raise ValueError("The length of the two vectors must be the same in "
                         "order to calculate Spearman's rho.")
    if len(x_items) < 2:
        raise ValueError("The two vectors must both contain at least 2 "
                "elements. The vectors are of length %d." % len(x_items))

    # Rank the two input vectors.
    rank1, ties1 = _get_rank(x_items)
    rank2, ties2 = _get_rank(y_items)

    if ties1 == 0 and ties2 == 0:
        n = len(rank1)
        sum_sqr = sum([(x-y)**2 for x,y in zip(rank1,rank2)])
        rho = 1 - (6*sum_sqr/(n*(n**2 - 1)))
    else:
        avg = lambda x: sum(x)/len(x)

        x_bar = avg(rank1)
        y_bar = avg(rank2)

        numerator = sum([(x-x_bar)*(y-y_bar) for x,y in zip(rank1, rank2)])
        denominator = sqrt(sum([(x-x_bar)**2 for x in rank1])*
                           sum([(y-y_bar)**2 for y in rank2]))

        # Calculate rho. Handle the case when there is no variation in one or
        # both of the input vectors.
        if denominator == 0.0:
            rho = 0.0
        else:
            rho = numerator/denominator
    return rho

def _get_rank(data):
    """Ranks the elements of a list. Used in Spearman correlation."""
    indices = range(len(data))
    ranks = range(1,len(data)+1)
    indices.sort(key=lambda index:data[index])
    ranks.sort(key=lambda index:indices[index-1])
    data_len = len(data)
    i = 0
    ties = 0
    while i < data_len:
        j = i + 1
        val = data[indices[i]]
        try:
            val += 0
        except TypeError:
            raise(TypeError)

        while j < data_len and data[indices[j]] == val:
            j += 1
        dup_ranks = j - i
        val = float(ranks[indices[i]]) + (dup_ranks-1)/2.0
        for k in range(i, i+dup_ranks):
            ranks[indices[k]] = val
        i += dup_ranks
        ties += dup_ranks-1
    return ranks, ties

def correlation(x_items, y_items):
    """Returns Pearson correlation between x and y, and its significance.

    WARNING: x_items and y_items must be same length!

    This function is retained for backwards-compatibility. Please use
    correlation_test() for more control over how the test is performed.
    """
    return correlation_test(x_items, y_items, method='pearson', tails=None,
                            permutations=0)[:2]

def correlation_test(x_items, y_items, method='pearson', tails=None,
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
    if tails is not None and tails != 'high' and tails != 'low':
        raise ValueError("Invalid tail type '%s'. Must be either None, "
                         "'high', or 'low'." % tails)
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
            parametric_p_val = t_tailed_prob(t, df, tails)
        except (ZeroDivisionError, FloatingPointError):
            # r/rho was presumably 1.
            parametric_p_val = 0

    # Perform the nonparametric test.
    permuted_corr_coeffs = []
    nonparametric_p_val = None
    better = 0
    for i in range(permutations):
        permuted_y_items = y_items[permutation(n)]
        permuted_corr_coeff = corr_fn(x_items, permuted_y_items)
        permuted_corr_coeffs.append(permuted_corr_coeff)

        if tails is None:
            if abs(permuted_corr_coeff) >= abs(corr_coeff):
                better += 1
        elif tails == 'high':
            if permuted_corr_coeff >= corr_coeff:
                better += 1
        elif tails == 'low':
            if permuted_corr_coeff <= corr_coeff:
                better += 1
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
            ci_low = tanh(arctanh(corr_coeff) - (z_crit / sqrt(n - 3)))
            ci_high = tanh(arctanh(corr_coeff) + (z_crit / sqrt(n - 3)))
        except (ZeroDivisionError, FloatingPointError):
            # r/rho was presumably 1 or -1. Match what R does in this case.
            ci_low, ci_high = corr_coeff, corr_coeff

    return (corr_coeff, parametric_p_val, permuted_corr_coeffs,
            nonparametric_p_val, (ci_low, ci_high))

def correlation_matrix(series, as_rows=True):
    """Returns pairwise correlations between each pair of series.
    """
    return corrcoef(series, rowvar=as_rows)
    #unused codes below
    if as_rows:
        return corrcoef(transpose(array(series)))
    else:
        return corrcoef(array(series))

def regress(x, y):
    """Returns coefficients to the regression line "y=ax+b" from x[] and y[].  

    Specifically, returns (slope, intercept) as a tuple from the regression of 
    y on x, minimizing the error in y assuming that x is precisely known.

    Basically, it solves 
        Sxx a + Sx b = Sxy
         Sx a +  N b = Sy
    where Sxy = \sum_i x_i y_i, Sx = \sum_i x_i, and Sy = \sum_i y_i.  The
    solution is
        a = (Sxy N - Sy Sx)/det
        b = (Sxx Sy - Sx Sxy)/det
    where det = Sxx N - Sx^2.  In addition,
        Var|a| = s^2 |Sxx Sx|^-1 = s^2 | N  -Sx| / det
           |b|       |Sx  N |          |-Sx Sxx|
        s^2 = {\sum_i (y_i - \hat{y_i})^2 \over N-2}
            = {\sum_i (y_i - ax_i - b)^2 \over N-2}
            = residual / (N-2)
        R^2 = 1 - {\sum_i (y_i - \hat{y_i})^2 \over \sum_i (y_i - \mean{y})^2}
            = 1 - residual/meanerror

    Adapted from the following URL:
    http://www.python.org/topics/scicomp/recipes_in_python.html
    """
    x, y = array(x,'Float64'), array(y,'Float64')
    N = len(x)
    Sx = sum(x)
    Sy = sum(y)
    Sxx = sum(x*x)
    Syy = sum(y*y)
    Sxy = sum(x*y)
    det = Sxx * N - Sx * Sx
    return (Sxy * N - Sy * Sx)/det, (Sxx * Sy - Sx * Sxy)/det

def regress_origin(x,y):
    """Returns coefficients to regression "y=ax+b" passing through origin.

    Requires vectors x and y of same length.
    See p. 351 of Zar (1999) Biostatistical Analysis.

    returns slope, intercept as a tuple.
    """
    x, y = array(x,'Float64'), array(y,'Float64')
    return sum(x*y)/sum(x*x), 0

def regress_R2(x, y):
    """Returns the R^2 value for the regression of x and y
    
    Used the method explained on pg 334 ofJ.H. Zar, Biostatistical analysis,
    fourth edition. 1999
    """
    slope, intercept = regress(x,y)
    coords = zip(x, y)
    Sx = Sy = Syy = SXY =  0.0
    n = float(len(y))
    for x, y in coords:
        SXY += x*y
        Sx += x
        Sy += y
        Syy += y*y
    Sxy = SXY - (Sx*Sy)/n
    regSS = slope * Sxy
    totSS = Syy - ((Sy*Sy)/n)
    return regSS/totSS

def regress_residuals(x, y):
    """reports the residual (error) for each point from the linear regression"""
    slope, intercept = regress(x, y)
    coords = zip(x, y)
    residuals = []
    for x, y in coords:
        e = y - (slope * x) - intercept
        residuals.append(e)
    return residuals
    
def stdev_from_mean(x):
    """returns num standard deviations from the mean of each val in x[]"""
    x = array(x)
    return (x - mean(x))/std(x)

def regress_major(x, y):
    """Returns major-axis regression line of y on x.

    Use in cases where there is error in both x and y.
    """
    x, y = array(x), array(y)
    N = len(x)
    Sx = sum(x)
    Sy = sum(y)
    Sxx = sum(x*x)
    Syy = sum(y*y)
    Sxy = sum(x*y)
    var_y = (Syy-((Sy*Sy)/N))/(N-1)
    var_x = (Sxx-((Sx*Sx)/N))/(N-1)
    cov = (Sxy-((Sy*Sx)/N))/(N-1)
    mean_y = Sy/N
    mean_x = Sx/N
    D = sqrt((var_y + var_x)*(var_y + var_x) - 4*(var_y*var_x - (cov*cov)))
    eigen_1 = (var_y + var_x + D)/2
    slope = cov/(eigen_1 - var_y)
    intercept = mean_y - (mean_x * slope)
    return (slope, intercept)
    

def z_test(a, popmean=0, popstdev=1, tails=None):
    """Returns z and probability score for a single sample of items.

Calculates the z-score on ONE sample of items with mean x, given a population 
mean and standard deviation (parametric).

Usage:   z, prob = z_test(a, popmean, popstdev, tails)

z is a float; prob is a probability.
a is a sample with Mean and Count.
popmean should be the parametric population mean; 0 by default.
popstdev should be the parametric population standard deviation, 1 by default.
tails should be None (default), 'high', or 'low'.
""" 
    try:
        z = (mean(a) - popmean)/popstdev*sqrt(len(a))
        return z, z_tailed_prob(z, tails)
    except (ValueError, TypeError, ZeroDivisionError, AttributeError, \
        FloatingPointError):
        return None

def z_tailed_prob(z, tails):
    """Returns appropriate p-value for given z, depending on tails."""
    if tails == 'high':
        return z_high(z)
    elif tails == 'low':
        return z_low(z)
    else:
        return zprob(z)

def t_tailed_prob(t, df, tails):
    """Return appropriate p-value for given t and df, depending on tails."""
    if tails == 'high':
        return t_high(t, df)
    elif tails == 'low':
        return t_low(t, df)
    else:
        return tprob(t,df)

def reverse_tails(tails):
    """Swaps high for low or vice versa, leaving other values alone."""
    if tails == 'high':
        return 'low'
    elif tails == 'low':
        return 'high'
    else:
        return tails

def tail(prob, test):
    """If test is true, returns prob/2. Otherwise returns 1-(prob/2).
    """
    prob /= 2
    if test:
        return prob
    else:
        return 1 - prob

def combinations(n, k):
    """Returns the number of ways of choosing k items from n.
    """
    return exp(lgam(n+1) - lgam(k+1) - lgam(n-k+1))

def multiple_comparisons(p, n):
    """Corrects P-value for n multiple comparisons.
    
    Calculates directly if p is large and n is small; resorts to logs
    otherwise to avoid rounding (1-p) to 1
    """
    if p > 1e-6:   #if p is large and n small, calculate directly
        return 1 - (1-p)**n
    else:
        return one_minus_exp(-n * p)

def multiple_inverse(p_final, n):
    """Returns p_initial for desired p_final with n multiple comparisons.
    
    WARNING: multiple_inverse is not very reliable when p_final is very close
    to 1 (say, within 1e-4) since we then take the ratio of two very similar
    numbers.
    """
    return one_minus_exp(log_one_minus(p_final)/n)

def multiple_n(p_initial, p_final):
    """Returns number of comparisons such that p_initial maps to p_final.
    
    WARNING: not very accurate when p_final is very close to 1.
    """
    return log_one_minus(p_final)/log_one_minus(p_initial)

def fisher(probs):
    """Uses Fisher's method to combine multiple tests of a hypothesis.

    -2 * SUM(ln(P)) gives chi-squared distribution with 2n degrees of freedom.
    """
    try:
        return chi_high(-2 * sum(map(log, probs)), 2 * len(probs))
    except OverflowError, e:
        return 0.0 

def f_value(a,b):
    """Returns the num df, the denom df, and the F value.

    a, b: lists of values, must have Variance attribute (recommended to
    make them Numbers objects.

    The F value is always calculated by dividing the variance of a by the 
    variance of b, because R uses the same approach. In f_two_value it's 
    decided what p-value is returned, based on the relative sizes of the
    variances. 
    """
    if not any(a) or not any(b) or len(a) <= 1 or len(b) <= 1:
        raise ValueError, "Vectors should contain more than 1 element"
    F = var(a)/var(b)
    dfn = len(a)-1
    dfd = len(b)-1
    return dfn, dfd, F


def f_two_sample(a, b, tails=None):
    """Returns the dfn, dfd, F-value and probability for two samples a, and b.

    a and b: should be independent samples of scores. Should be lists of
    observations (numbers).

    tails should be None(default, two-sided test), 'high' or 'low'.

    This implementation returns the same results as the F test in R.
    """
    dfn, dfd, F = f_value(a, b)
    if tails == 'low':
        return dfn, dfd, F, f_low(dfn, dfd, F)
    elif tails == 'high':
        return dfn, dfd, F, f_high(dfn, dfd, F)
    else:
        if var(a) >= var(b):
            side='right'
        else:
            side='left'
        return dfn, dfd, F, fprob(dfn, dfd, F, side=side)

def ANOVA_one_way(a):
    """Performs a one way analysis of variance

    a is a list of lists of observed values. Each list is the values
    within a category. The analysis must include 2 or more categories(lists).
    the lists must have a Mean and variance attribute. Recommende to make
    the Numbers objects
    
    An F value is first calculated as the variance of the group means
    divided by the mean of the within-group variances.
    """
    group_means = []
    group_variances = []
    num_cases = 0
    all_vals = []
    for i in a:
        num_cases += len(i)
        group_means.append(i.Mean)
        group_variances.append(i.Variance * (len(i)-1))
        all_vals.extend(i)
    group_means = Numbers(group_means)
    #get within group variances (denominator)
    group_variances = Numbers(group_variances)
    dfd = num_cases - len(group_means)
    within_MS = sum(group_variances)/dfd
    #get between group variances (numerator)
    grand_mean = Numbers(all_vals).Mean
    between_MS = 0
    for i in a:
        diff = i.Mean - grand_mean
        diff_sq = diff * diff
        x = diff_sq * len(i)
        between_MS += x
    dfn = len(group_means) - 1
    between_MS = between_MS/dfn
    F = between_MS/within_MS
    return dfn, dfd, F, between_MS, within_MS, group_means, f_high(dfn, dfd, F)

def MonteCarloP(value, rand_values, tail = 'high'):
    """takes a true value and a list of random values as
        input and returns a p-value

    tail indicates which side of the distribution to look at:
        low = look for smaller values than expected by chance
        high = look for larger values than expected by chance
    """
    pop_size= len(rand_values)
    rand_values.sort()
    if tail == 'high':
        num_better = pop_size
        for i, curr_val in enumerate(rand_values):
            if value <= curr_val:
                num_better = i
                break
        p_val = 1-(num_better / pop_size)
    elif tail == 'low':
        num_better = pop_size
        for i, curr_val in enumerate(rand_values):
            if value < curr_val:
                num_better = i
                break
        p_val = num_better / pop_size
    return p_val

def sign_test(success, trials, alt="two sided"):
    """Returns the probability for the sign test.
    
    Arguments:
        - success: the number of successes
        - trials: the number of trials
        - alt: the alternate hypothesis, one of 'less', 'greater', 'two sided'
          (default).
    """
    lo = ["less", "lo", "lower", "l"]
    hi = ["greater", "hi", "high", "h", "g"]
    two = ["two sided", "2", 2, "two tailed", "two"]
    alt = alt.lower().strip()
    if alt in lo:
        p = binomial_low(success, trials, 0.5)
    elif alt in hi:
        success -= 1
        p = binomial_high(success, trials, 0.5)
    elif alt in two:
        success = min(success, trials-success)
        hi = 1 - binomial_high(success, trials, 0.5)
        lo = binomial_low(success, trials, 0.5)
        p = hi+lo
    else:
        raise RuntimeError("alternate [%s] not in %s" % (lo+hi+two))
    return p

def ks_test(x, y=None, alt="two sided", exact = None, warn_for_ties = True):
    """Returns the statistic and probability from the Kolmogorov-Smirnov test.
    
    Arguments:
        - x, y: vectors of numbers whose distributions are to be compared.
        - alt: the alternative hypothesis, default is 2-sided.
        - exact: whether to compute the exact probability
        - warn_for_ties: warns when values are tied. This should left at True
          unless a monte carlo variant, like ks_boot, is being used.
    
    Note the 1-sample cases are not implemented, although their cdf's are
    implemented in ks.py"""
    # translation from R 2.4
    num_x = len(x)
    num_y = None
    x = zip(x, zeros(len(x), int))
    lo = ["less", "lo", "lower", "l", "lt"]
    hi = ["greater", "hi", "high", "h", "g", "gt"]
    two = ["two sided", "2", 2, "two tailed", "two", "two.sided"]
    Pval = None
    if y is not None: # in anticipation of actually implementing the 1-sample cases
        num_y = len(y)
        y = zip(y, ones(len(y), int))
        n = num_x * num_y / (num_x + num_y)
        combined = x + y
        if len(set(combined)) < num_x + num_y:
            ties = True
        else:
            ties = False
        
        combined = array(combined, dtype=[('stat', float), ('sample', int)])
        combined.sort(order='stat')
        cumsum = zeros(combined.shape[0], float)
        scales = array([1/num_x, -1/num_y])
        indices = combined['sample']
        cumsum = scales.take(indices)
        cumsum = cumsum.cumsum()
        if exact == None:
            exact = num_x * num_y < 1e4
        
        if alt in two:
            stat = max(fabs(cumsum))
        elif alt in lo:
            stat = -cumsum.min()
        elif alt in hi:
            stat = cumsum.max()
        else:
            raise RuntimeError, "Unknown alt: %s" % alt
        if exact and alt in two and not ties:
            Pval = 1 - psmirnov2x(stat, num_x, num_y)
    else:
        raise NotImplementedError
    
    if Pval == None:
        if alt in two:
            Pval = 1 - pkstwo(sqrt(n) * stat)
        else:
            Pval = exp(-2 * n * stat**2)
    
    if ties and warn_for_ties:
        warnings.warn("Cannot compute correct KS probability with ties")
    
    try: # if numpy arrays were input, the Pval can be an array of len==1
        Pval = Pval[0]
    except (TypeError, IndexError):
        pass
    return stat, Pval

def _get_bootstrap_sample(x, y, num_reps):
    """yields num_reps random samples drawn with replacement from x and y"""
    combined = array(list(x) + list(y))
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

def ks_boot(x, y, alt = "two sided", num_reps=1000):
    """Monte Carlo (bootstrap) variant of the Kolmogorov-Smirnov test. Useful
    for when there are ties.
    
    Arguments:
        - x, y: vectors of numbers
        - alt: alternate hypothesis, as per ks_test
        - num_reps: number of replicates for the  bootstrap"""
    # based on the ks_boot method in the R Matching package
    # see http://sekhon.berkeley.edu/matching/
    # One important difference is I preserve the original sample sizes
    # instead of making them equal
    tol = MACHEP * 100
    combined = array(list(x) + list(y))
    observed_stat, _p = ks_test(x, y, exact=False, warn_for_ties=False)
    total_obs = len(combined)
    num_x = len(x)
    num_greater = 0
    for sampled_x, sampled_y in _get_bootstrap_sample(x, y, num_reps):
        sample_stat, _p = ks_test(sampled_x, sampled_y, alt=alt, exact=False,
                                    warn_for_ties=False)
        if sample_stat >= (observed_stat - tol):
            num_greater += 1
    return observed_stat, num_greater / num_reps

def _average_rank(start_rank, end_rank):
    ave_rank = sum(range(start_rank, end_rank+1)) / (1+end_rank-start_rank)
    return ave_rank

def mw_test(x, y):
    """computes the Mann-Whitney U statistic and the probability using the
    normal approximation"""
    if len(x) > len(y):
        x, y = y, x
    
    num_x = len(x)
    num_y = len(y)
    
    x = zip(x, zeros(len(x), int), zeros(len(x), int))
    y = zip(y, ones(len(y), int), zeros(len(y), int))
    combined = x+y
    combined = array(combined, dtype=[('stat', float), ('sample', int),
                                      ('rank', float)])
    combined.sort(order='stat')
    prev = None
    start = None
    ties = False
    T = 0.0
    for index in range(combined.shape[0]):
        value = combined['stat'][index]
        sample = combined['sample'][index]
        if value == prev and start is None:
            start = index
            continue
        
        if value != prev and start is not None:
            ties = True
            ave_rank = _average_rank(start, index)
            num_tied = index - start + 1
            T += (num_tied**3 - num_tied)
            for i in range(start-1, index):
                combined['rank'][i] = ave_rank
            start = None
        combined['rank'][index] = index+1
        prev = value
    
    if start is not None:
        ave_rank = _average_rank(start, index)
        num_tied = index - start + 2
        T += (num_tied**3 - num_tied)
        for i in range(start-1, index+1):
            combined['rank'][i] = ave_rank
    
    total = combined.shape[0]
    x_ranks_sum = sum(combined['rank'][i] for i in range(total) if combined['sample'][i] == 0)
    prod = num_x * num_y
    U1 = prod + (num_x * (num_x+1) / 2) - x_ranks_sum
    U2 = prod - U1
    U = max([U1, U2])
    numerator = U - prod / 2
    denominator = sqrt((prod / (total * (total-1)))*((total**3 - total - T)/12))
    z = (numerator/denominator)
    p = zprob(z)
    return U, p

def mw_boot(x, y, num_reps=1000):
    """Monte Carlo (bootstrap) variant of the Mann-Whitney test.
    
    Arguments:
        - x, y: vectors of numbers
        - num_reps: number of replicates for the  bootstrap
    
    Uses the same Monte-Carlo resampling code as kw_boot
    """
    tol = MACHEP * 100
    combined = array(list(x) + list(y))
    observed_stat, obs_p = mw_test(x, y)
    total_obs = len(combined)
    num_x = len(x)
    num_greater = 0
    for sampled_x, sampled_y in _get_bootstrap_sample(x, y, num_reps):
        sample_stat, sample_p = mw_test(sampled_x, sampled_y)
        if sample_stat >= (observed_stat - tol):
            num_greater += 1
    return observed_stat, num_greater / num_reps

def permute_2d(m, p):
    """Performs 2D permutation of matrix m according to p."""
    return m[p][:, p]
    #unused below
    m_t = transpose(m)
    r_t = take(m_t, p, axis=0)
    return take(transpose(r_t), p, axis=0)

def mantel(m1, m2, n):
    """Compares two distance matrices. Reports P-value for correlation.
    
    The p-value is based on a two-sided test.

    WARNING: The two distance matrices must be symmetric, hollow distance
    matrices, as only the lower triangle (excluding the diagonal) will be used
    in the calculations (matching R's vegan::mantel function).

    This function is retained for backwards-compatibility. Please use
    mantel_test() for more control over how the test is performed.
    """
    return mantel_test(m1, m2, n)[0]

def mantel_test(m1, m2, n, alt="two sided",
                suppress_symmetry_and_hollowness_check=False):
    """Runs a Mantel test on two distance matrices.

    Returns the p-value, Mantel correlation statistic, and a list of Mantel
    correlation statistics for each permutation test.

    WARNING: The two distance matrices must be symmetric, hollow distance
    matrices, as only the lower triangle (excluding the diagonal) will be used
    in the calculations (matching R's vegan::mantel function).

    Arguments:
        m1  - the first distance matrix to use in the test (should be a numpy
            array or convertible to a numpy array)
        m2  - the second distance matrix to use in the test (should be a numpy
            array or convertible to a numpy array)
        n   - the number of permutations to test when calculating the p-value
        alt - the type of alternative hypothesis to test (can be either
            'two sided' for a two-sided test, 'greater' or 'less' for one-sided
            tests)
        suppress_symmetry_and_hollowness_check - by default, the input distance
            matrices will be checked for symmetry and hollowness. It is
            recommended to leave this check in place for safety, as the check
            is fairly fast. However, if you *know* you have symmetric and
            hollow distance matrices, you can disable this check for small
            performance gains on extremely large distance matrices
    """
    # Perform some sanity checks on our input.
    if alt not in ("two sided", "greater", "less"):
        raise ValueError("Unrecognized alternative hypothesis. Must be either "
                         "'two sided', 'greater', or 'less'.")
    m1, m2 = asarray(m1), asarray(m2)
    if m1.shape != m2.shape:
        raise ValueError("Both distance matrices must be the same size.")
    if n < 1:
        raise ValueError("The number of permutations must be greater than or "
                         "equal to one.")
    if not suppress_symmetry_and_hollowness_check:
        if not (is_symmetric_and_hollow(m1) and is_symmetric_and_hollow(m2)):
            raise ValueError("Both distance matrices must be symmetric and "
                             "hollow.")

    # Get a flattened list of lower-triangular matrix elements (excluding the
    # diagonal) in column-major order. Use these values to calculate the
    # correlation statistic.
    m1_flat, m2_flat = _flatten_lower_triangle(m1), _flatten_lower_triangle(m2)
    orig_stat = pearson(m1_flat, m2_flat)

    # Run our permutation tests so we can calculate a p-value for the test.
    size = len(m1)
    better = 0
    perm_stats = []
    for i in range(n):
        perm = permute_2d(m1, permutation(size))
        perm_flat = _flatten_lower_triangle(perm)
        r = pearson(perm_flat, m2_flat)

        if alt == 'two sided':
            if abs(r) >= abs(orig_stat):
                better += 1
        else:
            if ((alt == 'greater' and r >= orig_stat) or
                (alt == 'less' and r <= orig_stat)):
                better += 1
        perm_stats.append(r)
    return (better + 1) / (n + 1), orig_stat, perm_stats

def is_symmetric_and_hollow(matrix):
    return (matrix.T == matrix).all() and (trace(matrix) == 0)

def _flatten_lower_triangle(matrix):
    """Returns a list containing the flattened lower triangle of the matrix.

    The returned list will contain the elements in column-major order. The
    diagonal will be excluded.

    Arguments:
        matrix - numpy array containing the matrix data
    """
    matrix = asarray(matrix)
    flattened = []
    for col_num in range(matrix.shape[1]):
        for row_num in range(matrix.shape[0]):
            if col_num < row_num:
                    flattened.append(matrix[row_num][col_num])
    return flattened

def kendall_correlation(x, y, alt="two sided", exact=None, warn=True):
    """returns the statistic (tau) and probability from Kendall's non-parametric
    test of association that tau==0. Uses the large sample approximation when
    len(x) >= 50 or when there are ties, otherwise it computes the probability
    exactly.
    
    Based on the algorithm implemented in R v2.5
    
    Arguments:
        - alt: the alternate hypothesis (greater, less, two sided)
        - exact: when False, forces use of the large sample approximation
          (normal distribution). Not allowed for len(x) >= 50.
        - warn: whether to warn about tied values
    """
    
    assert len(x) == len(y), "data (x, y) not of same length"
    assert len(x) > 2, "not enough observations"
    
    # possible alternate hypotheses arguments
    lo = ["less", "lo", "lower", "l", "lt"]
    hi = ["greater", "hi", "high", "h", "g", "gt"]
    two = ["two sided", "2", 2, "two tailed", "two", "two.sided", "ts"]
    
    ties = False
    num = len(x)
    ties = len(set(x)) != num or len(set(y)) != num
    if ties and warn:
        warnings.warn("Tied values, using normal approximation")
    
    if not ties and num < 50:
        exact = True
    
    if num < 50 and not ties and exact:
        combs = int(num * (num-1) / 2)
        working = []
        for i in range(combs):
            row = [-1 for j in range(combs)]
            working.append(row)
        
        tau = kendalls_tau(x, y, False)
        q = round((tau+1)*num*(num-1) / 4)
        if alt in two:
            if q > num * (num - 1) / 4:
                p = 1 - pkendall(q-1, num, Gamma(num+1), working)
            else:
                p = pkendall(q, num, Gamma(num+1), working)
            p = min(2*p, 1)
        elif alt in hi:
            p = 1 - pkendall(q-1, num, Gamma(num+1), working)
        elif alt in lo:
            p = pkendall(q, num, Gamma(num+1), working)
    else:
        tau, p = kendalls_tau(x, y, True)
        if alt in hi:
            p /= 2
        elif alt in lo:
            p = 1 - p/2
    return tau, p

## Start functions for distance_matrix_permutation_test

def distance_matrix_permutation_test(matrix, cells, cells2=None,\
        f=t_two_sample, tails=None, n=1000, return_scores=False,\
        is_symmetric=True):
    """performs a monte carlo permutation test to determine if the 
    values denoted in cells are significantly different than the rest
    of the values in the matrix

    matrix: a numpy array
    cells: a list of indices of special cells to compare to the rest of the 
        matrix
    cells2: an optional list of indices to compare cells to. If set to None
        (default), compares cells to the rest of the matrix
    f: the statistical test used. Should take a "tails" parameter as input
    tails: can be None(default), 'high', or 'low'. Input into f.
    n: the number of replicates in the Monte Carlo simulations
    is_symmetric: corrects if the matrix is symmetric. Need to only look at
        one half otherwise the degrees of freedom value will be incorrect.
    """
    #if matrix is symmetric convert all indices to lower trangular
    if is_symmetric:
        cells = get_ltm_cells(cells)
        if cells2:
            cells2 = get_ltm_cells(cells2)
    # pull out the special values
    special_values, other_values = \
        get_values_from_matrix(matrix, cells, cells2, is_symmetric)
    # calc the stat and parameteric p-value for real data
    stat, p = f(special_values, other_values, tails)
    #calc for randomized matrices
    count_more_extreme = 0
    stats = []
    indices = range(len(matrix))
    for k in range(n):
        # shuffle the order of indices, and use those to permute the matrix
        permuted_matrix = permute_2d(matrix,permutation(indices))
        special_values, other_values = \
            get_values_from_matrix(permuted_matrix, cells,\
            cells2, is_symmetric)
        # calc the stat and p for a random subset (we don't do anything 
        # with these p-values, we only use the current_stat value)
        current_stat, current_p = f(special_values, other_values, tails)
        stats.append(current_stat)
        if tails == None:
            if abs(current_stat) > abs(stat): count_more_extreme += 1
        elif tails == 'low':
            if current_stat < stat: count_more_extreme += 1
        elif tails == 'high':
            if current_stat > stat: count_more_extreme += 1

    # pack up the parametric stat, parametric p, and empirical p; calc the
    # the latter in the process
    result = [stat, p, count_more_extreme/n]
    # append the scores of the n tests if requested
    if return_scores: result.append(stats)
    return tuple(result)

def get_values_from_matrix(matrix, cells, cells2=None, is_symmetric=True):
    """get values from matrix positions in cells and cells2

        matrix: the numpy array from which values should be taken
        cells: indices of first set of requested values
        cells2: indices of second set of requested values or None
         if they should be randomly selected
        is_symmetric: True if matrix is symmetric 

    """

    # pull cells values
    cells_values = [matrix[i] for i in cells]
    # pull cells2 values
    if cells2:
        cells2_values = [matrix[i] for i in cells2]
    # or generate the indices and grab them if they
    # weren't passed in
    else:
        cells2_values = []
        for i, val_i in enumerate(matrix):
            for j, val in enumerate(val_i):
                if is_symmetric:
                    if (i,j) not in cells and i > j:
                        cells2_values.append(val)
                else:
                    if (i,j) not in cells:
                        cells2_values.append(val)
    return cells_values, cells2_values

def get_ltm_cells(cells):
    """converts matrix indices so all are below the diagonal

        cells: list of indices into a 2D integer-indexable object
         (typically a list or lists of array of arrays)
    
    """
    new_cells = []
    for cell in cells:
        if cell[0] < cell[1]:
            new_cells.append((cell[1], cell[0]))
        elif cell[0] > cell[1]:
            new_cells.append(cell)
    #remove duplicates
    new_cells = set(new_cells)
    return list(new_cells)

## End functions for distance_matrix_permutation_test
