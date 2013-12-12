#!/usr/bin/env python
from __future__ import division
from cogent.maths.stats.special import lgam
from cogent.maths.optimisers import minimise
from math import ceil, e
from numpy import array, zeros, concatenate, arange, log, sqrt, exp, asarray
from numpy.random import gamma, shuffle
from cogent.maths.scipy_optimize import fmin_powell
import cogent.maths.stats.rarefaction as rarefaction
from functools import reduce

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2012, The Cogent Project"
__credits__ = ["Rob Knight", "Justin Kuczynski, William Van Treuren",
               "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.5.3-dev"
__maintainer__ = "Rob Knight"
__email__ = "rob@spot.colorado.edu"
__status__ = "Production"


def expand_counts(counts):
    """Converts vector of counts at each index to vector of indices."""
    result = []
    for i, c in enumerate(counts):
        result.append(zeros(c, int) + i)
    return concatenate(result)


def counts(indices, result=None):
    """Converts vector of indices to counts of each index.

    WARNING: does not check that 'result' array is big enough to store new
    counts, suggest preallocating based on whole dataset if doing cumulative
    analysis.
    """
    if result is None:
        max_val = indices.max()
        result = zeros(max_val + 1)
    for i in indices:
        result[i] += 1
    return result


def singles(counts):
    """Returns count of single occurrences."""
    return (counts == 1).sum()


def doubles(counts):
    """Returns count of double occurrences."""
    return (counts == 2).sum()


def observed_species(counts):
    """Calculates number of distinct species."""
    return (counts != 0).sum()


def osd(counts):
    """Returns observed, singles and doubles from counts.

    Handy for diversity calculations."""
    return (counts != 0).sum(), (counts == 1).sum(), (counts == 2).sum()


def margalef(counts):
    """Margalef's index, assumes log accumulation. Magurran 2004, p 77."""
    return (observed_species(counts) - 1) / log(counts.sum())


def menhinick(counts):
    """Menhinick's index, assumes sqrt accumulation. Magurran 2004, p 77."""
    return observed_species(counts) / sqrt(counts.sum())


def dominance(counts):
    """Dominance = 1 - Simpson's index, sum of squares of probs.

    Briefly, gives probability that two species sampled are the same."""
    freqs = counts / float(counts.sum())
    return (freqs * freqs).sum()


def simpson(counts):
    """Simpson's index = 1-dominance."""
    return 1 - dominance(counts)


def reciprocal_simpson(counts):
    """1/Simpson's index"""
    return 1.0 / simpson(counts)


def simpson_reciprocal(counts):
    """1/D (1/dominance)"""
    return 1.0 / dominance(counts)


def shannon(counts, base=2):
    """Returns Shannon entropy of counts, default in bits."""
    freqs = counts / float(counts.sum())
    nonzero_freqs = freqs[freqs.nonzero()]
    return -sum(nonzero_freqs * log(nonzero_freqs)) / log(base)


def equitability(counts, base=2):
    """Returns Shannon index corrected for # species, pure evenness."""
    return shannon(counts, base) / (log((counts != 0).sum()) / log(base))


def berger_parker_d(counts):
    """Fraction of the sample that belongs to the most abundant species.

    Berger & Parker 1970, by way of SDR-IV online help.
    """
    return counts.max() / float(counts.sum())


def mcintosh_d(counts):
    """McIntosh index of alpha diversity (McIntosh 1967, by way of SDR-IV)."""
    u = sqrt((counts * counts).sum())
    n = counts.sum()
    return (n - u) / (n - sqrt(n))


def brillouin_d(counts):
    """Brilloun index of alpha diversity: Pielou 1975, by way of SDR-IV."""
    nz = counts[counts.nonzero()]
    n = nz.sum()
    return (lgam(n + 1) - array(map(lgam, nz + 1)).sum()) / n


def kempton_taylor_q(counts, lower_quantile=.25, upper_quantile=.75):
    """Kempton-Taylor (1976) q index of alpha diversity, by way of SDR-IV.

    Estimates the slope of the cumulative abundance curve in the interquantile
    range. By default, uses lower and upper quartiles, rounding inwards.

    Note: this differs slightly from the results given in Magurran 1998.
    Specifically, we have 14 in the numerator rather than 15. Magurran
    recommends counting half of the species with the same # counts as the
    point where the UQ falls and the point where the LQ falls, but the
    justification for this is unclear (e.g. if there were a very large #
    species that just overlapped one of the quantiles, the results would
    be considerably off). Leaving the calculation as-is for now, but consider
    changing.
    """
    n = len(counts)
    lower = int(ceil(n * lower_quantile))
    upper = int(n * upper_quantile)
    sorted = sorted(counts.copy())
    return (upper - lower) / log(sorted[upper] / sorted[lower])


def strong(counts):
    """Strong's 2002 dominance index, by way of SDR-IV."""
    cc = sorted(counts.copy())
    cc = cc[::-1]
    sorted_sum = cc.cumsum()
    n = counts.sum()
    s = (counts != 0).sum()
    i = arange(1, len(counts) + 1)
    return (sorted_sum / float(n) - (i / float(s))).max()


def fisher_alpha(counts, bounds=(1e-3, 1e12)):
    """Fisher's alpha: S = alpha ln(1+N/alpha) where S=species, N=individuals

    bounds are guess for Powell optimizer bracket search.

    WARNING: may need to adjust bounds for some datasets.
    """
    n = counts.sum()
    s = (counts != 0).sum()

    def f(alpha):
        return (alpha * log(1 + (n / alpha)) - s) ** 2

    alpha = minimise(f, 1.0, bounds, local=True)

    if f(alpha) > 1.0:
        raise RuntimeError("optimizer failed to converge (error > 1.0)," +
                           " so no fisher alpha returned")
    return alpha


def mcintosh_e(counts):
    """McIntosh's evenness measure: Heip & Engels 1974 p 560 (wrong in SDR-IV)."""
    numerator = sqrt((counts * counts).sum())
    n = counts.sum()
    s = (counts != 0).sum()
    denominator = sqrt((n - s + 1) ** 2 + s - 1)
    return numerator / denominator


def heip_e(counts):
    """Heip's evenness measure: Heip & Engels 1974."""
    return exp(shannon(counts, base=e) - 1) / ((counts != 0).sum() - 1)


def simpson_e(counts):
    """Simpson's evenness.
    Previous formula updated based on discussion at
    https://groups.google.com/forum/#!searchin/qiime-forum/simpsons%2420e/qiime-forum/F7GvTT2Elgw/rU8fnNfLoNMJ"""
    return (1. / dominance(counts)) / (counts != 0).sum()


def robbins(counts):
    """Robbins 1968 estimator for Pr(unobserved) at n trials.

    H. E. Robbins (1968, Ann. of Stats. Vol 36, pp. 256-257)
    probability_of_unobserved_colors = S/(n+1),
    (where s = singletons).

    Note that this is the estimate for (n-1) counts, i.e. x axis is off by 1.
    """
    return float(singles(counts)) / counts.sum()


def robbins_confidence(counts, alpha=0.05):
    """Robbins 1968 confidence interval for counts given alpha.

    Note: according to Manuel's equation, if alpha=0.05, we get a 95% CI.
    """
    s = singles(counts)
    n = counts.sum()
    k = sqrt((n + 1) / alpha)
    return (s - k) / (n + 1), (s + k) / (n + 1)


# TODO: SDR-IV also implements NHC, Carmago, Smith & Wilson, Gini indices.
# possibly worth adding.
# http://www.pisces-conservation.com/sdrhelp/index.html

def michaelis_menten_fit(counts, num_repeats=1, params_guess=None,
                         return_b=False):
    """Michaelis-Menten fit to rarefaction curve of observed species

    Note: there is some controversy about how to do the fitting. The ML model
    givem by Raaijmakers 1987 is based on the assumption that error is roughly
    proportional to magnitude of observation, reasonable for enzyme kinetics
    but not reasonable for rarefaction data. Here we just do a nonlinear
    curve fit for the parameters using least-squares.


    S = Smax*n/(B + n) . n: number of individuals, S: # of species
    returns Smax

    inputs:
    num_repeats: will perform rarefaction (subsampling without replacement)
    this many times at each value of n
    params_guess: intial guess of Smax, B (None => default)
    return_b: if True will return the estimate for Smax, B. Default is just Smax

    the fit is made to datapoints where n = 1,2,...counts.sum(),
    S = species represented in random sample of n individuals

    """
    counts = asarray(counts)
    if params_guess is None:
        params_guess = array([100, 500])

    # observed # of species vs # of individuals sampled, S vs n
    xvals = arange(1, counts.sum() + 1)
    ymtx = []
    for i in range(num_repeats):
        ymtx.append(array([observed_species(rarefaction.subsample(counts, n))
                           for n in xvals]))
    ymtx = asarray(ymtx)
    yvals = ymtx.mean(0)

    # fit to obs_sp = max_sp * num_idiv / (num_indiv + B)
    # return max_sp
    def fitfn(p, n):  # works with vectors of n, returns vector of S
        return p[0] * n / (p[1] + n)

    def errfn(p, n, y):  # vectors of actual vals y and number of individuals n
        return ((fitfn(p, n) - y) ** 2).sum()

    p1 = fmin_powell(errfn, params_guess, args=(xvals, yvals), disp=0)
    if return_b:
        return p1
    else:
        return p1[0]  # return only S_max, not the K_m (B) param


def chao1_uncorrected(observed, singles, doubles):
    """Calculates chao1 given counts. Eq. 1 in EstimateS manual.

    Formula: chao1 = S_obs + N_1^2/(2*N_2) where N_1 and N_2 are
    count of singletons and doubletons respectively.

    Note: this is the original formula from Chao 1984, not bias-corrected,
    and is Equation 1 in the EstimateS manual.
    """
    return observed + singles ** 2 / float(doubles * 2)


def chao1_bias_corrected(observed, singles, doubles):
    """Calculates bias-corrected chao1 given counts: Eq. 2 in EstimateS manual.

    Formula: chao1 = S_obs + N_1(N_1-1)/(2*(N_2+1)) where N_1 and N_2 are
    count of singletons and doubletons respectively.

    Note: this is the bias-corrected formulat from Chao 1987, Eq. 2 in the
    EstimateS manual.
    """
    return observed + singles * (singles - 1) / (2.0 * (doubles + 1))


def chao1(counts, bias_corrected=True):
    """Calculates chao1 according to table in EstimateS manual.

    Specifically, uses bias-corrected version unless bias_corrected is set
    to False _and_ there are both singletons and doubletons."""
    o, s, d = osd(counts)
    if not bias_corrected:
        if s:
            if d:
                return chao1_uncorrected(o, s, d)
    return chao1_bias_corrected(o, s, d)


def chao1_var_uncorrected(singles, doubles):
    """Calculates chao1, uncorrected.

    From EstimateS manual, equation 5.
    """
    r = float(singles) / doubles
    return doubles * (.5 * r ** 2 + r ** 3 + .24 * r ** 4)


def chao1_var_bias_corrected(singles, doubles):
    """Calculates chao1 variance, bias-corrected.

    From EstimateS manual, equation 6.
    """
    s, d = float(singles), float(doubles)
    return s * (s - 1) / (2 * (d + 1)) + (s * (2 * s - 1) ** 2) / (4 * (d + 1) ** 2) + \
        (s ** 2 * d * (s - 1) ** 2) / (4 * (d + 1) ** 4)


def chao1_var_no_doubletons(singles, chao1):
    """Calculates chao1 variance in absence of doubletons.

    From EstimateS manual, equation 7.

    chao1 is the estimate of the mean of Chao1 from the same dataset.
    """
    s = float(singles)
    return s * (s - 1) / 2 + s * (2 * s - 1) ** 2 / 4 - s ** 4 / (4 * chao1)


def chao1_var_no_singletons(n, observed):
    """Calculates chao1 variance in absence of singletons. n = # individuals.

    From EstimateS manual, equation 8.
    """
    o = float(observed)
    return o * exp(-n / o) * (1 - exp(-n / o))


def chao1_var(counts, bias_corrected=True):
    """Calculates chao1 variance using decision rules in EstimateS."""
    o, s, d = osd(counts)
    if not d:
        c = chao1(counts, bias_corrected)
        return chao1_var_no_doubletons(s, c)
    if not s:
        n = counts.sum()
        return chao1_var_no_singletons(n, o)
    if bias_corrected:
        return chao1_var_bias_corrected(s, d)
    else:
        return chao1_var_uncorrected(s, d)


def chao_confidence_with_singletons(chao, observed, var_chao, zscore=1.96):
    """Calculates confidence bounds for chao1 or chao2.

    Uses Eq. 13 of EstimateS manual.

    zscore = score to use for confidence, default = 1.96 for 95% confidence.
    """
    T = float(chao - observed)
    # if no diff betweeh chao and observed, CI is just point estimate of
    # observed
    if T == 0:
        return observed, observed
    K = exp(abs(zscore) * sqrt(log(1 + (var_chao / T ** 2))))
    return observed + T / K, observed + T * K


def chao_confidence_no_singletons(n, observed, zscore=1.96):
    """Calculates confidence bounds for chao1/chao2 in absence of singletons.

    Uses Eq. 14 of EstimateS manual.

    n = number of individuals, observed = number of species.
    """
    s = float(observed)
    P = exp(-n / s)
    return max(s, s / (1 - P) - zscore * sqrt((s * P / (1 - P)))), \
        s / (1 - P) + zscore * sqrt(s * P / (1 - P))


def chao1_confidence(counts, bias_corrected=True, zscore=1.96):
    """Returns chao1 confidence (lower, upper) from counts."""
    o, s, d = osd(counts)
    if s:
        chao = chao1(counts, bias_corrected)
        chaovar = chao1_var(counts, bias_corrected)
        return chao_confidence_with_singletons(chao, o, chaovar, zscore)
    else:
        n = counts.sum()
        return chao_confidence_no_singletons(n, o, zscore)


def chao1_lower(*args, **kwargs):
    """Convenience wrapper for chao1_confidence to fit overall diversity API."""
    return chao1_confidence(*args, **kwargs)[0]


def chao1_upper(*args, **kwargs):
    """Convenience wrapper for chao1_confidence to fit overall diversity API."""
    return chao1_confidence(*args, **kwargs)[1]


def ACE(count, rare_threshold=10):
    """Implements the ACE metric from EstimateS. Based on the equations
    given under ACE:Abundance-based Coverage Estimator.

    count = an OTU by sample vector
    rare_threshold = threshold at which a species containing as many or
    fewer individuals will be considered rare.

    IMPORTANT NOTES:

    Raises a value error if every rare species is a singleton.

    if no rare species exist, just returns the number of abundant species

    rare_threshold default value is 10. Based on Chao 2000 in Statistica
    Sinica pg. 229 citing empirical observations by Chao, Ma, Yang 1993.

    If the count vector contains 0's, indicating species which are known
    to exist in the environment but did not appear in the sample, they
    will be ignored for the purpose of calculating s_rare."""

    def frequency_counter(count):
        """Creates a frequency count array to beused by every other function."""
        return counts(count)

    def species_rare(freq_counts, rare_threshold):
        """freq_counts number of rare species. Default value of rare is 10 or
        fewer individuals. Based on Chao 2000 in Statistica Sinica pg. 229
        citing empirical observations by Chao, Ma and Yang in 1993."""
        return freq_counts[1:rare_threshold + 1].sum()

    def species_abundant(freq_counts, rare_threshold):
        """freq_counts number of abundant species. Default value of abundant is
        greater than 10 individuals. Based on Chao 2000 in Statistica Sinica
        pg.229  citing observations by Chao, Ma and Yang in 1993."""
        return freq_counts[rare_threshold + 1:].sum()

    def number_rare(freq_counts, gamma=False):
        """Number of individuals in rare species. gamma=True generates the
        n_rare used for the variation coefficient."""

        n_rare = 0
        if gamma == True:
            for i, j in enumerate(freq_counts[:rare_threshold + 1]):
                n_rare = n_rare + (i * j) * (i - 1)
            return n_rare

        for i, j in enumerate(freq_counts[:rare_threshold + 1]):
            n_rare = n_rare + (i * j)
        return n_rare

    # calculations begin

    freq_counts = frequency_counter(count)

    if freq_counts[1:rare_threshold].sum() == 0:
        return species_abundant(freq_counts, rare_threshold)

    if freq_counts[1] == freq_counts[1:rare_threshold].sum():
        raise ValueError("only rare species are singletons, ACE " +
                         "metric is undefined. EstimateS suggests using bias corrected Chao1")

    s_abun = species_abundant(freq_counts, rare_threshold)

    s_rare = species_rare(freq_counts, rare_threshold)

    n_rare = number_rare(freq_counts)

    c_ace = 1 - (freq_counts[1]).sum() / float(n_rare)

    top = s_rare * number_rare(freq_counts, gamma=True)
    bottom = c_ace * n_rare * (n_rare - 1.0)

    gamma_ace = (top / bottom) - 1.0

    if 0 > gamma_ace:
        gamma_ace = 0

    return s_abun + (s_rare / c_ace) + ((freq_counts[1] / c_ace) * gamma_ace)


def diversity(indices, f=chao1, step=1, start=None, verbose=False):
    """Calculates diversity index (default:chao1) for each window of size step.

    indices: vector of indices of species
    f: f(counts) -> diversity measure (default: chao1)
    start: first index to sum up to (default: step)
    step: step size (default:1)

    Note: use rarefaction module if you need more versatility/speed.
    """
    result = []
    if start is None:
        start = step
    freqs = zeros(max(indices) + 1)
    i = 0
    for j in range(start, len(indices) + 1, step):
        freqs = counts(indices[i:j], freqs)
        try:
            curr = f(freqs)
        except (ZeroDivisionError, FloatingPointError, e):
            curr = 0
        if verbose:
            print curr
        result.append(curr)
        i = j
    return array(result)


def lorenz_curve(data):
    """Calculates the Lorenz curve for input data.
    Notes:
     Formula available on wikipedia.
    Inputs:
     data - list or 1d arr, counts/abundances/proportions etc. All entries must
     be non-negative."""
    if any(array(data) < 0):
        raise ValueError('Lorenz curves aren\'t meaningful for non-positive ' +
                         'data.')
    # dont wan't to change input, copy and sort
    sdata = array(sorted((data[:])))
    n = float(len(sdata))
    Sn = sdata.sum()
    # ind+1 because must sum first point, eg. x[:0] = []
    lorenz_points = [((ind + 1) / n, sdata[:ind + 1].sum() / Sn)
                     for ind in range(int(n))]
    return lorenz_points


def lorenz_curve_integrator(lc_pts, method):
    """Calculates the area under a lorenz curve.
    Notes:
     Could be utilized for integrating other simple, non-pathological
     'functions' where width of the trapezoids is constant.
     Two methods are available.
     1. Trapezoids, connecting the lc_pts by linear segments between them.
        Basically assumes that given sampling is accurate and that more features
        of given data would fall on linear gradients between the values of this
        data. formula is: dx[(h_0+h_n)/2 + sum(i=1,i=n-1,h_i)]
     2. Rectangles, connecting lc_pts by lines parallel to x axis. This is the
        correct method in my opinion though trapezoids might be desirable in
        some circumstances. forumla is : dx(sum(i=1,i=n,h_i))
    Inputs:
     lc_pts - list of tuples, output of lorenz_curve.
     method - str, either 'rectangles' or 'trapezoids'
    """
    if method is 'trapezoids':
        dx = 1. / len(lc_pts)  # each point differs by 1/n
        h_0 = 0.0  # 0 percent of the population has zero percent of the goods
        h_n = lc_pts[-1][1]
        sum_hs = sum([pt[1] for pt in lc_pts[:-1]])  # the 0th entry is at x=
        # 1/n
        return dx * ((h_0 + h_n) / 2. + sum_hs)
    elif method is 'rectangles':
        dx = 1. / len(lc_pts)  # each point differs by 1/n
        return dx * sum([pt[1] for pt in lc_pts])
    else:
        return ValueError(('Method \'%s\' not implemeneted.' +
                           ' Available methods \'rectangles\', \'trapezoids\'.') % method)


def gini_index(data, method='rectangles'):
    """Calculates the gini index of data.
    Notes:
     formula is G = A/(A+B) where A is the area between y=x and the Lorenz curve
     and B is the area under the Lorenz curve. Simplifies to 1-2B since A+B=.5
     Formula available on wikipedia.
    Inputs:
     data - list or 1d arr, counts/abundances/proportions etc. All entries must
     be non-negative.
     method - str, either 'rectangles' or 'trapezoids'. see
     lorenz_curve_integrator for details.
    """
    lorenz_points = lorenz_curve(data)
    B = lorenz_curve_integrator(lorenz_points, method)
    return 1 - 2 * B


def enspie(counts):
    '''Calculates ENS_pie alpha diversity measure

    ENS_pie = 1/sum(pi^2) with the sum occuring over all S species in the pool.
    pi is the proportion of the entire community that species i represents.

    The forumla is taken from "Scale-dependent effect sizes of ecological
    drivers on biodiversity: why standardised sampling is not enough" Chase and
    Knight. Ecology Letters, Volume 16, Issue Supplement s1, pgs 17-26 May
    2013. http://onlinelibrary.wiley.com/doi/10.1111/ele.12112/full

    It appears that ENS_pie is just 1/dominance(counts).
    '''
    return 1. / dominance(counts)


def lladser_point_estimates(sample, r=10):
    """Series of point estimates of the conditional uncovered probability

    sample: series of random observations
    r: Number of new colors that are required for the next prediction

    This is the point estimator described in Theorem 2 (i):
    Lladser, Gouet, and Reeder, "Extrapolation of Urn Models via Poissonization:
    Accurate Measurements of the Microbial Unknown" PLoS 2011.
    Returns: Each new color yields 3-tuple:
         - point estimate
         - position in sample of prediction
         - random variable from poisson process (mostly to make testing easier)
    """

    if(r <= 2):
        raise ValueError("r must be >=3")
    for count, seen, cost, i in get_interval_for_r_new_species(sample, r):
        t = gamma(count, 1)
        point_est = (r - 1) / t
        yield(point_est, i, t)


def get_interval_for_r_new_species(seq, r):
    """For seq of samples compute interval between r new species.

    seq: series of observations (the actual sample, not the frequencies)
    r: Number of new colors to that need to be observed for a new interval

    Imagine an urn with colored balls. Given a drawing of balls from the urn,
    compute how many balls need to be looked at to discover r new colors.
    Colors can be repeated.

    Returns: for each new color seen for the first time, yield:
             - length of interval, i.e. number of observations looked at
             - the set of seen colors
             - position in seq after seeing the last new color (end of interval)
             - position in seq where interval is started
    """

    seen = set()
    seq_len = len(seq)
    for i in range(seq_len):
        curr = seq[i]  # note: first iteration is after looking at first char
        # bail out if there's nothing new
        if curr in seen:
            continue
        else:
            seen.add(curr)

        # otherwise, need to see distance to get k colors
        unseen = 0
        j = i + 1
        while unseen < r and j < seq_len:
            if seq[j] not in seen:
                unseen += 1
            j += 1  # note: increments after termination condition

        count = j - i - 1  # the interval to see r new colors
        cost = j  # the position in seq after seeing r new ones
        if (not count) or (unseen < r):  # bail out if not enough unseen
            raise StopIteration
        yield count, set(seen), cost, i


def esty_ci(counts, **args):
    """Esty's CI for (1-m).

    counts: Vector of counts (NOT the sample)

    Esty's CI is defined in
    Esty WW (1983) A Normal limit law for a nonparametric estimator of the
    coverage of a random sample. Ann Statist 11: 905-912.

    n1 / n  +/- z * square-root(W);

    where
    n1 = number of species observed once in n samples;
    n = sample size;
    z = a constant that depends on the targeted confidence and based on
        the Normal distribution. For a 95% CI, z=1.959963985;
    n2 = number of species observed twice in n samples;
    W = [ n1*(n - n1)  +  2*n*n2 ] / (n**3).

    Note: for other confidence levels we first need the appropriate z,
          Not yet hooked up to CLI.

    Returns: (upper bound, lower bound)
    """

    n1 = singles(counts)
    n2 = doubles(counts)
    n = counts.sum()
    z = 1.959963985
    W = (n1 * (n - n1) + 2 * n * n2) / (n ** 3)

    return n1 / n + z * sqrt(W), n1 / n - z * sqrt(W)


def starr_est(sample, m=1):
    """Series of Starr estimates for the unseen mass.

    sample: the series of observations
    m: Starr's lookahead

    From: Starr N (1979) Linear estimation of the probability of discovering a new species.
          Ann Stat 7: 644-652.

    Note: No test code, thus not hooked up to CLI yet.
    Note: for m=1 starr equals robbins

    TODO: If we ever only are interested in the last estimate,
          we can speed up the code considerably

    Returns: list of estimates for each position of sample
    """
    estimates = []
    counts = defaultdict(int)
    # precompute since it is constant for chosen m
    fact_m_1 = factorial(m - 1)
    for n in range(1, len(sample)):
        sum = 0
        counts[sample[n - 1]] += 1
        count = array(counts.values())
        filtered_count = array([x for x in count if x <= (m + 1)])
        for k in range(1, m + 1):
            k_tons = (filtered_count == k).sum()

            if k_tons != 0:
                # computing factorials for large number is inefficient
                # according to manuel (and some tests) the above formula is
                # equal to:
                # increment = binomial_coeff(m-1, k-1)/binomial_coeff(n+m,k)
                increment = k * fact_m_1 / factorial(m - k)  \
                              * (1 / reduce(mul, range(n + m - k + 1, n + m + 1)))
                sum += increment * k_tons

        estimates.append(sum)
    return estimates


def lladser_ci_series(seq, r, alpha=0.95, f=10, ci_type='ULCL'):
    """Constructs r-color confidence intervals for the uncovered conditional prob.

    seq: Input is a sequence of colors (the actual sample, not the counts)
    r  : Number of new colors that are required for the next prediction
    alpha: desired confidence level
    f: ratio between upper and lower bound
    ci_type: type of confidence interval. One of:
             ULCL: upper and lower bounds with conservative lower bound
             ULCU: upper and lower woth conservative upper bound
             U: Upper bound only, lower bound fixed to 0
             L: Lower bound only, upper bound fixed to 1

    Returns: One CI prediction for each new color that is detected and where.
    Lladser, Gouet, and Reeder, "Extrapolation of Urn Models via Poissonization:
    Accurate Measurements of the Microbial Unknown" PLoS 2011.
    """

    for count, seen, cost, i in get_interval_for_r_new_species(seq, r):
        t = gamma(count, 1)
        yield lladser_ci_from_r(r, t, alpha, f, ci_type)


def lladser_ci_from_r(r, t, alpha=0.95, f=10, ci_type='ULCL'):
    """Constructs r-color confidence intervals for the uncovered conditional prob.

    r: Number of new colors that are required for the next prediction
    t: A value from the gamma distribution gamma (count,1)
    alpha: desired confidence level
    f: ratio between upper and lower bound
    ci_type: type of confidence interval. One of:
             ULCL: upper and lower bounds with conservative lower bound
             ULCU: upper and lower woth conservative upper bound
             U: Upper bound only, lower bound fixed to 0
             L: Lower bound only, upper bound fixed to 1

    This is the formula that is described in Theorem 2 iii
    Lladser, Gouet, and Reeder, "Extrapolation of Urn Models via Poissonization:
    Accurate Measurements of the Microbial Unknown" PLoS 2011.

    Returns: A confidence interval that contains the true conditional
             uncovered probability with a probability of 100% * alpha.
    """

    if ci_type == 'U':
        upper_bound = upper_confidence_bound(r, alpha) / t
        return (0, upper_bound)
    elif ci_type == 'L':
        lower_bound = lower_confidence_bound(r, alpha) / t
        return (lower_bound, 1)

    bound_params = ul_confidence_bounds(f, r, alpha)
    if ci_type == 'ULCL':
        bound_param = bound_params[0]
    elif ci_type == 'ULCU':
        bound_param = bound_params[1]
    else:
        raise ValueError("Unknown ci_type: %s" % (ci_type))

    upper_bound = bound_param * f / t
    lower_bound = bound_param / t

    # make sure upper bound is at most 1
    if (upper_bound > 1):
        upper_bound = 1

    return lower_bound, upper_bound


def lower_confidence_bound(r, alpha):
    """Get constant for confidence interval with upper bound fixed to 1

    r: number of new colors
    alpha: Confidence interval (for 95% conf use 0.95)

    Compute constant b according to Theorem 2 iii with b=1
    aka c_3 from Table 3
    Lladser, Gouet, and Reeder, "Extrapolation of Urn Models via Poissonization:
    Accurate Measurements of the Microbial Unknown" PLoS 2011.
    Returns: Constant c such that the confidence interval is [c/T_r, 1]
    """

    alpha = round(alpha, 2)
    if not (alpha == 0.95):
        raise ValueError("alpha must be 0.95")
    data = {1: 0.051293294,
            2: 0.355361510,
            3: 0.817691447,
            4: 1.366318397,
            5: 1.970149568,
            6: 2.613014744,
            7: 3.285315692,
            8: 3.980822786,
            9: 4.695227540,
            10: 5.425405697,
            11: 6.169007289,
            12: 6.924212514,
            13: 7.689578292,
            14: 8.463937522,
            15: 9.246330491,
            16: 10.03595673,
            17: 10.83214036,
            18: 11.63430451,
            19: 12.44195219,
            20: 13.25465160,
            21: 14.07202475,
            22: 14.89373854,
            23: 15.71949763,
            24: 16.54903871,
            25: 17.38212584
            }
    try:
        return (data[r])
    except KeyError:
        raise ValueError("r must be between 1,..,25 or 50")


def upper_confidence_bound(r, alpha):
    """Get constant for confidence interval with lower bound fixed to 0

    r: number of new colors
    alpha: Confidence interval (for 95% conf use 0.95)

    Compute constant b according to Theorem 2 iii with a=0
    aka c_0 from Table 3
    Lladser, Gouet, and Reeder, "Extrapolation of Urn Models via Poissonization:
    Accurate Measurements of the Microbial Unknown" PLoS 2011.

    Returns: Constant c such that the confidence interval is [0,c/T_r]
    """

    alpha = round(alpha, 2)
    if not (alpha == 0.95):
        raise ValueError("alpha must be 0.95")

    data = {1: 2.995732274,
            2: 4.743864518,
            3: 6.295793622,
            4: 7.753656528,
            5: 9.153519027,
            6: 10.51303491,
            7: 11.84239565,
            8: 13.14811380,
            9: 14.43464972,
            10: 15.70521642,
            11: 16.96221924,
            12: 18.20751425,
            13: 19.44256933,
            14: 20.66856908,
            15: 21.88648591,
            16: 23.09712976,
            17: 24.30118368,
            18: 25.49923008,
            19: 26.69177031,
            20: 27.87923964,
            21: 29.06201884,
            22: 30.24044329,
            23: 31.41481021,
            24: 32.58538445,
            25: 33.75240327,
            50: 62.17105670,
            }
    try:
        return (data[r])
    except KeyError:
        raise ValueError("r must be between 1,..,25 or 50")


def ul_confidence_bounds(f, r, alpha):
    """returns confidence bounds based on ratio f and alpha

    f: desired ratio of upper to lower bound
    r: number of new colors
    alpha: Confidence interval (for 95% conf use 0.95)

    This function is just a lookup of some precomputed values.

    Returns: Constant c_1,c_2 such that the confidence interval is:
             [c_1/T_r, c_1*f/T_r] for conservative lower bound intervals and
             [c_2/T_r, c_2*f/T_r] for conservative upper bound intervals
    """
    alpha = round(alpha, 2)
    a = None
    b = None

    if (f, r, alpha) in precomputed_table:
        return precomputed_table[(f, r, alpha)]

    # all others combination are only computed for f=10
    # and alpha = 0.90, 0.95 and 0.99
    if f == 10 and r <= 50:
        if alpha in cbs and r < len(cbs[alpha]):
            (a, b) = cbs[alpha][r]

    if a is None or b is None:
        raise ValueError("No constants are precomputed for the combination of "
                         "f:%f, r:%d, and alpha:%.2f" % (f, r, alpha))

    return (a, b)

# Hack in some special values we used for the paper.
# Since Manuel needs to compute those semi-automatically
# using Maple, we pre-calculate only a few common ones

# precomputed table is {(f, r, alpha):(c_1, c_2)}
precomputed_table = {
    (2, 50, 0.95): (31.13026306, 38.94718565),
    (2, 33, 0.95): (22.3203508, 23.4487304),
    (1.5, 100, 0.95): (79.0424349, 83.22790086),
    (1.5, 94, 0.95): (75.9077267, 76.5492088),
    (2.5, 19, 0.95): (11.26109001, 11.96814857),

    # In the next block for each f, we report the smallest possible value
    # of r from table 4 in the paper
    (80, 2, 0.95): (0.0598276655, 0.355361510),
    (48, 2, 0.95): (0.1013728884, 0.355358676),
    (40, 2, 0.95): (0.1231379857, 0.355320458),
    (24, 2, 0.95): (0.226833483, 0.346045204),
    (20, 3, 0.95): (0.320984257, 0.817610455),
    (12, 3, 0.95): (0.590243030, 0.787721610),
    (10, 4, 0.95): (0.806026244, 1.360288674),
    (6, 6, 0.95): (1.8207383, 2.58658608),
    (5, 7, 0.95): (2.48303930, 3.22806682),
    (3, 14, 0.95): (7.17185045, 8.27008349),
    (2.5, 19, 0.95): (11.26109001, 11.96814857),
    (1.5, 94, 0.95): (75.9077267, 76.5492088),
    (1.25, 309, 0.95): (275.661191, 275.949782)
}

######
# Below are the values used for Theorem 3 iii
# Values hand computed by Manuel Lladser using Maple

cb_90 = [
    (None, None),  # 0, makes indexing easier
    (None, None),  # no feasible solution
    (None, None),  # no feasible solution
    (.5635941995, 1.095834700),
    (.6764656264, 1.744588615),
    (.8018565594, 2.432587343),
    (.9282215025, 3.151897973),
    (1.053433716, 3.894766804),
    (1.177158858, 4.656118177),
    (1.299491033, 5.432468058),
    (1.420604842, 6.221304605),  # 10
    (1.540665805, 7.020746595),
    (1.659812701, 7.829342026),
    (1.778158703, 8.645942495),
    (1.895796167, 9.469621185),
    (2.012801198, 10.29961731),
    (2.129237257, 11.13529724),
    (2.245157877, 11.97612664),
    (2.360608695, 12.82164994),
    (2.475628991, 13.67147502),
    (2.590252861, 14.52526147),  # 20
    (2.704510123, 15.38271151),
    (2.818427036, 16.24356290),
    (2.932026869, 17.10758326),
    (3.045330351, 17.97456551),
    (3.158356050, 18.84432420),
    (None, None),  # not computed
    (None, None),
    (None, None),
    (None, None),
    (3.719850286, 23.22944415),  # 30
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (4.828910181, 32.13892224),  # 40
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (5.924900191, 41.17906791)  # 50
]

cb_95 = [
    (None, None),  # 0
    (None, None),
    (None, None),
    (None, None),
    (.8060262438, 1.360288674),  # 4
    (.9240311584, 1.969902537),
    (1.053998892, 2.613007253),
    (1.185086998, 3.285315518),
    (1.315076337, 3.980822783),
    (4.695227540, 4.695227541),
    (1.570546801, 5.425405698),  # 10
    (1.696229569, 6.169007289),
    (1.820753729, 6.924212513),
    (1.944257622, 7.689578291),
    (2.066857113, 8.463937522),
    (2.188648652, 9.246330491),
    (2.309712994, 10.03595673),
    (2.430118373, 10.83214036),
    (2.549923010, 11.63430451),
    (2.669177032, 12.44195219),
    (2.787923964, 13.25465160),  # 20
    (2.906201884, 14.07202475),
    (3.024044329, 14.89373854),
    (3.141481020, 15.71949763),
    (3.258538445, 16.54903871),
    (3.375240327, 17.38212584),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (3.954097220, 21.59397923),  # 30
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (5.093973695, 30.19573919),  # 40
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (6.217105673, 38.96473258)  # 50
]

cb_99 = [
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (1.360316290, 1.768978323),
    (1.470856924, 2.329171347),
    (1.604478487, 2.906049304),
    (1.741759456, 3.507452949),
    (1.878809285, 4.130199076),  # 10
    (2.014632329, 4.771246173),
    (2.149044735, 5.428180734),
    (2.282101533, 6.099073460),
    (2.413917374, 6.782354878),
    (2.544610844, 7.476728267),
    (2.674289153, 8.181107778),
    (2.803045614, 8.894573463),
    (2.930960779, 9.616337916),
    (3.058104355, 10.34572103),
    (3.184536992, 11.08213063),  # 20
    (3.310311816, 11.82504734),
    (3.435475649, 12.57401269),
    (3.560070013, 13.32861956),
    (3.684131925, 14.08850448),
    (3.807694563, 14.85334135),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (4.41897094, 18.7424258),  # 30
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (5.61643962, 26.7700386),  # 40
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (6.79033616, 35.0324474)  # 50
]

cbs = {0.90: cb_90,
       0.95: cb_95,
       0.99: cb_99}


# CLI driver adopted from alph_diversity.py
# Fake common interface by shoving unused args into **args

def starr(counts, r=1, **args):
    """The Starr estimator for the unseen mass

    Shorthand notation to be called as AlphaDiversityCalc
    """
    sample = expand_counts(counts)
    shuffle(sample)

    # misuse r param for starr lookahead
    return starr_est(sample, m=r)[-1]


def lladser_pe(counts, r=10, **args):
    """Single point estimate of the conditional uncovered probability

    This function is just a wrapper around the full point estimator,
    intended to be called fo a single best estimate on a complete sample.
    """
    sample = expand_counts(counts)
    shuffle(sample)
    try:
        pe = list(lladser_point_estimates(sample, r))[-1][0]
    except IndexError:
        pe = 'NaN'
    return pe


def lladser_ci(counts, r, alpha=0.95, f=10, ci_type='ULCL', **args):
    """Single CI of the conditional uncovered probability

    This function is just a wrapper around the full point estimator,
    intended to be called for a single best estimate on a complete sample.
    """

    sample = expand_counts(counts)
    shuffle(sample)
    try:
        pe = list(lladser_ci_series(sample, r))[-1]
    except IndexError:
        pe = ('NaN', 'NaN')
    return pe

# this is important for one of the tests to pass, unknown why
esty_ci.return_names = ('esty_lower_bound', 'esty_upper_bound')
lladser_ci.return_names = ('lladser_lower_bound', 'lladser_upper_bound')


def goods_coverage(counts):
    """Returns Good's Coverage of counts.

    C = 1 - (n1/N)
    n1 = number of OTUs with abundance of 1
    N = number of individuals (sum of abundances for all OTUs)"""
    n1 = (asarray(counts) == 1).sum()
    N = (asarray(counts)).sum()
    return (1 - (n1 / float(N)))
