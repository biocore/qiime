#!/usr/bin/env python
from __future__ import division
from cogent.maths.stats.special import lgam
from cogent.maths.optimisers import minimise
from math import ceil, e
from numpy import array, zeros, concatenate, arange, log, sqrt, exp, asarray
from cogent.maths.scipy_optimize import fmin_powell
import cogent.maths.stats.rarefaction as rarefaction

__author__ = "Rob Knight"
__copyright__ = "Copyright 2007-2011, The Cogent Project"
__credits__ = ["Rob Knight", "Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.6.0dev"
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
        result = zeros(max_val+1)
    for i in indices:
        result[i] += 1
    return result

def singles(counts):
    """Returns count of single occurrences."""
    return (counts==1).sum()

def doubles(counts):
    """Returns count of double occurrences."""
    return (counts==2).sum()

def observed_species(counts):
    """Calculates number of distinct species."""
    return (counts!=0).sum()

def osd(counts):
    """Returns observed, singles and doubles from counts.

    Handy for diversity calculations."""
    return (counts!=0).sum(), (counts==1).sum(), (counts==2).sum()

def margalef(counts):
    """Margalef's index, assumes log accumulation. Magurran 2004, p 77."""
    return (observed_species(counts)-1)/log(counts.sum())

def menhinick(counts):
    """Menhinick's index, assumes sqrt accumulation. Magurran 2004, p 77."""
    return observed_species(counts)/sqrt(counts.sum())

def dominance(counts):
    """Dominance = 1 - Simpson's index, sum of squares of probs.

    Briefly, gives probability that two species sampled are the same."""
    freqs = counts/float(counts.sum())
    return (freqs*freqs).sum()

def simpson(counts):
    """Simpson's index = 1-dominance."""
    return 1 - dominance(counts)

def reciprocal_simpson(counts):
    """1/Simpson's index"""
    return 1.0/simpson(counts)

def simpson_reciprocal(counts):
    """1/D (1/dominance)"""
    return 1.0/dominance(counts)

def shannon(counts, base=2):
    """Returns Shannon entropy of counts, default in bits."""
    freqs = counts/float(counts.sum())
    nonzero_freqs = freqs[freqs.nonzero()]
    return -sum(nonzero_freqs*log(nonzero_freqs))/log(base)

def equitability(counts, base=2):
    """Returns Shannon index corrected for # species, pure evenness."""
    return shannon(counts, base)/(log((counts!=0).sum())/log(base))

def berger_parker_d(counts):
    """Fraction of the sample that belongs to the most abundant species.

    Berger & Parker 1970, by way of SDR-IV online help.
    """
    return counts.max()/float(counts.sum())

def mcintosh_d(counts):
    """McIntosh index of alpha diversity (McIntosh 1967, by way of SDR-IV)."""
    u = sqrt((counts*counts).sum())
    n = counts.sum()
    return (n-u)/(n-sqrt(n))

def brillouin_d(counts):
    """Brilloun index of alpha diversity: Pielou 1975, by way of SDR-IV."""
    nz = counts[counts.nonzero()]
    n = nz.sum()
    return (lgam(n+1) - array(map(lgam, nz+1)).sum())/n

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
    lower = int(ceil(n*lower_quantile))
    upper = int(n*upper_quantile)
    sorted = counts.copy()
    sorted.sort()
    return (upper-lower)/log(sorted[upper]/sorted[lower])

def strong(counts):
    """Strong's 2002 dominance index, by way of SDR-IV."""
    cc = counts.copy()
    cc.sort()
    cc = cc[::-1]
    sorted_sum = cc.cumsum()
    n = counts.sum()
    s = (counts != 0).sum()
    i = arange(1,len(counts)+1)
    return (sorted_sum/float(n) - (i/float(s))).max()

def fisher_alpha(counts, bounds=(1e-3,1e12)):
    """Fisher's alpha: S = alpha ln(1+N/alpha) where S=species, N=individuals
    
    bounds are guess for Powell optimizer bracket search.

    WARNING: may need to adjust bounds for some datasets.
    """
    n = counts.sum()
    s = (counts!=0).sum()
    def f(alpha):
        return (alpha * log(1 + (n/alpha)) - s)**2
    
    alpha = minimise(f, 1.0, bounds, local=True)
    
    if f(alpha) > 1.0:
        raise RuntimeError("optimizer failed to converge (error > 1.0)," +\
            " so no fisher alpha returned")
    return alpha

def mcintosh_e(counts):
    """McIntosh's evenness measure: Heip & Engels 1974 p 560 (wrong in SDR-IV)."""
    numerator = sqrt((counts*counts).sum())
    n = counts.sum()
    s = (counts!=0).sum()
    denominator = sqrt((n-s+1)**2 + s - 1)
    return numerator/denominator

def heip_e(counts):
    """Heip's evenness measure: Heip & Engels 1974."""
    return exp(shannon(counts, base=e)-1)/((counts!=0).sum()-1)

def simpson_e(counts):
    """Simpson's evenness, from SDR-IV."""
    return reciprocal_simpson(counts)/(counts!=0).sum()

def robbins(counts):
    """Robbins 1968 estimator for Pr(unobserved) at n trials.

    H. E. Robbins (1968, Ann. of Stats. Vol 36, pp. 256-257)
    probability_of_unobserved_colors = S/(n+1),
    (where s = singletons).

    Note that this is the estimate for (n-1) counts, i.e. x axis is off by 1.
    """
    return float(singles(counts))/counts.sum()

def robbins_confidence(counts, alpha=0.05):
    """Robbins 1968 confidence interval for counts given alpha.
    
    Note: according to Manuel's equation, if alpha=0.05, we get a 95% CI.
    """
    s = singles(counts)
    n = counts.sum()
    k = sqrt((n+1)/alpha)
    return (s-k)/(n+1), (s+k)/(n+1)


#TODO: SDR-IV also implements NHC, Carmago, Smith & Wilson, Gini indices.
#possibly worth adding.
#http://www.pisces-conservation.com/sdrhelp/index.html

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
        params_guess = array([100,500])

    # observed # of species vs # of individuals sampled, S vs n
    xvals = arange(1,counts.sum()+1)
    ymtx = []
    for i in range(num_repeats):
        ymtx.append( array([observed_species(rarefaction.subsample(counts,n)) \
        for n in xvals]))
    ymtx = asarray(ymtx)
    yvals = ymtx.mean(0)
    
    # fit to obs_sp = max_sp * num_idiv / (num_indiv + B)
    # return max_sp
    def fitfn(p,n): # works with vectors of n, returns vector of S
        return p[0]*n/(p[1] + n)
    
    def errfn(p,n,y): # vectors of actual vals y and number of individuals n
        return ((fitfn(p,n) - y)**2).sum()

    p1 = fmin_powell(errfn, params_guess, args=(xvals,yvals), disp=0)
    if return_b:
        return p1
    else:
        return p1[0] # return only S_max, not the K_m (B) param

def chao1_uncorrected(observed, singles, doubles):
    """Calculates chao1 given counts. Eq. 1 in EstimateS manual.

    Formula: chao1 = S_obs + N_1^2/(2*N_2) where N_1 and N_2 are
    count of singletons and doubletons respectively.

    Note: this is the original formula from Chao 1984, not bias-corrected,
    and is Equation 1 in the EstimateS manual.
    """
    return observed + singles**2/float(doubles*2)

def chao1_bias_corrected(observed, singles, doubles):
    """Calculates bias-corrected chao1 given counts: Eq. 2 in EstimateS manual.

    Formula: chao1 = S_obs + N_1(N_1-1)/(2*(N_2+1)) where N_1 and N_2 are
    count of singletons and doubletons respectively.

    Note: this is the bias-corrected formulat from Chao 1987, Eq. 2 in the
    EstimateS manual.
    """
    return observed + singles*(singles-1) / (2.0*(doubles+1))

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
    r = float(singles)/doubles
    return doubles*(.5*r**2 + r**3 + .24*r**4)

def chao1_var_bias_corrected(singles, doubles):
    """Calculates chao1 variance, bias-corrected.
    
    From EstimateS manual, equation 6.
    """
    s, d = float(singles), float(doubles)
    return s*(s-1)/(2*(d+1)) + (s*(2*s-1)**2)/(4*(d+1)**2) + \
        (s**2 * d * (s-1)**2)/(4*(d+1)**4)

def chao1_var_no_doubletons(singles, chao1):
    """Calculates chao1 variance in absence of doubletons.

    From EstimateS manual, equation 7.

    chao1 is the estimate of the mean of Chao1 from the same dataset.
    """
    s = float(singles)
    return s*(s-1)/2 + s*(2*s-1)**2/4 - s**4/(4*chao1)

def chao1_var_no_singletons(n, observed):
    """Calculates chao1 variance in absence of singletons. n = # individuals.

    From EstimateS manual, equation 8.
    """
    o = float(observed)
    return o*exp(-n/o)*(1-exp(-n/o))
    
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
    #if no diff betweeh chao and observed, CI is just point estimate of observed
    if T == 0:
        return observed, observed
    K = exp(abs(zscore)*sqrt(log(1+(var_chao/T**2))))
    return observed + T/K, observed + T*K

def chao_confidence_no_singletons(n, observed, zscore=1.96):
    """Calculates confidence bounds for chao1/chao2 in absence of singletons.

    Uses Eq. 14 of EstimateS manual.

    n = number of individuals, observed = number of species.
    """
    s = float(observed)
    P = exp(-n/s)
    return max(s, s/(1-P)-zscore*sqrt((s*P/(1-P)))), \
        s/(1-P) + zscore*sqrt(s*P/(1-P))

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

def chao1_lower(*args,**kwargs):
    """Convenience wrapper for chao1_confidence to fit overall diversity API."""
    return chao1_confidence(*args,**kwargs)[0]

def chao1_upper(*args,**kwargs):
    """Convenience wrapper for chao1_confidence to fit overall diversity API."""
    return chao1_confidence(*args,**kwargs)[1]

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
        return freq_counts[1:rare_threshold+1].sum()
        
    def species_abundant(freq_counts, rare_threshold):
        """freq_counts number of abundant species. Default value of abundant is
        greater than 10 individuals. Based on Chao 2000 in Statistica Sinica 
        pg.229  citing observations by Chao, Ma and Yang in 1993."""
        return freq_counts[rare_threshold+1:].sum()

    def number_rare(freq_counts, gamma=False):
        """Number of individuals in rare species. gamma=True generates the
        n_rare used for the variation coefficient."""
        
        n_rare=0 
        if gamma == True:
            for i, j in enumerate(freq_counts[:rare_threshold+1]):
                n_rare = n_rare + (i*j)*(i-1)
            return n_rare
        
        for i, j in enumerate(freq_counts[:rare_threshold+1]):
            n_rare = n_rare + (i*j)
        return n_rare
    
    # calculations begin
    
    freq_counts = frequency_counter(count)
    
    if freq_counts[1:rare_threshold].sum() == 0:
        return species_abundant(freq_counts, rare_threshold)

    if freq_counts[1] == freq_counts[1:rare_threshold].sum():
        raise ValueError("only rare species are singletons, ACE "+\
            "metric is undefined. EstimateS suggests using bias corrected Chao1")
    
    s_abun = species_abundant(freq_counts, rare_threshold) 

    
    s_rare = species_rare(freq_counts, rare_threshold)

    n_rare = number_rare(freq_counts)

    c_ace = 1 - (freq_counts[1]).sum()/float(n_rare)

    top = s_rare*number_rare(freq_counts, gamma=True)
    bottom = c_ace*n_rare*(n_rare-1.0)
    
    gamma_ace = (top/bottom) - 1.0
    
    if 0 > gamma_ace:
        gamma_ace = 0

    return s_abun + (s_rare/c_ace) + ((freq_counts[1]/c_ace)*gamma_ace)


def diversity(indices,f=chao1,step=1,start=None,verbose=False):
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
    for j in range(start, len(indices)+1, step):
        freqs = counts(indices[i:j], freqs)
        try:
            curr = f(freqs)
        except (ZeroDivisionError, FloatingPointError, e):
            curr = 0
        if verbose:
            print curr
        result.append(curr)
        i=j
    return array(result)

