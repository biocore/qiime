#!/usr/bin/env python

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project" 
__credits__ = ["Jens Reeder", "Rob Knight"]
#remember to add yourself if you make changes
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"
__status__ = "Development"

"""
Contains methods to compute the conditional uncovered probability
from Lladser, Gouet, and Reeder, 
"Extrapolation of Urn Models via Poissonization: Accurate
Measurements of the Microbial Unknown" PLoS 2011.

A lot of this might migrate into cogent at some point, 
e.g. robbins is already in there.

"""
from sys import exit, stderr
from numpy import array, sqrt
from numpy.random import gamma, random, shuffle
from math import factorial
from operator import mul
from collections import defaultdict

from qiime.alpha_diversity import AlphaDiversityCalc, AlphaDiversityCalcs

import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')

from cogent.maths.stats.alpha_diversity import expand_counts, singles, doubles
import cogent.maths.stats.alpha_diversity as alph 

def lladser_point_estimates(sample, r=10):
    """Series of point estimates of the conditional uncovered probability

    sample: series of random observations  
    r: Number of new colors that are required for the next prediction

    This is the point estimator described in Theorem 2 (i):

    Returns: Each new color yields 3-tuple:
         - point estimate
         - position in sample of prediction
         - random variable from poisson process (mostly to make testing easier)
    """

    if(r<=2):
        raise ValueError, "r must be >=3"
    for count, seen, cost, i in get_interval_for_r_new_species(sample, r):
        t = gamma(count, 1)
        point_est = (r-1) / t
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
        curr = seq[i]   #note: first iteration is after looking at first char
        # bail out if there's nothing new
        if curr in seen:
            continue
        else:
            seen.add(curr)
      
        #otherwise, need to see distance to get k colors
        unseen = 0
        j = i + 1
        while unseen < r and j < seq_len:
            if seq[j] not in seen:
                unseen += 1
            j += 1 #note: increments after termination condition

        count = j-i-1 # the interval to see r new colors
        cost = j # the position in seq after seeing r new ones 
        if (not count) or (unseen < r): #bail out if not enough unseen
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
    n  = counts.sum()
    z  = 1.959963985
    W  = (n1*(n-n1) + 2*n*n2)/(n**3)

    return  n1/n + z*sqrt(W), n1/n - z*sqrt(W) 

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
    #precompute since it is constant for chosen m
    fact_m_1 = factorial(m-1)
    for n in range(1,len(sample)):
        sum=0
        counts[sample[n-1]]+=1
        count = array(counts.values())
        filtered_count = array([x for x in count if x<=(m+1)])
        for k in range(1,m+1):
            k_tons = (filtered_count==k).sum()

            if k_tons!=0:
                #computing factorials for large number is inefficient
                #according to manuel (and some tests) the above formula is
                # equal to:
                # increment = binomial_coeff(m-1, k-1)/binomial_coeff(n+m,k)   
                increment = k * fact_m_1/factorial(m-k)  \
                              * (1/reduce(mul, range(n+m-k+1,n+m+1)))
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
    """

    for count,seen,cost,i in get_interval_for_r_new_species(seq, r):
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

    Returns: A confidence interval that contains the true conditional
             uncovered probability with a probability of 100% * alpha.
    """

    if ci_type=='U':
        upper_bound = upper_confidence_bound(r, alpha) / t
        return (0,upper_bound)
    elif ci_type=='L':
        lower_bound = lower_confidence_bound(r, alpha) / t
        return (lower_bound,1)

    bound_params = ul_confidence_bounds(f, r, alpha)
    if ci_type=='ULCL':
        bound_param = bound_params[0]
    elif ci_type=='ULCU':
        bound_param = bound_params[1]
    else:
        raise ValueError, "Unknown ci_type: %s" % (ci_type)

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

    Returns: Constant c such that the confidence interval is [c/T_r, 1]
    """
   
    alpha = round(alpha,2)
    if not (alpha == 0.95):
        raise ValueError, "alpha must be 0.95"
    data = { 1:  0.051293294,
             2:  0.355361510,
             3:  0.817691447,
             4:  1.366318397,
             5:  1.970149568,
             6:  2.613014744,
             7:  3.285315692,
             8:  3.980822786,
             9:  4.695227540,
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
        raise ValueError, "r must be between 1,..,25 or 50"

def upper_confidence_bound(r, alpha):
    """Get constant for confidence interval with lower bound fixed to 0
    
    r: number of new colors
    alpha: Confidence interval (for 95% conf use 0.95)

    Compute constant b according to Theorem 2 iii with a=0
    aka c_0 from Table 3

    Returns: Constant c such that the confidence interval is [0,c/T_r]
    """
   
    alpha = round(alpha,2)
    if not (alpha == 0.95):
        raise ValueError, "alpha must be 0.95"

    data = { 1:  2.995732274,
             2:  4.743864518,
             3:  6.295793622,
             4:  7.753656528,
             5:  9.153519027,
             6:  10.51303491,
             7:  11.84239565,
             8:  13.14811380,
             9:  14.43464972,
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
        raise ValueError, "r must be between 1,..,25 or 50"

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
    
    alpha = round(alpha,2)
    # Hack in some special values we used for the paper.
    # Since Manuel needs to compute those semi-automatically
    # using Maple, we pre-calculate only a few common ones

    if f==2 and r==50 and alpha==0.95:
        return (31.13026306, 38.94718565)
    elif f==2 and r==33 and alpha==0.95:
        return (22.3203508, 23.4487304)
    elif f==1.5 and r==100 and alpha==0.95:
        return (79.0424349, 83.22790086)
    elif f==1.5 and r==94 and alpha==0.95:
        return (75.9077267, 76.5492088)
    elif f==2.5 and r==19 and alpha==0.95:
        return (11.26109001, 11.96814857)

    #In the next block for each f, we report the smallest possible value of r
    #from table 4 in the paper
    if f==80 and r==2 and alpha==0.95:
        return(0.0598276655, 0.355361510)
    elif f==48 and r==2 and alpha==0.95:
        return(0.1013728884, 0.355358676)
    elif f==40 and r==2 and alpha==0.95:
        return(0.1231379857, 0.355320458)
    elif f==24 and r== 2 and alpha==0.95:
        return( 0.226833483, 0.346045204)
    elif f==20 and r==3 and alpha==0.95:
        return(0.320984257, 0.817610455)
    elif f==12 and r==3 and alpha==0.95:
        return(0.590243030, 0.787721610)
    elif f==10 and r==4 and alpha==0.95:
        return(0.806026244, 1.360288674)
    elif f==6 and r==6 and alpha==0.95:
        return(1.8207383, 2.58658608)
    elif f==5 and r==7 and alpha==0.95:
        return(2.48303930, 3.22806682)
    elif f==3 and r==14 and alpha==0.95:
        return(7.17185045, 8.27008349)
    elif f==2.5 and r==19 and alpha==0.95:
        return(11.26109001, 11.96814857)
    elif f==1.5 and r==94 and alpha==0.95:
        return(75.9077267, 76.5492088)
    elif f==1.25 and r==309 and alpha==0.95:
        return(275.661191, 275.949782)

    # all others combination are only computed for f=10
    # and alpha = 0.90, 0.95 and 0.99
    elif f==10 and r<=50:
        try:
            (a,b) = cbs[round(alpha,2)][r]
        except KeyError:
            raise ValueError, "No constants are precomputed for the combination of " \
                + "f:%f, r:%d, and alpha:%.2f"% (f,r,alpha)
        if (a!=None and b!=None):
           return(a,b) 
    
    raise ValueError, "No constants are precomputed for the combination of " \
         + "f:%f, r:%d, and alpha:%.2f"% (f,r,alpha)

###### 
# Below are the values used for Theorem 3 iii
# Values hand computed by Manuel Lladser using Maple

cb_90 = [
    (None, None), # 0, makes indexing easier 
    (None, None), # no feasible solution
    (None, None), # no feasible solution
    (.5635941995, 1.095834700),
    (.6764656264, 1.744588615),
    (.8018565594, 2.432587343),
    (.9282215025, 3.151897973),
    (1.053433716, 3.894766804),
    (1.177158858, 4.656118177),
    (1.299491033, 5.432468058),
    (1.420604842, 6.221304605), # 10
    (1.540665805, 7.020746595),
    (1.659812701, 7.829342026),
    (1.778158703, 8.645942495),
    (1.895796167, 9.469621185),
    (2.012801198, 10.29961731),
    (2.129237257, 11.13529724),
    (2.245157877, 11.97612664),
    (2.360608695, 12.82164994),
    (2.475628991, 13.67147502),
    (2.590252861, 14.52526147), # 20
    (2.704510123, 15.38271151),
    (2.818427036, 16.24356290),
    (2.932026869, 17.10758326),
    (3.045330351, 17.97456551),
    (3.158356050, 18.84432420),
    (None, None), # not computed
    (None, None),
    (None, None),
    (None, None),
    (3.719850286, 23.22944415), #30
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (4.828910181, 32.13892224), #40
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (5.924900191, 41.17906791) #50
]

cb_95 = [
    (None, None), # 0
    (None, None), 
    (None, None),
    (None, None),
    (.8060262438, 1.360288674), # 4
    (.9240311584, 1.969902537),
    (1.053998892, 2.613007253),
    (1.185086998, 3.285315518),
    (1.315076337, 3.980822783),
    (4.695227540, 4.695227541),
    (1.570546801, 5.425405698), # 10
    (1.696229569, 6.169007289),
    (1.820753729, 6.924212513),
    (1.944257622, 7.689578291),
    (2.066857113, 8.463937522),
    (2.188648652, 9.246330491),
    (2.309712994, 10.03595673),
    (2.430118373, 10.83214036),
    (2.549923010, 11.63430451),
    (2.669177032, 12.44195219),
    (2.787923964, 13.25465160), # 20
    (2.906201884, 14.07202475),
    (3.024044329, 14.89373854),
    (3.141481020, 15.71949763),
    (3.258538445, 16.54903871),
    (3.375240327, 17.38212584),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (3.954097220, 21.59397923), #30
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (5.093973695, 30.19573919), #40
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (6.217105673, 38.96473258) #50
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
    (1.878809285, 4.130199076), #10
    (2.014632329, 4.771246173),
    (2.149044735, 5.428180734),
    (2.282101533, 6.099073460),
    (2.413917374, 6.782354878),
    (2.544610844, 7.476728267),
    (2.674289153, 8.181107778),
    (2.803045614, 8.894573463),
    (2.930960779, 9.616337916),
    (3.058104355, 10.34572103),
    (3.184536992, 11.08213063), #20
    (3.310311816, 11.82504734),
    (3.435475649, 12.57401269),
    (3.560070013, 13.32861956),
    (3.684131925, 14.08850448),
    (3.807694563, 14.85334135),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (4.41897094, 18.7424258), #30
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (5.61643962, 26.7700386),#  40
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (None, None),
    (6.79033616, 35.0324474) #50
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
        pe = ('NaN','NaN')
    return pe

def robbins(counts, **argv):
    return(alph.robbins(counts))

esty_ci.return_names = ('esty_lower_bound', 'esty_upper_bound')
lladser_ci.return_names = ('lladser_lower_bound', 'lladser_upper_bound')
cup_metrics = [lladser_pe,
               lladser_ci,
               #starr, not yet needs tests
               esty_ci,
               robbins
               ]

def get_cup_metric(name):
    """Gets metric by name from list in this module
    """
    for metric in cup_metrics:
        if metric.__name__.lower() == name.lower():
            return metric    
    raise AttributeError

def list_known_metrics():
    """Show the names of available metrics."""
    return [ metric.__name__ for metric in cup_metrics ]

def single_file_cup(otu_filepath, metrics, outfilepath, r, alpha, f, ci_type):
    """Compute variations of the conditional uncovered probability.

    otufilepath: path to otu_table file
    metrics: comma separated list of required metrics
    outfilepath: path to output file

    r: Number of new colors that are required for the next prediction
    alpha: desired confidence level
    f: ratio between upper and lower bound
    ci_type: type of confidence interval. One of:
             ULCL: upper and lower bounds with conservative lower bound
             ULCU: upper and lower woth conservative upper bound
             U: Upper bound only, lower bound fixed to 0
             L: Lower bound only, upper bound fixed to 1

    The opposite of uncovered probability is sometimes called coverage.
    """
    metrics_list = metrics.split(',')
    calcs = []

    params = {'r': r,
              'alpha': alpha,
              'f':f,
              'ci_type':ci_type}
                  
    for metric in metrics_list:
        try:
            metric_f = get_cup_metric(metric)
        except AttributeError:
            stderr.write(
                "Could not find metric %s.\n Known metrics are: %s\n" \
                    % (metric, ', '.join(list_known_metrics())))
            exit(1)
            
        c = AlphaDiversityCalc(metric_f, params=params)
        calcs.append(c)
    
    all_calcs = AlphaDiversityCalcs(calcs)

    try:
        result = all_calcs(data_path=otu_filepath,
            result_path=outfilepath, log_path=None)
        if result:  #can send to stdout instead of file
            print all_calcs.formatResult(result)
    except IOError, e:
        stderr.write("Failed because of missing files.\n")
        stderr.write(str(e)+'\n')
        exit(1)
