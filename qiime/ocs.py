#!/usr/bin/env python
# File created on 13 Aug 2013
from __future__ import division

__author__ = "Luke Ursell"
__copyright__ = "Copyright 2013, The QIIME project"
__credits__ = ["Will Van Treuren, Luke Ursell"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Luke Ursell"
__email__ = "lkursell@gmail.com"
__status__ = "Development"

from biom.parse import parse_biom_table
from qiime.parse import parse_mapping_file_to_dict
from numpy import array, argsort
from cogent.maths.stats.util import Numbers

"""
Library for otu_category_significance.
"""

two_group_tests = ['parametric_t_test', 'nonparametric_t_test', 
    'mann_whitney_u', 'bootstrap_mann_whitney_u']

def sync_biom_and_mf(pmf, bt):
    """Reduce mapping file dict and biom table to shared samples."""
    mf_samples = set(pmf.keys())
    bt_samples = set(bt.SampleIds)
    if mf_samples == bt_samples:
        # agreement, can continue without fear of breaking code
        pass
    else: 
        shared_samples = mf_samples.intersection(bt_samples)
        # check that we shared something
        assert len(shared_samples)!=0, \
            "sync_biom_and_mf: No shared samples, no point in continuing."
        # tell the user which samples we are excluding
        print "The following samples were not shared, and will not be "+\
            "considered in the analysis:\n" + \
            ', '.join(shared_samples-(mf_samples.union(bt_samples)))
        # remove samples that were in the mapping file but not biom file
        npmf = {k:v for k,v in pmf.items() if k in shared_samples}
        # remove samples in the biom table that were not in the mapping file
        def _f(sv, sid, smd):
            if sid in shared_samples:
                return True
            else:
                return False
        nbt = bt.filterSamples(_f)
    return npmf, nbt

def get_sample_cats(pmf, category):
    """Create {SampleID:category_value} for samples in parsed mf dict."""
    # ignore samples where the value in the mapping file is empty
    return {k:pmf[k][category] for k in pmf.keys() if pmf[k][category] != ""}

def get_cat_sample_groups(sam_cats):
    """Create {category_value:[samples_with_that_value} dict."""
    cat_sam_groups = {group:[] for group in set(sam_cats.values())}
    [cat_sam_groups[v].append(k) for k,v in sam_cats.items()]
    return cat_sam_groups

def get_sample_indices(cat_sam_groups, bt):
    """Create {category_value:index_of_sample_with_that_value} dict."""
    return {k:[bt.SampleIds.index(i) for i in v] for k,v in cat_sam_groups.items()}

def row_generator(bt, cat_sam_indices):
    """Produce a generator that feeds lists of arrays to any test."""
    data = array([bt.observationData(i) for i in bt.ObservationIds])
    return ([row[cat_sam_indices[k]] for k in cat_sam_indices] for row in data)

def run_ocs_test(data_generator, test, test_choices, *args):
    """Run any of the implemented tests."""
    pvals, test_stats, means = [], [], []
    # test choices defined in the ocs.py script
    for row in data_generator:
        if test == 'nonparametric_t_test':
            test_stat, _, _, pval = test_choices[test](row[0], row[1], *args)
        elif test == 'bootstrap_mann_whitney_u':
            test_stat, pval = test_choices[test](row[0], row[1], *args)
        elif test in ['parametric_t_test', 'mann_whitney_u']:
            test_stat, pval = test_choices[test](row[0], row[1])
        else:
            test_stat, pval = test_choices[test](row)
        test_stats.append(test_stat)
        pvals.append(pval)
        means.append([i.mean() for i in row])
    return test_stats, pvals, means

def fdr_correction(probs):
    """corrects a list of probs using the false discovery rate method

    ranks the p-values from low to high. multiplies each p-value by the #
    of comparison divided by the rank.
    """
    corrected_probs = [None] * len(probs)
    for rank, index in enumerate(argsort(probs)):
        correction = len(probs) / float(rank + 1)
        if probs[index]:
            fdr_p = probs[index] * correction
        else:
            fdr_p = 'NA'
        corrected_probs[index] = fdr_p
    return corrected_probs

def bonferroni_correction(probs):
    """Make Bonferroni correction to probs."""
    return array(probs)*len(probs)

def output_formatter(bt, test_stats, pvals, fdr_pvals, bon_pvals, means, 
    cat_sample_indices):
    """Format the output for all tests so it can be easily written."""
    header = ['OTU', 'Test-Statistic', 'P', 'FDR_P', 'Bonferroni_P']
    header += ['%s_mean' % i for i in cat_sample_indices.keys()]
    # find out if bt came with taxonomy. this could be improved
    if bt.ObservationMetadata is None:
        include_taxonomy = False
    else:
        include_taxonomy = True
        header += ['Taxonomy']
    num_lines = len(pvals)
    lines = ['\t'.join(header)]
    for i in range(num_lines):
        tmp = [bt.ObservationIds[i], test_stats[i], pvals[i], fdr_pvals[i], 
            bon_pvals[i]] + means[i] 
        if include_taxonomy:
            tmp += ';'.join(bt.ObservationMetadata[i].values())
        lines.append('\t'.join(map(str, tmp)))
    return lines

def sort_by_pval(lines):
    """Sort lines with pvals in descending order."""
    # output_formatter will always put pvals in index 2
    return [lines[0]] + sorted(lines[1:], key=lambda x: float(x.split('\t')[2]))

##########
##########
# Functions for gradient correlation testing
##########
##########

def correlation_row_generator(bt, pmf, category):
    """Produce a generator which will feed correlation tests rows."""
    data = array([bt.observationData(i) for i in bt.ObservationIds])
    # ensure that the order of the category vector is the same as the order of
    # the samples, otherwise will have hard to diagnose correspondence issues 
    category_vector = \
        array([pmf[s][category] for s in bt.SampleIds]).astype(float)
    return ((row,category_vector) for row in data)

def parametric_correlation_significance(test_stat, n):
    """Calculate the significance of a Pearson or Spearman r, rho.

    Notes: ignoring the possibility of using something other than a two tailed
    test since that makes little sense in the context of +- correlation values.
    Our alternate hypothesis is that uncorrelated bivariate data could generate
    a r or rho value as extreme or more extreme, thus we have two tails.

    """
    df = n-2 #degrees of freedom
    if n<3: #need at least 3 samples for students t parametric p value calc
        p_pval = 1.0 #p_pval = parametric p value
    else: 
        try:
            t_stat = test_stat/sqrt((1.-test_stat**2)/float(df))
            p_pval = t_prob(t_stat, df) #we force it to be two tailed
        except (ZeroDivisionError, FloatingPointError):
            # something unpleasant happened, most likely r or rho where +- 1 
            # which means the parametric p val should be 1 or 0 or nan
            p_pval = nan
    return pval 

def nonparametric_correlation_significance(test_stat, test, v1, v2,
    permutations=1000, confidence_level=.95):
    """Calculate the significance of a Pearson or Spearman r, rho.

    Notes: ignoring the possibility of using something other than a two tailed
    test since that makes little sense in the context of +- correlation values.
    Our alternate hypothesis is that uncorrelated bivariate data could generate
    a r or rho value as extreme or more extreme, thus we have two tails.

    """
    perm_corr_vals = []
    for i in range(permutations):
        perm_corr_vals.append(test(v1, permutation(v2)))
    # calculate number of bootstrapped statistics which were greater than or 
    # equal to passed test_stat
    return (array(perm_corr_vals) >= test_stat).sum()/float(permutations)

def fisher_confidence_intervals(test_stat, n):
    """Compute the confidence intervals around the test statistic."""
    # compute confidence intervals using fishers z transform
    z_crit = abs(ndtri((1 - confidence_level) / 2.))
    ci_low, ci_high = None, None
    if n > 3:
        try:
            ci_low = tanh(arctanh(test_stat) - (z_crit / sqrt(n - 3)))
            ci_high = tanh(arctanh(test_stat) + (z_crit / sqrt(n - 3)))
        except (ZeroDivisionError, FloatingPointError):
            # r or rho was presumably 1 or -1. Match what R does in this case.
            # feel like nan should be returned here given that we can't make 
            # the calculation
            ci_low, ci_high = test_stat, test_stat
    return ci_low, ci_high

def run_correlation_test(data_generator, test, test_choices):
    """Run correlation tests."""
    corr_coefs, p_pvals, np_pvals, ci_highs, ci_lows = [], [], [], [], []
    for row in data_generator:
        # kendalls tau calculates its own paramteric p value
        if test == 'kendall':
            test_stat, p = test_choices[test](row[0], row[1], return_p=True)
            p_pval = p
        else:
            test_stat = test_choices[test](row[0], row[1])
            p_pval = parametric_correlation_significance(test_stat, len(row[0]))
        np_pval = nonparametric_correlation_significance(test_stat, test, 
            row[0], row[1])
        ci_low, ci_high = fisher_confidence_intervals(test_stat,len(row[0]))
        corr_coefs.append(test_stat)
        p_pvals.append(p)
        np_pvals.append(np_pval)
        ci_lows.append(ci_low)
        ci_highs.append(ci_high)
    return corr_coefs, p_pvals, np_pvals, ci_highs, ci_lows



