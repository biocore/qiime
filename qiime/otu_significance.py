#!/usr/bin/env python
# File created on 13 Aug 2013
from __future__ import division

__author__ = "Will Van Treuren, Luke Ursell"
__copyright__ = "Copyright 2013, The QIIME project"
__credits__ = ["Will Van Treuren", "Luke Ursell", "Catherine Lozupone"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Luke Ursell"
__email__ = "lkursell@gmail.com"

from qiime.parse import parse_mapping_file_to_dict
from numpy import (array, argsort, vstack, isnan, inf, nan, apply_along_axis,
                   mean, zeros, isinf, logical_or)

from qiime.stats import (fisher_population_correlation,
                                  pearson, spearman, g_fit, ANOVA_one_way, 
                                  kruskal_wallis, mw_t, mw_boot, t_paired, 
                                  mc_t_two_sample, t_two_sample, fisher, 
                                  kendall, assign_correlation_pval, cscore)

from qiime.util import biom_taxonomy_formatter
from collections import defaultdict
from itertools import izip

"""
Library for group_significance.py.
The code in this library is based around two central frameworks. For the group
significance tests the framework is as follows:
The OTU table is a rowXcolumn (otuXsample) matrix. The mapping file specifies
certain groups of samples based on metadata (eg. samples 1,2,6,9 are obese mice
and samples 3,4 5,12 are lean mice). The code slices the OTU matrix into rows
(individual otus) and groups together columns (samples) which have the same
metadata value based on the passed metadata field. These groupings are then
compared using the specified test.
For the gradient correlation tests the framework is as follows:
Each row of the OTU table is correlated with a numeric value found in the some
metadata category of the mapping file. No grouping within categories occurs.
Some abbreviations that are used in this code are:
pmf - parsed mapping file. Nested dict created by parse_mapping_file_to_dict
 which has as top level keys the sample IDs and then all assocaited metadata
 in a dictionary with the mapping file column headers as keys.
bt - biom table object.
row - row is used in several places (eg. row_generator). The 'row' being
 returned is actually a list of arrays that comes from a single row, and a
 collection of grouped columns from the OTU table (for the group significance
 tests). For the gradient correlations a row is actually a full row of the OTU.
"""

# Pursuant to cogent/qiime coding guidelines, globals are uppercase. These dicts
# map the script interface names to the actual functions running these tests.
CORRELATION_TEST_CHOICES = {'pearson': pearson, 'spearman': spearman,
                            'kendall': kendall, 'cscore': cscore}

GROUP_TEST_CHOICES = {'ANOVA': ANOVA_one_way, 'g_test': g_fit,
                      'kruskal_wallis': kruskal_wallis,
                      'parametric_t_test': t_two_sample,
                      'nonparametric_t_test': mc_t_two_sample,
                      'mann_whitney_u': mw_t,
                      'bootstrap_mann_whitney_u': mw_boot}

TWO_GROUP_TESTS = ['parametric_t_test', 'nonparametric_t_test',
                   'mann_whitney_u', 'bootstrap_mann_whitney_u']

# these are the available correlation pvalue calculation methods. kendall is
# appropriate only for kendall, while the other methods are appropriate for
# any metric.
CORRELATION_PVALUE_CHOICES = [
    'parametric_t_distribution', 'fisher_z_transform',
    'bootstrapped', 'kendall']

# Functions for group significance testing


def get_sample_cats(pmf, category):
    """Create {SampleID:category_value} for samples in parsed mf dict.

    Inputs:
     pmf - parsed mapping file. Described at top of library.
     category - string, key in the pmf.
    """
    # ignore samples where the value in the mapping file is empty
    return {k: pmf[k][category] for k in pmf if pmf[k][category] != ""}


def get_cat_sample_groups(sam_cats):
    """Create {category_value:[samples_with_that_value]} dict.

    Inputs:
     sam_cats - dict, output of get_sample_cats.
    """
    cat_sam_groups = defaultdict(list)
    for k, v in sam_cats.iteritems():
        cat_sam_groups[v].append(k)
    return cat_sam_groups


def get_sample_indices(cat_sam_groups, bt):
    """Create {category_value:index_of_sample_with_that_value} dict.

    Inputs:
     cat_sam_groups - dict, output of get_cat_sample_groups.
     bt - biom table object. Described at top of library.
    """
    cat_sam_indices = defaultdict(list)
    for k, v in cat_sam_groups.iteritems():
        cat_sam_indices[k] = [bt.index(i, axis='sample') for i in v]
    return cat_sam_indices


def group_significance_row_generator(bt, cat_sam_indices):
    """Produce generator that feeds lists of arrays to group significance tests.

    Read library documentation for description of what a 'row' is.
    Inputs:
     bt - biom table object. Described at top of library.
     cat_sam_indices - dict, output of get_sample_indices.
    """
    data = array([i for i in bt.iter_data(axis='observation')])
    indices = cat_sam_indices.values()  # list of lists of column indices
    return izip(*[data.take(i, axis=1) for i in indices])


def run_group_significance_test(data_generator, test, test_choices, reps=1000):
    """Run any of the group significance tests.

    Inputs:
     data_generator - generator object, output of row_generator. The output of
      each iter of the data_generator is a list of arrays which is fed to one
      of the tests.
     test - string, key of group_test_choices. the script interface name for the
      functions.
     test_choices - dictionary, defined as global at top of library.
     reps - int, number of reps or permutations to do for the bootstrapped
      tests.
    Ouputs are lists of test statistics, p values, and means of each group.
    """
    pvals, test_stats, means = [], [], []
    for row in data_generator:
        if test == 'nonparametric_t_test':
            test_stat, _, _, pval = test_choices[test](row[0], row[1],
                                                       permutations=reps)
        elif test == 'bootstrap_mann_whitney_u':
            test_stat, pval = test_choices[test](row[0], row[1], num_reps=reps)
        elif test in ['parametric_t_test', 'mann_whitney_u']:
            test_stat, pval = test_choices[test](row[0], row[1])
        else:
            # ANOVA, kruskal_wallis, g_fit will get caught here
            try:
                test_stat, pval = test_choices[test](row)
            except (ValueError):
                test_stat = nan
                pval = nan
        test_stats.append(test_stat)
        pvals.append(pval)
        means.append([i.mean() for i in row])
    return test_stats, pvals, means


def group_significance_output_formatter(bt, test_stats, pvals, fdr_pvals,
                                        bon_pvals, means, cat_sample_indices, md_key):
    """Format the output for gradient tests so it can be easily written.

    Inputs are lists of test statistics, pvalues, fdr corrected pvalues,
    bonferonni corrected pvalues, group means, and the dict of
    {category:sample_index}, and the key to use to extract metadata from the
    biom table. Output is a list of lines.
    """
    header = ['OTU', 'Test-Statistic', 'P', 'FDR_P', 'Bonferroni_P'] + \
        ['%s_mean' % i for i in cat_sample_indices] + [md_key]
    num_lines = len(pvals)
    lines = ['\t'.join(header)]
    for i in range(num_lines):
        tmp = [bt.ids(axis='observation')[i], test_stats[i], pvals[i], fdr_pvals[i],
               bon_pvals[i]] + means[i]
        lines.append('\t'.join(map(str, tmp)))
    # attempt to add metadata
    nls = _add_metadata(bt, md_key, lines)
    return nls


# Functions for gradient correlation testing

def grouped_correlation_row_generator(bt, pmf, category, gc_to_samples):
    """Create generator for grouped correlation tests.

    Inputs:
     bt - biom table object.
     pmf - parsed mapping file (nested dict).
     category - string, header in pmf that will be used to extract continuous
      metadata for gradient correlations.
     gc_to_samples - dict, output of get_cat_sample_groups. Keys are grouping
      categories for a sample group, values are those samples.
    Output:
     tuple of category values (in order of computation), the metadata values,
      and the otu values.
    """
    data = array([i for i in bt.iter_data(axis='observation')])
    category_values = gc_to_samples.keys()
    samples = gc_to_samples.values()
    sample_inds = [[bt.index(i, axis='sample') for i in group] for group in samples]
    try:
        md_vals = []
        for grp in samples:
            md_vals.append(array([pmf[s][category] for s in grp], dtype=float))
    except ValueError:
        raise ValueError("Couldn't convert sample metadata to float.")
    otu_vals = [data.take(inds, 1) for inds in sample_inds]
    return category_values, md_vals, otu_vals


def run_grouped_correlation(md_vals, otu_arrays, test, test_choices,
                            pval_assignment_method, permutations=None):
    """Run grouped correlation test

    This function runs the grouped correlation test. Briefly, it ingests the
    metadata values, the arrays of otu values that are to be correlated with
    them, and the test and pvalue assignment method to use. It calculates the
    individual correlation coefficients for each group (specified implicitly
    by the grouping and ordering of md_vals and otu_arrays) and then it combines
    the corrcoeffs and the pvalues with methods by Fisher.
    Inputs:
     md_vals - list of 1d arrays, continuous metadata to be correlated.
     otu_arrays - list of 1d, otu abundances to be correlated.
     test - str, one of CORRELATION_TEST_CHOICES keys.
     test_choices - dict, CORRELATION_TEST_CHOICES.
     pval_assignment_method - str, one of CORRELATION_PVALUE_CHOICES.
     permutations - int or None, number of permutations to use for bootstrapped
      methods.
    """
    test_fn = test_choices[test]
    sample_sizes = map(len, md_vals)

    def _rho(otu_vals, md_vals):
        return test_fn(otu_vals, md_vals)
    # find the correlations. rhos is list of 1D arrays.
    rhos = []
    for i in range(len(md_vals)):
        rhos.append(apply_along_axis(_rho, 1, otu_arrays[i], md_vals[i]))
    pvals = []
    for i, group_rhos in enumerate(rhos):
        pvals_i = zeros(len(group_rhos))
        for j, rho in enumerate(group_rhos):
            pvals_i[j] = assign_correlation_pval(rho, sample_sizes[i],
                                                 pval_assignment_method, permutations, test_fn, otu_arrays[
                                                     i][j],
                                                 md_vals[i])
        pvals.append(array(pvals_i))
    # calculate combined stats
    fisher_pvals = apply_along_axis(fisher, 0, array(pvals))
    fisher_rho_and_h = apply_along_axis(fisher_population_correlation, 0,
                                        array(rhos), sample_sizes)
    return (
        (rhos, pvals, fisher_pvals, fisher_rho_and_h[0], fisher_rho_and_h[1])
    )


def grouped_correlation_formatter(bt, rhos, pvals, f_rhos, f_pvals, f_hs,
                                  grouping_category, category_values, md_key):
    """Format output from longitudinal tests to be written.

    Inputs are biom table, list of test statistics.
    """
    header = ['OTU'] + \
        ['Rho_%s:%s' % (grouping_category, c) for c in category_values] + \
        ['Pval_%s:%s' % (grouping_category, c) for c in category_values] + \
        ['Fisher population correlation', 'Fisher combined p',
         'Homogeneity pval', md_key]
    num_lines = len(f_rhos)
    lines = ['\t'.join(header)]
    for i in range(num_lines):
        tmp_rhos = [x[i] for x in rhos]
        tmp_pvals = [x[i] for x in pvals]
        tmp = [bt.ids(axis='observation')[i]] + tmp_rhos + tmp_pvals + \
            [f_rhos[i]] + [f_pvals[i]] + [f_hs[i]]
        lines.append('\t'.join(map(str, tmp)))
    nls = _add_metadata(bt, md_key, lines)
    return nls


def correlation_row_generator(bt, pmf, category):
    """Produce a generator that feeds lists of arrays to any gradient test.

    In this function, a row is a full row of the OTU table, a single 1D array.
    Inputs:
     bt - biom table object. Described at top of library.
     cat_sam_indices - dict, output of get_sample_indices.
     category - str, category to pull continuous sample metadata from.
    """
    data = array([i for i in bt.iter_data(axis='observation')])
    # ensure that the order of the category vector sample values is the same
    # as the order of the samples in data. otherwise will have hard to
    # diagnose correspondence issues
    try:
        cat_vect = array([pmf[s][category] for s in bt.ids()], dtype=float)
        return ((row, cat_vect) for row in data)
    except ValueError:
        raise ValueError("Mapping file category contained data that couldn't " +
                         "be converted to float. Can't continue.")


def run_correlation_test(data_generator, test, test_choices,
                         pval_assignment_method, permutations=None):
    """Run correlation tests.

    Inputs:
     data_generator - generator from correlation_row_generator, basically a list
      of tuples where each tuple contains two arrays.
     test - str, one of CORRELATION_TEST_CHOICES keys.
     test_choices - dict, CORRELATION_TEST_CHOICES.
     pval_assignment_method - str, one of CORRELATION_PVALUE_CHOICES.
     permutations - int or None, number of permutations to use for bootstrapped
      methods.
    """
    corr_coefs, pvals = [], []
    test_fn = test_choices[test]
    for otu_vals, md_vals in data_generator:
        r = test_fn(otu_vals, md_vals)
        if pval_assignment_method == 'bootstrapped':
            pval = assign_correlation_pval(r, len(otu_vals),
                                           pval_assignment_method, permutations, test_fn, otu_vals,
                                           md_vals)
        else:
            pval = assign_correlation_pval(r, len(otu_vals),
                                           pval_assignment_method)
        corr_coefs.append(r)
        pvals.append(pval)
    return corr_coefs, pvals


def correlate_output_formatter(bt, test_stats, pvals, fdr_pvals, bon_pvals,
                               md_key):
    '''Produce lines for a tab delimited text file for correlations.py.

    Paramaters
    ----------
    bt : biom table object
    test_stats : array-like
        Floats representing correlation coefficients or paired t test
        statistics.
    pvals : array-like
        Floats representing pvalues for given correlation coefficients.
    fdr_pvals : array-like
        Floats representing FDR corrected pvals.
    bon_pvals : array-like
        Floats representing Bonferroni corrected pvals.
    md_key : str or None
        Key for extracting feature metadata from biom table.

    Returns
    -------
    list of strs
    '''
    header = ['Feature ID', 'Test stat.', 'pval', 'pval_fdr', 'pval_bon',
              md_key]
    num_lines = len(test_stats)
    lines = ['\t'.join(header)]
    for i in range(num_lines):
        tmp = [bt.ids(axis='observation')[i], test_stats[i], pvals[i],
               fdr_pvals[i], bon_pvals[i]]
        lines.append('\t'.join(map(str, tmp)))
    nls = _add_metadata(bt, md_key, lines)
    return nls

def run_paired_t(bt, s1, s2):
    '''Perform a paired t test between samples.

    Parameters
    ----------
    bt : biom table object
        Table containing data for samples in s1 and s2.
    s1 : list
        List of sample ids found in bt (strs).
    s2 : list
        List of sample ids found in bt (strs).

    Returns
    -------
    test_stats : list
        Floats representing the test statistic for each paired comparison. 
    pvals : list
        Floats representing the p-values of the reported test statistics.
    '''
    test_stats = []
    pvals = []
    s1_indices = [bt.index(i, axis='sample') for i in s1]
    s2_indices = [bt.index(i, axis='sample') for i in s2]

    for data in bt.iter_data(axis='observation'):
        test_stat, pval = t_paired(data.take(s1_indices), data.take(s2_indices))
        test_stats.append(test_stat)
        pvals.append(pval)
    return test_stats, pvals

def sort_by_pval(lines, ind):
    """Sort lines with pvals in descending order.

    ind is the index of each line, split on \t, that is to be used for sorting.
    """
    def _nan_safe_sort(line):
        """Sort lines based on pvals and force nans to have inf pval."""
        val = float(line.split('\t')[ind])
        if not isnan(val):
            return val
        else:
            return inf
    return [lines[0]] + sorted(lines[1:], key=_nan_safe_sort)


def _add_metadata(bt, md_key, lines):
    """Add metadata to formatted correlation output lines."""
    taxonomy_md = biom_taxonomy_formatter(bt, md_key)
    if taxonomy_md is not None:
        for i in range(len(lines) - 1):  # one more line than OTU
            # skip header line in lines
            lines[i + 1] = lines[i + 1] + '\t' + taxonomy_md[i]
        return lines
    else:  # remove md_header from the first line
        nls = ['\t'.join(lines[0].split('\t')[:-1])] + lines[1:]
        return nls

def is_computable_float(v):
    '''Return float if v can be converted to float excluding nan and inf.

    Parameters
    ----------
    v : variable
        Value to be converted to float if possible.

    Returns
    -------
    If v can be converted to a float that is not nan or inf, return v.
    Otherwise return False.
    '''
    tmp = float(v)
    if not logical_or(isnan(tmp), isinf(tmp)):
        return tmp 
    else:
        return False
