#!/usr/bin/env python
# File created on 13 Aug 2013
from __future__ import division

__author__ = "Luke Ursell"
__copyright__ = "Copyright 2013, The QIIME project"
__credits__ = ["Will Van Treuren", "Luke Ursell", "Catherine Lozupone"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Luke Ursell"
__email__ = "lkursell@gmail.com"
__status__ = "Development"

from biom.parse import parse_biom_table
from qiime.parse import parse_mapping_file_to_dict
from numpy import array, argsort, vstack
from cogent.maths.stats.util import Numbers
from qiime.pycogent_backports.test import (parametric_correlation_significance,
    nonparametric_correlation_significance, fisher_confidence_intervals,
    pearson, spearman, kendall_correlation, G_fit, ANOVA_one_way, 
    kruskal_wallis, mw_test, mw_boot, t_paired)
from cogent.maths.stats.test import (t_two_sample, mc_t_two_sample)

"""
Library for test_group_significance.py and test_gradient_correlation.py. The 
code in this library is all based around a central framework. The OTU table is 
a rowXcolumn (otuXsample) matrix. The mapping file specifies certain groupings 
of samples based on metadata (eg. samples 1,2,6,9 are obese mice and samples 3,4
5,12 are lean mice). This code slices the OTU matrix into rows (individual otus)
and groups together columns (samples) which have the same metadata based on the 
passed metadata field. These groupings are then compared using the specified 
test. Some abbreviations that are used in this code are:
pmf - parsed mapping file. Nested dict created by parse_mapping_file_to_dict 
 which has as top level keys the sample IDs and then all assocaited metadata 
 in a dictionary with the mapping file column headers as keys. 
bt - biom table object. Created with parse_biom_table. 
row - row is used in several places (eg. row_generator). The 'row' being 
 returned is actually a list of arrays that comes from a single row, and a 
 collection of grouped columns from the OTU table. 
"""

# Pursuant to cogent/qiime coding guidelines, globals are uppercase. These dicts
# map the script interface names to the actual functions running these tests.
CORRELATION_TEST_CHOICES = {'pearson': pearson, 'spearman': spearman,
    'kendall': kendall_correlation}

GROUP_TEST_CHOICES = {'ANOVA': ANOVA_one_way, 'g_test': G_fit, 
    'kruskal_wallis': kruskal_wallis, 'parametric_t_test': t_two_sample,
    'nonparametric_t_test': mc_t_two_sample, 'mann_whitney_u': mw_test, 
    'bootstrap_mann_whitney_u': mw_boot}

TWO_GROUP_TESTS = ['parametric_t_test', 'nonparametric_t_test', 
    'mann_whitney_u', 'bootstrap_mann_whitney_u']

# Functions for group significance testing

def sync_biom_and_mf(pmf, bt, verbose=True):
    """Reduce mapping file dict and biom table to shared samples.

    Inputs: 
     pmf - parsed mapping file from parse_mapping_file_to_dict (nested dict).
     bt - parse biom table from parse_biom_table (biom table object).
    Outputs are a bt and pmf that contain only shared samples.
    """
    mf_samples = set(pmf.keys())
    bt_samples = set(bt.SampleIds)
    if mf_samples == bt_samples:
        # agreement, can continue without fear of breaking code
        return pmf, bt
    else: 
        shared_samples = mf_samples.intersection(bt_samples)
        # check that we shared something
        assert len(shared_samples)!=0, \
            "sync_biom_and_mf: No shared samples, no point in continuing."
        if verbose:
            print "The following samples were not shared, and will not be "+\
                "considered in the analysis:\n" + \
                ', '.join(mf_samples.union(bt_samples)-shared_samples)
        # remove samples that were in the mapping file but not biom file
        npmf = {k:v for k,v in pmf.items() if k in shared_samples}
        # remove samples in the biom table that were not in the mapping file
        def _f(sv, sid, smd):
            return sid in shared_samples
        nbt = bt.filterSamples(_f)
    return npmf, nbt

def get_sample_cats(pmf, category):
    """Create {SampleID:category_value} for samples in parsed mf dict.

    Inputs:
     pmf - parsed mapping file. Described at top of library.
     category - string, key in the pmf.
    """
    # ignore samples where the value in the mapping file is empty
    return {k:pmf[k][category] for k in pmf.keys() if pmf[k][category] != ""}

def get_cat_sample_groups(sam_cats):
    """Create {category_value:[samples_with_that_value} dict.

    Inputs:
     sam_cats - dict, output of get_sample_cats."""
    cat_sam_groups = {group:[] for group in set(sam_cats.values())}
    [cat_sam_groups[v].append(k) for k,v in sam_cats.items()]
    return cat_sam_groups

def get_sample_indices(cat_sam_groups, bt):
    """Create {category_value:index_of_sample_with_that_value} dict.

    Inputs: 
     cat_sam_groups - dict, output of get_cat_sample_groups.
     bt - biom table object. Described at top of library.
    """
    return {k:[bt.SampleIds.index(i) for i in v] for k,v in cat_sam_groups.items()}

def row_generator(bt, cat_sam_indices):
    """Produce a generator that feeds lists of arrays to any test.

    Read library documentation for description of what a 'row' is. 
    Inputs: 
     bt - biom table object. Described at top of library.
     cat_sam_indices - dict, output of get_sample_indices.
    """
    data = array([bt.observationData(i) for i in bt.ObservationIds])
    return ([row[cat_sam_indices[k]] for k in cat_sam_indices] for row in data)

def run_group_significance_test(data_generator, test, test_choices, *args):
    """Run any of the group significance tests.

    Inputs:
     data_generator - generator object, output of row_generator. The output of 
      each iter of the data_generator is a list of arrays which is fed to one 
      of the tests.
     test - string, key of group_test_choices. the script interface name for the
      functions.
     test_choices - dictionary, defined as global at top of library.
    Ouputs are lists of test statistics, p values, and means of each group.
    """
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
            # ANOVA, kruskal_wallis, G_fit will get caught here
            test_stat, pval = test_choices[test](row)
        test_stats.append(test_stat)
        pvals.append(pval)
        means.append([i.mean() for i in row])
    return test_stats, pvals, means

def fdr_correction(pvals):
    """Adjust pvalues for multiple tests using the false discovery rate method.

    In short: ranks the p-values in ascending order and multiplies each p-value 
    by the number of comparisons divided by the rank of the p-value in the 
    sorted list. Input is list of floats.
    """
    tmp = array(pvals)
    return tmp*tmp.size()/(1.+argsort(tmp).astype(float))

def bonferroni_correction(pvals):
    """Adjust pvalues for multiple tests using the Bonferroni method.

    In short: multiply all pvals by the number of comparisons."""
    return array(pvals)*len(pvals)

def output_formatter(bt, test_stats, pvals, fdr_pvals, bon_pvals, means, 
    cat_sample_indices):
    """Format the output for gradient tests so it can be easily written.

    Inputs are lists of test statistics, pvalues, fdr corrected pvalues, 
    bonferonni corrected pvalues, group means, and the dict of
    {category:sample_index}.
    """
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
            tmp.append(biom_taxonomy_formatter(bt.ObservationMetadata[i]))
        lines.append('\t'.join(map(str, tmp)))
    return lines

def sort_by_pval(lines, ind):
    """Sort lines with pvals in descending order.

    ind is the index of each line, split on \t, that is to be used for sorting.
    """
    return [lines[0]] + \
        sorted(lines[1:], key=lambda x: float(x.split('\t')[ind]))

# Functions for gradient correlation testing

def correlation_row_generator(bt, pmf, category, ref_sample=None):
    """Produce a generator that feeds lists of arrays to any gradient test.

    Read library documentation for description of what a 'row' is. 
    Inputs: 
     bt - biom table object. Described at top of library.
     cat_sam_indices - dict, output of get_sample_indices.
    """
    data = array([bt.observationData(i) for i in bt.ObservationIds])
    if ref_sample is not None:
        # user passed a ref sample to adjust all the other sample OTU values
        # we subtract the reference sample from all data as Cathy did in the 
        # original implementation. 
        data = (data.T - bt.sampleData(ref_sample)).T
    try:
        # ensure that the order of the category vector sample values is the same 
        # as the order of the samples in data. otherwise will have hard to 
        # diagnose correspondence issues 
        category_vector = \
            array([pmf[s][category] for s in bt.SampleIds]).astype(float)
        return ((row,category_vector) for row in data)
    except ValueError:
        raise ValueError("Mapping file category contained data that couldn't "+\
            "be converted to float. Can't continue.")

def run_correlation_test(data_generator, test, test_choices):
    """Run correlation tests."""
    corr_coefs, p_pvals, np_pvals, ci_highs, ci_lows = [], [], [], [], []
    test_fn = test_choices[test]
    for row in data_generator:
        # kendalls tau calculates its own paramteric p value
        if test == 'kendall':
            test_stat, p = test_fn(row[0], row[1], return_p=True)
            p_pval = p
        else: # spearman, pearson executed here
            test_stat = test_fn(row[0], row[1])
            p_pval = parametric_correlation_significance(test_stat, len(row[0]))
        
        np_pval = nonparametric_correlation_significance(test_stat, test_fn, 
            row[0], row[1])
        ci_low, ci_high = fisher_confidence_intervals(test_stat,len(row[0]))
        corr_coefs.append(test_stat)
        p_pvals.append(p_pval)
        np_pvals.append(np_pval)
        ci_lows.append(ci_low)
        ci_highs.append(ci_high)
    return corr_coefs, p_pvals, np_pvals, ci_highs, ci_lows

def correlation_output_formatter(bt, corr_coefs, p_pvals, p_pvals_fdr, 
    p_vals_bon, np_pvals, np_pvals_fdr, np_pvals_bon, ci_highs, 
    ci_lows):
    """Format the output of the correlations for easy writing."""
    header = ['OTU', 'Correlation_Coef', 'parametric_P', 'parametric_P_FDR',
        'parametric_P_Bon', 'nonparametric_P', 'nonparametric_P_FDR',
        'nonparametric_P_Bon', 'confidence_low', 'confidence_high']
    # find out if bt came with taxonomy. this could be improved
    if bt.ObservationMetadata is None:
        include_taxonomy = False
    else:
        include_taxonomy = True
        header += ['Taxonomy']
    num_lines = len(corr_coefs)
    lines = ['\t'.join(header)]
    for i in range(num_lines):
        tmp = [bt.ObservationIds[i], corr_coefs[i], p_pvals[i], p_pvals_fdr[i], 
            p_vals_bon[i], np_pvals[i], np_pvals_fdr[i], np_pvals_bon[i], ci_highs[i], 
            ci_lows[i]]
        if include_taxonomy:
            tmp.append(biom_taxonomy_formatter(bt.ObservationMetadata[i]))
        lines.append('\t'.join(map(str, tmp)))
    return lines

def paired_t_generator(bt, s_before, s_after):
    """Produce a generator to run paired t tests on each OTU."""
    b_data = vstack([bt.sampleData(i) for i in s_before]).T
    a_data = vstack([bt.sampleData(i) for i in s_after]).T
    return ((b_data[i], a_data[i]) for i in range(len(bt.ObservationIds)))

def run_paired_t(data_generator):
    """Run paired t test on data."""
    test_stats, pvals = [], []
    for b_data, a_data in data_generator:
        test_stat, pval = t_paired(b_data, a_data)
        test_stats.append(test_stat)
        pvals.append(pval)
    return test_stats, pvals

def paired_t_output_formatter(bt, test_stats, pvals, fdr_pvals, bon_pvals):
    """Format the output for all tests so it can be easily written."""
    header = ['OTU', 'Test-Statistic', 'P', 'FDR_P', 'Bonferroni_P']
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
            bon_pvals[i]]
        if include_taxonomy:
            tmp.append(biom_taxonomy_formatter(bt.ObservationMetadata[i]))
        lines.append('\t'.join(map(str, tmp)))
    return lines

def biom_taxonomy_formatter(data):
    """Figure out what type of metadata the biom table has, create string."""
    try:
        keys = data.keys()
        if len(keys) > 1:
            raise ValueError("1 < metadata keys. Can't decide which to use.")
        else: 
            md_data = data[keys[0]]
        if type(md_data) == dict:
            return ''.join(['%s_%s' % (k,v) for k,v in md_data.items()])
        elif type(md_data) == list:
            return ';'.join(md_data)
        elif type(md_data) == str:
            return md_data
    except AttributeError:
        raise ValueError('metadata not formatted in a dictionary.')


