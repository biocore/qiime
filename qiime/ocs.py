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
        [pmf.pop(i) for i in pmf if i not in shared_samples]
        # remove samples in the biom table that were not in the mapping file
        def _f(sv, sid, smd):
            if sid in shared_samples:
                return True
            else:
                return False
        bt = bt.filterSamples(_f)
    return pmf, bt

def get_sample_cats(pmf, category):
    """Create {SampleID:category_value} for samples in parsed mf dict."""
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




