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
from numpy import array
from qiime.pycogent_backports.test import ANOVA_one_way, correlation_test
from cogent.maths.stats.util import Numbers

from collections import defaultdict

"""To Do:

"""

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
    return {k:[bt.SampleIds.index(i) for i in v] for k,v in cat_sam_groups.items}

def row_generator(bt, cat_sam_indices):
    """Produce a generator that feeds lists of arrays to any test."""
    data = array([bt.observationData(i) for i in bt.ObservationIds])
    return ([row[cat_sam_indices[k]] for k in cat_sam_indices] for row in data)



###################
# output formatters 



def G_fit_formatter(bt, gs, ps, means, cat_sample_indices):
    """Create list of lines for writing G_fit results."""
    # find out if bt came with taxonomy. this could be improved
    if bt.ObservationMetadata is None:
        header =['OTU', 'G-Stat', 'Pvalue', 'FDR-Pvalue', 'Bonferroni-Pvalue']+\
            ['%i_mean' % i for i in cat_sample_indices.keys()
        include_taxonomy = False
    else: 
        header = header + ['Taxonomy']
        include_taxonomy = True
    # avoid zip; creates a new object, wastes memory and compute. all items of
    # equal length, index is sufficient.
    num_lines = len(gs)
    lines = ['\t'.join(header)]
    for i in range(num_lines):
        nl = '\t'.join(map(str,[bt.ObservationIds[i], gs[i], ps[i]]+means[i]))
        if include_taxonomy == True:
            nl+='\t%s' % ';'.join(bt.ObservationMetadata[i])
        lines.append(nl)
     return lines

def ANOVA_formatter(bt, ps, means, cat_sample_indices):
    """Create list of lines for writing simple ANOVA results."""
    # find out if bt came with taxonomy. this could be improved
    if bt.ObservationMetadata is None:
        header =['OTU', 'Pvalue', 'FDR-Pvalue', 'Bonferroni-Pvalue']+\
            ['%i_mean' % i for i in cat_sample_indices.keys()
        include_taxonomy = False
    else: 
        header = header + ['Taxonomy']
        include_taxonomy = True
    # avoid zip; creates a new object, wastes memory and compute. all items of
    # equal length, index is sufficient.
    num_lines = len(gs)
    lines = ['\t'.join(header)]
    for i in range(num_lines):
        nl = '\t'.join(map(str,[bt.ObservationIds[i], ps[i]]+means[i]))
        if include_taxonomy == True:
            nl+='\t%s' % ';'.join(bt.ObservationMetadata[i])
        lines.append(nl)
     return lines



# def get_group_means(data):
#     """Return mean of each group in data. Format is output of row_generator."""
#     return [i.mean() for i in data]

# def get_group_names(cat_sample_indices):
#     """Return names of groups in order of generation by row_generator."""
#     return cat_sample_indices.keys()

# # order of means is cat_sample_indices.keys()
# def add_means(bt.ObservationIds, taxonomy_key):

#     if include_means:

#     if include_taxonomy:

# G_fit_header = ['OTU', 'G-Stat', 'Pvalue', 'FDR-Pvalue', 'Bonferroni-Pvalue'] +\
#     cat_sample_indices.keys()

# ANOVA_header = ['OTU', 'Pvalue', 'FDR-Pvalue', 'Bonferroni-Pvalue']



# kruskal_wallis(data)
# G_fit(data, williams=True)




def run_ANOVA(listed_category_data):
    """Compute the ANOVA for inputed lists of data

    Take in the list of array data produced by get_category_arrays.
    Runs ANOVA_one_way on each OTU.

    Returns a list of group means and the propability from ANOVA

        [([0.2, 0.0], 0.40708382206558885),
        ([0.0, 0.25], 0.29235199244023835)]
    """
    result = []
    counter = 0

    while counter < len(listed_category_data[0]):
        a = []
        
        # Create a list of lists where contents are in Numbers format
        a.append(Numbers(listed_category_data[0][counter].tolist()))
        a.append(Numbers(listed_category_data[1][counter].tolist()))
        
        # Compute ANOVA
        dfn, dfd, F, between_MS, within_MS, group_means, prob = ANOVA_one_way(a)
        
        # Append ANOVA output to result
        output = group_means, prob
        result.append(output)
        counter += 1

    return result

def run_correlation(bt_data, category_values, corr_method, tails=None, 
                    permutations=999, confidence_level=0.95):
    """Compute specified correlation between x and y data

    Take in bt_data (an array of the biom table observations), a list of the
    category_values (obtained from the mapping file), and the desired
    correlation method (either 'spearman' or 'pearson')

    Returns list of results, each entry containing corr_coeff, parametric_p_val,
    and nonparametric_p_val
    """
    result = []
    
    #Convert continuous variables from the mapping file to integers
    cat_vals = [int(i) for i in category_values]

    for otu in bt_data:
        # compute correlation
        (corr_coeff, parametric_p_val, permuted_corr_coeffs, \
        nonparametric_p_val, (ci_low, ci_high)) = \
        correlation_test(otu, cat_vals, corr_method, tails, permutations, confidence_level)

        # grab selected correlation info and append to results
        output = corr_coeff, parametric_p_val, nonparametric_p_val
        result.append(output)

    return result







