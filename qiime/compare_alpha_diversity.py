#!/usr/bin/env python
# File created on 19 May 2011
from __future__ import division

__author__ = "William Van Treuren"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["William Van Treuren", "Greg Caporaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.5.0-dev"
__maintainer__ = "William Van Treuren"
__email__ = "vantreur@colorado.edu"
__status__ = "Development"

from numpy import array, isnan, min as np_min
from qiime.format import format_p_value_for_num_iters
from qiime.parse import parse_mapping_file_to_dict, parse_rarefaction
from cogent.maths.stats.test import mc_t_two_sample, t_two_sample
from itertools import combinations
from collections import defaultdict
from qiime.otu_category_significance import fdr_correction

test_types = ['parametric', 'nonparametric']
correction_types = ['bonferroni', 'fdr', 'none']

def sampleId_pairs(mapping_data, rarefaction_data, category):
    """Returns list of sampleId tuples.
    Notes:
     The function first constructs pairwise combinations of values found under
     the passed category. It then finds which sampleId's belong to each
     category and returns a list of tuples where each tuple is a treatment pair
     and the entries in the tuples are lists of the sampleId's which have one of
     those treatments.
    Inputs:
     mapping_data - nested dict, output of parse_mapping_file_to_dict.
     rarefaction_data - list of 4 elements, output of parse_rarefaction.
     category - str, category that samples have in the mapping file. eg 'Weight'
    """
    # make treatment pairs list eg. 
    # [('Obese','Fat'),('Obese','notFat'),('Fat','notFat')]


    # frequently the user will have a filtered their table before rarefactions 
    # such that the mapping file contains sampleId's that are not found in the 
    # filtered table. this caused an error with the script because it wasn't 
    # able to find the id's. to avoid this, make treatment covering set from 
    # only the sampleIds found in the rarefaction file. 
    sids = rarefaction_data[0][3:] # 0-2 are header strings
    categories = []
    check = 0
    for k in sids:
        try:
            categories.append(mapping_data[k][category])
        except KeyError: #a sample did not have the category
            check+=1
    if check==len(sids): #no samples had the category passed
        raise ValueError('No samples had a value for category: \'%s\'.' % 
            category)
    combos = list(combinations(set(categories),2))
    # make sampleId combos list eg. 
    #[(['Id4'],['Id1','Id2']), (['Id1','Id2'],['Id3]), (['Id4'],['Id3])]
    sid_pairs = []
    for pair0, pair1 in combos:
        pair0_sids = [sid for sid in sids if \
            mapping_data[sid].has_key(category) and \
            mapping_data[sid][category] == pair0]
        pair1_sids = [sid for sid in sids if \
            mapping_data[sid].has_key(category) and \
            mapping_data[sid][category] == pair1]
        sid_pairs.append((pair0_sids,pair1_sids))
    return sid_pairs, combos

def _correct_compare_alpha_results(result, method):
    """Correct compare_alpha_diversities for multiple comps based on method.
    Inputs:
     result - dict, output of compare_alpha_diversities.
     method - str, in ['FDR','Bonferroni','None']
    """
    method = method.lower()
    if method not in correction_types:
        raise ValueError('You must specify a method to correct for multiple '+\
            'comparisons. You may pass \'bonferroni\' or \'fdr\' or \'none\'.')

    corrected_result = {}
    if method == 'bonferroni':
        num_comps = float(len(result))
        for k,v in result.items():
            corrected_result[k] = (v[0],min(v[1]*num_comps,1.0))
    elif method == 'fdr':
        # pull out the uncorrected pvals and apply fdr correction
        tmp_pvals = [v[1] for k,v in result.items()]
        fdr_corr_vals = fdr_correction(tmp_pvals)
        # create corrected results by going through result in same order as 
        # pvalues were removed and replacing them with fdr corrected
        for i,k in enumerate(result): #steps through in same order as items 
            t,p = result[k]
            corrected_result[k] = (t, min(fdr_corr_vals[i],1.0)) #same as above
    elif method == 'none':
        corrected_result = result
    return corrected_result

def compare_alpha_diversities(rarefaction_lines, mapping_lines, category, depth,
    test_type='nonparametric', num_permutations=999):
    """Compares alpha diversity values for differences per category treatment.
    Notes: 
     Returns a defaultdict which as keys has the pairs of treatments being 
     compared, and as values, lists of (pval,tval) tuples for each comparison at
     for a given iteration.     
    Inputs:
     rarefaction_lines - list of lines, result of multiple rarefactions.
     mapping_lines - list of lines, mapping file lines. 
     category - str, the category to be compared, eg 'Treatment' or 'Age'.
     depth - int, depth of the rarefaction file to use.
     test_type - str, the type of t-test to perform. Must be either
     'parametric' or 'nonparametric'.
     num_permutations - int, the number of Monte Carlo permutations to use if
     test_type is 'nonparametric'.    
    """
    if test_type == 'nonparametric' and num_permutations < 1:
        raise ValueError("Invalid number of permutations: %d. Must be greater "
                         "than zero." % num_permutations)
     
    rarefaction_data = parse_rarefaction(rarefaction_lines)
    mapping_data = parse_mapping_file_to_dict(mapping_lines)[0]
    # samid_pairs, treatment_pairs are in the same order
    samid_pairs, treatment_pairs = sampleId_pairs(mapping_data, 
        rarefaction_data, category)
    
    # extract only rows of the rarefaction data that are at the given depth
    rare_mat = array([row for row in rarefaction_data[3] if row[0]==depth])
    # average each column of the rarefaction matrix because we are computing 
    # the t_test on the average scores for a given comparison. this avoids 
    # a much larger number of tests which kills significance.
    rare_mat = rare_mat.sum(0)/rare_mat.shape[0]
    sids = rarefaction_data[0][3:] # 0-2 are header strings
    results = {}

    for sid_pair, treatment_pair in zip(samid_pairs, treatment_pairs):
        # first two cols of rare_mat are depth, iteration, don't want to grab 
        # those for values so we add 2 to all indices.
        pair0_indices = [sids.index(i)+2 for i in sid_pair[0]]
        pair1_indices = [sids.index(i)+2 for i in sid_pair[1]]
        t_key = '%s,%s' % (treatment_pair[0], treatment_pair[1])
        i = rare_mat.take(pair0_indices)
        j = rare_mat.take(pair1_indices)
        # found discussion of how to quickly check an array for nan here:
        # http://stackoverflow.com/questions/6736590/fast-check-for-nan-in-numpy
        if isnan(np_min(i)) or isnan(np_min(j)):
            raise ValueError("nan present in alpha diversities at depth %d."
             " Re-run alpha diversity comparison at lower sampling depth where"
             " diversity was computed for all samples, or recompute alpha"
             " diversity after filtering samples with too few sequences." % depth)
        if test_type == 'parametric':
            obs_t, p_val = t_two_sample(i,j)
        elif test_type == 'nonparametric':
            obs_t, _, _, p_val = mc_t_two_sample(i,j, 
                permutations=num_permutations)
            p_val = float(format_p_value_for_num_iters(p_val, 
                num_iters=num_permutations))
        else:
            raise ValueError("Invalid test type '%s'." % test_type)
        results[t_key]= (obs_t,p_val)
    return results
