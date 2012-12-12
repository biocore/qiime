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

from numpy import array, isnan
from qiime.format import format_p_value_for_num_iters
from qiime.parse import parse_mapping_file_to_dict, parse_rarefaction
from cogent.maths.stats.test import mc_t_two_sample, t_two_sample
from itertools import combinations
from collections import defaultdict
from qiime.otu_category_significance import fdr_correction

test_types = ['parametric', 'nonparametric']

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
    print sids, mapping_data.keys()
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
    print sid_pairs, combos
    return sid_pairs, combos

def _correct_compare_alpha_results(result, method):
    """Correct compare_alpha_diversities for multiple comps based on method.
    Inputs:
     result - dict, output of compare_alpha_diversities.
     method - str, in ['FDR','Bonferroni','None']
    """
    if method not in ['Bonferroni','FDR','None']:
        raise ValueError('You must specify a method to correct for multiple '+\
            'comparisons. You may pass \'Bonferroni\' or \'FDR\' or \'None\'.')
    
    if method == 'Bonferroni':
        # the number of comparisons is the number of iterations times the number
        # of different groups compared. keys=groups compared, and the length
        # of any item from any key will be the same, is the number of iters
        num_comps = float(len(result.keys())*len(result.items()[0][1]))
        corrected_result = {}
        for k,v in result.items():
            corrected_result[k] = [(tval,pval*num_comps) for tval,pval in v]
        # this returns bizarre results like pval=35 since the bonferroni isn't
        # really designed to correct individual values, but rather to correct
        # the level of your test. 
    elif method == 'FDR':
        # pull out the uncorrected pvals and apply fdr correction
        tmp_pvals = []
        [tmp_pvals.extend([k[1] for k in j]) for i,j in result.items()]
        fdr_corr_vals = fdr_correction(tmp_pvals)
        # create corrected results by going through result in same order as 
        # pvalues were removed and replacing them with fdr corrected
        corrected_result = defaultdict(list)
        count = 0
        for k in result:
            for t,p in result[k]: #t=tval, p=pval
                corrected_result[k].append((t,fdr_corr_vals[count]))
                count+=1
    elif method == 'None':
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
    sids = rarefaction_data[0][3:] # 0-2 are header strings
    results = defaultdict(list)

    for sid_pair, treatment_pair in zip(samid_pairs, treatment_pairs):
        # first two cols of rare_mat are depth, iteration, don't want to grab 
        # those for values so we add 2 to all indices.
        pair0_indices = [sids.index(i)+2 for i in sid_pair[0]]
        pair1_indices = [sids.index(i)+2 for i in sid_pair[1]]
        t_key = '%s,%s' % (treatment_pair[0], treatment_pair[1])
        for row in rare_mat: 
            i = row.take(pair0_indices)
            j = row.take(pair1_indices)
            if test_type == 'parametric':
                obs_t, p_val = t_two_sample(i,j)
            elif test_type == 'nonparametric':
                obs_t, _, _, p_val = mc_t_two_sample(i,j, 
                    permutations=num_permutations)
                #format_p_value returns a string for some reason
                p_val = float(format_p_value_for_num_iters(p_val, 
                    num_iters=num_permutations))
            else:
                raise ValueError("Invalid test type '%s'." % test_type)
            results[t_key].append((obs_t,p_val))
    return results


