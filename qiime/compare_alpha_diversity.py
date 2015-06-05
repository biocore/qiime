#!/usr/bin/env python
# File created on 19 May 2011
from __future__ import division

__author__ = "William Van Treuren"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["William Van Treuren", "Greg Caporaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Van Treuren"
__email__ = "vantreur@colorado.edu"

from itertools import combinations, izip
from collections import defaultdict
from numpy import array, isnan, min as np_min
from skbio.draw import boxplots
from qiime.format import format_p_value_for_num_iters
from qiime.parse import (parse_mapping_file_to_dict,
                         parse_rarefaction,
                         group_by_field,
                         parse_mapping_file)
from qiime.stats import (mc_t_two_sample, t_two_sample,
                         benjamini_hochberg_step_down,
                         bonferroni_correction, fdr_correction)
from itertools import combinations
from collections import defaultdict

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
    sids = rarefaction_data[0][3:]  # 0-2 are header strings
    categories = []
    check = 0
    for k in sids:
        try:
            categories.append(mapping_data[k][category])
        except KeyError:  # a sample did not have the category
            check += 1
    if check == len(sids):  # no samples had the category passed
        raise ValueError('No samples had a value for category: \'%s\'.' %
                         category)
    combos = list(combinations(set(categories), 2))
    # make sampleId combos list eg.
    #[(['Id4'],['Id1','Id2']), (['Id1','Id2'],['Id3]), (['Id4'],['Id3])]
    sid_pairs = []
    for pair0, pair1 in combos:
        pair0_sids = [sid for sid in sids if
                      category in mapping_data[sid] and
                      mapping_data[sid][category] == pair0]
        pair1_sids = [sid for sid in sids if
                      category in mapping_data[sid] and
                      mapping_data[sid][category] == pair1]
        sid_pairs.append((pair0_sids, pair1_sids))
    return sid_pairs, combos


def _correct_compare_alpha_results(result, method):
    """Correct compare_alpha_diversities for multiple comps based on method.
    Inputs:
     result - dict, output of compare_alpha_diversities.
     method - str, in ['FDR','Bonferroni','None']
    Notes:
     Nones are ignored for purposes of calculating how many comparisons have
     been done.
    """
    method = method.lower()
    if method not in correction_types:
        raise ValueError('You must specify a method to correct for multiple ' +
                         'comparisons. You may pass \'bonferroni\' or \'fdr\' or \'none\'.')

    corrected_result = {}
    if method == 'bonferroni':
        # correct for comparisons which get tval,pval=None,None
        num_comps = float(len(result)) - result.values().count((None, None))
        for k, v in result.items():
            if v[0] is None:
                corrected_result[k] = v
            else:
                corrected_result[k] = (v[0], min(v[1] * num_comps, 1.0))
    elif method == 'fdr':
        # pull out the uncorrected pvals and apply fdr correction. skip vals
        # which are (None,None). If left in, FDR correction will count them as
        # tests for the purpose of the correction factor.
        tmp_pvals = [v[1] for k, v in result.items() if v != (None, None)]
        fdr_corr_vals = fdr_correction(tmp_pvals)
        # Place values in corrected_result in same order as they were extracted
        # skipping nones.
        i = 0
        for k in result:
            t, p = result[k]
            if (t, p) == (None, None):
                corrected_result[k] = (None, None)
            else:
                corrected_result[k] = (t, min(fdr_corr_vals[i], 1.0))
                i += 1
    elif method == 'none':
        corrected_result = result
    return corrected_result


def get_category_value_to_sample_ids(mapping_lines, category):
    mapping_data, headers, _ = parse_mapping_file(mapping_lines)
    return group_by_field([headers] + mapping_data, category)


def collapse_sample_diversities_by_category_value(category_value_to_sample_ids,
                                                  per_sample_average_diversities):
    result = defaultdict(list)
    for cat, sids in category_value_to_sample_ids.items():
        for sid in sids:
            try:
                sid_average_diversity = per_sample_average_diversities[sid]
            except KeyError:
                pass
            else:
                result[cat].append(sid_average_diversity)
    return result


def get_per_sample_average_diversities(rarefaction_data, depth=None):
    """Extract data rows from rarefaction data matrix.

    Notes: if depth is not given default to the deepest rarefaction available.
    Rarefaction file is not guaranteed to be in order of rarefaction depth so we
    traverse the entire thing looking for the max.

    Inputs:
     rarefaction_data - tuple of lists, results from parse_rarefaction. First
     entry is header line (split on tabs), second is comments, third is file
     names, 4th is list/array of values where first two columns ar depth of
     rarefaction and iteration at that depth (respectively).
     depth - int, depth to use or None if max depth should be used.
    Outputs:
     dict, {sampleID:avg_score_at_given_depth}
    """
    if depth is None:
        depth = array(rarefaction_data[3])[:, 0].max()
    rare_mat = array([row for row in rarefaction_data[3] if row[0] == depth])
    # Average each col of the rarefaction mtx. Computing t test on averages over
    # all iterations. Avoids more comps which kills significance.
    rare_mat = rare_mat.mean(0)[2:]  # remove depth,iter cols
    sids = rarefaction_data[0][3:]  # 0-2 are header strings
    return dict(izip(sids, rare_mat))


def generate_alpha_diversity_boxplots(rarefaction_lines,
                                      mapping_lines,
                                      category,
                                      depth=None):
    rarefaction_data = parse_rarefaction(rarefaction_lines)

    category_value_to_sample_ids = \
        get_category_value_to_sample_ids(mapping_lines,
                                         category)

    per_sample_average_diversities = \
        get_per_sample_average_diversities(rarefaction_data,
                                           depth)

    per_category_value_average_diversities = \
        collapse_sample_diversities_by_category_value(
            category_value_to_sample_ids,
            per_sample_average_diversities)

    # sort the data alphabetically
    sorted_per_category_value_average_diversities = \
        sorted(per_category_value_average_diversities.items())

    x_tick_labels = []
    distributions = []
    for cat, avg_diversities in sorted_per_category_value_average_diversities:
        x_tick_labels.append("%s (n=%d)" % (cat, len(avg_diversities)))
        distributions.append(avg_diversities)

    return boxplots(distributions, x_tick_labels=x_tick_labels)


def compare_alpha_diversities(rarefaction_lines, mapping_lines, category,
                              depth=None, test_type='nonparametric', num_permutations=999):
    """Compares alpha diversity values for differences per category treatment.

    Notes:
     Returns a defaultdict which as keys has the pairs of treatments being
     compared, and as values, lists of (pval,tval) tuples for each comparison at
     for a given iteration.
    Inputs:
     rarefaction_lines - list of lines, result of multiple rarefactions.
     mapping_lines - list of lines, mapping file lines.
     category - str, the category to be compared, eg 'Treatment' or 'Age'.
     depth - int, depth of the rarefaction file to use. if None, then will use
     the deepest available in the file.
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

    ps_avg_div = get_per_sample_average_diversities(rarefaction_data, depth)

    ttest_results, ad_avgs = {}, {}
    for sid_pair, treatment_pair in zip(samid_pairs, treatment_pairs):
        # if there is only 1 sample for each treatment in a comparison, and mc
        # using mc method, will error (e.g. mc_t_two_sample([1],[1]).
        if len(sid_pair[0]) == 1 and len(sid_pair[1]) == 1:
            ttest_results[treatment_pair] = (None, None)
            # add alpha diversity averages and standard deviations. since their
            # is only a single sample if we are in this part of the loop, we can
            # just record the sample value as the avg and 0 as the std.
            ad_avgs[treatment_pair[0]] = (sid_pair[0][0], 0.)
            ad_avgs[treatment_pair[1]] = (sid_pair[1][0], 0.)
        else:
            i = array([ps_avg_div[x] for x in sid_pair[0]])
            j = array([ps_avg_div[x] for x in sid_pair[1]])
            # add alpha diversity averages and standard deviations.
            ad_avgs[treatment_pair[0]] = (i.mean(), i.std())
            ad_avgs[treatment_pair[1]] = (j.mean(), j.std())
            # conduct tests
            if isnan(np_min(i)) or isnan(np_min(j)):
                ttest_results[treatment_pair] = (None, None)
                continue
            if test_type == 'parametric':
                obs_t, p_val = t_two_sample(i, j)
            elif test_type == 'nonparametric':
                obs_t, _, _, p_val = mc_t_two_sample(i, j,
                                                     permutations=num_permutations)
                if p_val is not None:
                    p_val = float(format_p_value_for_num_iters(p_val,
                                                               num_iters=num_permutations))
                elif p_val is None:  # None will error in format_p_val
                    obs_t, p_val = None, None
            else:
                raise ValueError("Invalid test type '%s'." % test_type)
            ttest_results[treatment_pair] = (obs_t, p_val)

    return ttest_results, ad_avgs
