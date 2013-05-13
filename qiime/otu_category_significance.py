#!/usr/bin/env python
# File created on 07 Oct 2009.
from __future__ import division

__author__ = "Catherine Lozupone"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Catherine Lozupone", "Jesse Stombaugh", "Dan Knights",
               "Jai Ram Rideout", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.6.0-dev"
__maintainer__ = "Catherine Lozupone"
__email__ = "lozupone@colorado.edu"
__status__ = "Development"


from optparse import OptionParser
from os.path import split, splitext, join
from os import listdir
from string import strip
import sys
import numpy as np
from numpy import argsort, mean, array
from cogent.util.dict2d import Dict2D
from qiime.pycogent_backports.test import calc_contingency_expected, \
    G_fit_from_Dict2D, ANOVA_one_way, correlation, t_paired
from cogent.maths.stats.util import Numbers
from qiime.longitudinal_otu_category_significance import get_sample_individual_info
from qiime.parse import parse_mapping_file
from biom.exception import TableException

"""Look for OTUs that are associated with a category. Currently can do:
    1) perform g-test of independence to determine whether OTU presence
    absence is associated with a category in the mapping file. Can
    be used to look for co-occurrence using qPCR data or to look for 
    associations with a categorical environmental variable.
    2) perform ANOVA to determine whether OTU abundance is significantly
    different across a category
    3) perform a pearson correlation to determine whether OTU abundance is
    associated with a continuous variable in the category mapping file (e.g. pH)
"""

def filter_OTUs(otu_table, filt, all_samples=True,\
                category_mapping_info=None, max_filter=None):
    """Get the list of OTUs found in at least <filter> samples. 

    The filter is supplied as a fraction of samples.

    optionally filters out those that are found in all_samples if True:                                                              
        G_test: set to True since can't show presence/absence patterns if                                                            
            present in all (also causes an error)                                                                                    
        ANOVA: set to False since can still see differences in abundance                                                             

    optionally takes a category mapping file as input and only considers                                                             
        inclusion in samples that are included therein. This is a method                                                             
        of making sure that the samples represented in the category                                                                  
        mapping file and the OTU table are in sync         

    If max_filter is set, it will exclude OTUs from the list that are 
        found above that value.
    """
    result = []
    # Only consider samples that are present in mapping file AND otu table.
    if category_mapping_info:
        mapping_file_samples = set(category_mapping_info.keys())
        otu_table_samples = set(otu_table.SampleIds)
        included_samples = mapping_file_samples & otu_table_samples
        otu_table = otu_table.filterSamples(lambda v,i,m: i in included_samples)

    # min and max number of samples per observation to keep
    min_filter = int(round(filt * len(otu_table.SampleIds)))
    if max_filter:
        max_filter = int(round(max_filter * len(otu_table.SampleIds)))

    def filter_f(values, id_, md):
        """filter observations based on how many samples are covered"""
        # determine the number of samples the observation occurs in
        n_samples = sum(values > 0)

        # we don't have enough samples represented
        if n_samples < min_filter:
            return False

        # we have to many samples represented
        if max_filter and n_samples > max_filter:
            return False

        # if all samples have this observation
        if all_samples and n_samples == len(values):
            return False

        return True
    
    try:
        filtered = otu_table.filterObservations(filter_f)
    except TableException:
        # all observations filtered out
        return []

    return list(filtered.ObservationIds)

def sync_mapping_to_otu_table(otu_table, mapping):
    """removes samples from the mapping file that are not in the otu table

    otu_table and mapping are each parsed with the standard parsers first
    returns a new parsed mapping and a list of samples removed
    """
    mapping_data, header, comments = mapping
    otu_table_samples = otu_table.SampleIds

    new_mapping_data = []
    removed_samples = []
    for i in mapping_data:
        sample = i[0]
        if sample not in otu_table_samples:
            removed_samples.append(sample)
        else:
            new_mapping_data.append(i)
    return [new_mapping_data, header, comments], removed_samples

def run_single_G_test(OTU_name, category_info, otu_table, category_values,
                      suppress_warnings=True):
    """run the G test on a single OTU
       If suppress_warnings=True, doesn't warn when sample in map but not otu table
    """
    contingency_matrix = make_contingency_matrix(OTU_name, category_info, \
        otu_table, category_values, suppress_warnings=suppress_warnings)
    contingency_matrix = calc_contingency_expected(contingency_matrix)
    g_val, prob = G_fit_from_Dict2D(contingency_matrix)
    return g_val, prob, contingency_matrix

def make_contingency_matrix(OTU_name, category_info, otu_table,
                            category_values, suppress_warnings=True):
    """make the contingency table for running the G test of independence
    
    counts OTU as present (count > 1) or absent
    makes a column in the matrix for each category value
    If suppress_warnings=True, doesn't warn when sample in map but not otu table
    """
    result = {'OTU_pos':{},
                'OTU_neg':{}}
    for category in category_values:
        result['OTU_pos'][category + '_pos'] = 0
        result['OTU_neg'][category + '_pos'] = 0
    for sample in category_info:
        category = category_info[sample]
        try:
            OTU_count = float(otu_table.getValueByIds(OTU_name, sample))
            OTU_count = int(OTU_count)
            worked = True
        except Exception: 
            if not suppress_warnings:
                print "Warning:",sample,"is in the  sample mapping file but not the OTU table" 
            worked = False
        if worked:
            if OTU_count == 0:
                result['OTU_neg'][category + '_pos'] += 1
            elif OTU_count > 0: 
                result['OTU_pos'][category + '_pos'] += 1
    return Dict2D(result, Default=0, Pad=True)

def run_G_test_OTUs(OTU_list, category_info, otu_table, category_values,
                    suppress_warnings=True):
    """run the G test on all OTUs in the OTU list
       If suppress_warnings=True, doesn't warn when sample in map but not otu table
    """
    result = {}
    num_comparisons = len(OTU_list)
    for OTU in OTU_list:
        g_val, g_prob, contingency_matrix = run_single_G_test(OTU,
             category_info, otu_table, category_values,
             suppress_warnings=suppress_warnings)
        Bonf_p_val = g_prob * num_comparisons
        result[OTU] = [g_val, g_prob, contingency_matrix, Bonf_p_val]
    return result

def run_ANOVA_OTUs(OTU_list, category_info, otu_table, category_values):
    """run ANOVA on all OTUs in the OTU list"""
    result = {}
    num_comparisons = len(OTU_list)
    for OTU in OTU_list:
        if (OTU in otu_table.ObservationIds and
            sum(map(float, otu_table.observationData(OTU))) > 0):
            group_means, prob = run_single_ANOVA(OTU, category_info, \
                otu_table, category_values)
        else:
            group_means = [0.0] * len(category_values)
            prob = 0.5
        Bonf_p_val = prob * num_comparisons
        result[OTU] = [group_means, prob, Bonf_p_val]
    return result

def run_single_ANOVA(OTU, category_info, otu_table, category_values):
    """runs ANOVA on the designated OTU"""
    result = {}
    #get a list of values for each category
    values = []
    for category in category_values:
        values.append(Numbers([]))
    sample_data = otu_table.observationData(OTU)
    for sample in category_info:
        if sample in otu_table.SampleIds:
            sample_index = otu_table.SampleIds.index(sample)
            count = sample_data[sample_index]
            category = category_info[sample]
            index = category_values.index(category)
            values[index].append(count)
    #    else:
    #        print "Warning " + sample + "is in the category mapping file " +\
    #            "but not the OTU table"
    try:
        dfn, dfd, F, between_MS, within_MS, group_means, prob = ANOVA_one_way(values)
        return group_means, prob
    except ValueError:
        #set the p-value to 'diff' if the variances are 0.0 (within rounding 
        #error) and the means are not all the same. If the means are all
        #the same and the variances are 0.0, set the p-value to 1
        group_means = []
        group_variances = []
        for i in values:
            group_means.append(i.Mean)
            group_variances.append(i.Variance)
        group_means = set(group_means)
        if sum(group_variances) < 1e-21 and len(group_means) > 1:
            prob = 0.0
        else:
            prob = 1.0
        return group_means, prob

def run_paired_T_test_OTUs(OTU_list, mapping_data, header, individual_column,\
    timepoint_zero_column, otu_table, ignore_val, filt):
    """
    Runs a paired T test on all of the OTUs in the OTU list
    """
    all_results = {}
    for OTU in OTU_list:
        result = run_single_paired_T_test(OTU, mapping_data, header,
            individual_column, timepoint_zero_column, otu_table, ignore_val,
            filt)
        if result:
            t, prob, after_treatment_vals, num_pairs = result
            average_diff = sum(after_treatment_vals)/float(len(after_treatment_vals))

            all_results[OTU] = [t, prob, average_diff, num_pairs]
    return all_results

def run_single_paired_T_test(OTU, mapping_data, header, individual_column,\
    timepoint_zero_column, otu_table, ignore_val, filt):
    """run the paired T test on a single OTU
    """
    timepoint_zero_vals, after_treatment_vals = get_single_paired_T_values(OTU,
        mapping_data, header, individual_column, timepoint_zero_column,
        otu_table, ignore_val)
    #get the number of samples from the category mapping and multiply by the filter
    #two samples per individual so divide by 2
    filt = (filt*len(mapping_data))/2
    if len(timepoint_zero_vals) >= round(filt):
        t, prob = t_paired(timepoint_zero_vals, after_treatment_vals,
                           tails=None, exp_diff=0)
        return t, prob, after_treatment_vals, len(after_treatment_vals)
    else:
        return None

def get_single_paired_T_values(OTU, mapping_data, header, individual_column,
    timepoint_zero_column, otu_table, ignore_val):
    """gets the values for running a paired T test
   
    the timepoint_zero column should all be zero. The other column from that
    individual should be either negative or positive depending. If the
    value is the ignore value, it should be ignored. If the number of values
    do not pass the filter, the OTU should be ignored
    """
    timepoint_zero_vals = []
    after_treatment_vals = []
    #sample_to_subtract is the sample_names mapped to the timepoint zero sample
    #names
    samples_from_subject, sample_to_subtract = get_sample_individual_info(\
        mapping_data, header, individual_column, timepoint_zero_column)
    for sample in samples_from_subject:
        num_samples = len(samples_from_subject[sample])
        if num_samples != 2:
            raise ValueError("Each individual/site must contain exactly 2 samples")
    timepoint_zeros = []
    for sample in sample_to_subtract:
        #if the sample is mapped to itself it is timepoint zero
        if sample == sample_to_subtract[sample]:
            timepoint_zeros.append(sample)
    for sample in timepoint_zeros:
        if sample in otu_table.SampleIds:
            timepoint_zero_val = otu_table.getValueByIds(OTU, sample)
        else:
            timepoint_zero_val = None
        if timepoint_zero_val != ignore_val:
            for sample2 in sample_to_subtract:
                if sample2 != sample and sample_to_subtract[sample2] == sample:
                    after_treatment_sample_id = sample2
                    if after_treatment_sample_id in otu_table.SampleIds:
                        if timepoint_zero_val != None:
                            after_treatment_vals.append(
                                    otu_table.getValueByIds(
                                        OTU, after_treatment_sample_id))
                            #print "OTU: %s" % OTU
                            #print "sample id: %s" % after_treatment_sample_id
                            #print "value: %s" % otu_table.getValueByIds(OTU, after_treatment_sample_id)
                            timepoint_zero_vals.append(timepoint_zero_val)
    return timepoint_zero_vals, after_treatment_vals

def run_correlation_OTUs(OTU_list, category_info, otu_table, ignore_val=None,
                         filt=1):
    """runs pearson correlation on all OTUs in the OTU list
    """
    result = {}
    num_samples = len(category_info)
    filt = filt * num_samples
    for OTU in OTU_list:
        if OTU in otu_table.ObservationIds:
            OTU_abundance_vals, category_vals = \
                get_single_correlation_values(OTU, category_info, \
                otu_table, ignore_val=ignore_val)
            if len(category_vals) >= round(filt):
                r, prob = run_single_correlation(OTU_abundance_vals, \
                    category_vals)
            else:
                r, prob = None, None
        else:
            r = 0.0
            prob = 0.5
            OTU_abundance_vals = 'NA'
            category_vals = 'NA'
        if r:
            result[OTU] = [r, prob, str(OTU_abundance_vals), str(category_vals)]
    return result

def run_single_correlation(OTU_abundance_values, category_values):
    """runs pearson correlation  on the designated OTU
    """
    return correlation(Numbers(category_values), Numbers(OTU_abundance_values))

def get_single_correlation_values(OTU, category_info, otu_table,
                                  ignore_val=None):
    """gets the x and y values for running pearson correlation
    
    ignore_val: A value that should be ignored in the calculations. This
        was made for a special application (Ley humans for obesity metaanalysis
        where I wanted to ignore individuals in which the otu was never
        observed across a time series. If the OTU value is the ignore_val, it
        does not include that sample in the regression.
    """
    #get a list of values for each category
    OTU_abundance_values = []
    category_values = []
    for sample in category_info:
        # even if this OTU is not observed, we can use count=0
        if sample in otu_table.SampleIds:
            count = float(otu_table.getValueByIds(OTU, sample))
        #else:
        #    count = 0
            add = True
            if ignore_val:
                if count == ignore_val:
                    add=False
            if add:
                try:
                    cat_val = float(category_info[sample])
                    category_values.append(cat_val)
                    OTU_abundance_values.append(count)
                except ValueError:
                    raise ValueError("The category values must be numeric to use the correlation option")
    return OTU_abundance_values, category_values

def output_results_G_test(G_test_results, taxonomy_info=None):
    """creates the results output using result of run_G_test_OTUs"""
    header = ['OTU', 'g_val', 'g_prob', 'Bonferroni_corrected', 'FDR_corrected']
    sample_cont_matrix = G_test_results[G_test_results.keys()[0]][2]
    for i in sample_cont_matrix:
        for j in sample_cont_matrix[i]:
            header.append(i + '##' + j)
    if taxonomy_info:
            header.append('Consensus Lineage')
    output = ['\t'.join(header)]
    #perform fdr correction
    G_test_results = add_fdr_correction_to_results(G_test_results)
    for OTU in G_test_results:
        g_val = str(G_test_results[OTU][0])
        g_prob = str(G_test_results[OTU][1])
        Bonf_p = str(G_test_results[OTU][3])
        fdr_p = str(G_test_results[OTU][4])
        line = [OTU, g_val, g_prob, Bonf_p, fdr_p]
        for i in G_test_results[OTU][2]:
            for j in G_test_results[OTU][2][i]:
                line.append(str(G_test_results[OTU][2][i][j]))
        if taxonomy_info:
            line.append(taxonomy_info[OTU])
        output.append('\t'.join(line))
    return sort_rows(output, 2)

def output_results_ANOVA(ANOVA_results, category_values, taxonomy_info=None):
    """creates the results output using result of run_ANOVA_OTUs"""
    header = ['OTU', 'prob', 'Bonferroni_corrected', 'FDR_corrected']
    for i in category_values:
        header.append(i + '_mean')
    if taxonomy_info:
            header.append('Consensus Lineage')
    output = ['\t'.join(header)]
    #perform fdr correction
    ANOVA_results = add_fdr_correction_to_results(ANOVA_results)
    for OTU in ANOVA_results:
        prob = str(ANOVA_results[OTU][1])
        Bonf_p = str(ANOVA_results[OTU][2])
        fdr_p = str(ANOVA_results[OTU][3])
        line = [OTU, prob, Bonf_p, fdr_p]
        for i in ANOVA_results[OTU][0]:
            line.append(str(i))
        if taxonomy_info:
            line.append(taxonomy_info[OTU])
        output.append('\t'.join(line))
    return sort_rows(output, 1)

def output_results_correlation(correlation_results, taxonomy_info=None):
    """creates the results output using result of run_correlation_OTUs"""
    header = ['OTU', 'prob', 'otu_values_y', 'cat_values_x', 'Bonferroni_corrected', 'FDR_corrected', 'r']
    if taxonomy_info:
            header.append('Consensus Lineage')
    output = ['\t'.join(header)]
    correlation_results = add_bonferroni_to_results(correlation_results)
    correlation_results = add_fdr_correction_to_results(correlation_results)
    for OTU in correlation_results:
        prob = str(correlation_results[OTU][1])
        otu_values = str(correlation_results[OTU][2]) 
        cat_values = str(correlation_results[OTU][3]) 
        Bonf_p = str(correlation_results[OTU][4])
        fdr_p = str(correlation_results[OTU][5])
        r = str(correlation_results[OTU][0])
        line = [OTU, prob, otu_values, cat_values, Bonf_p, fdr_p, r]
        if taxonomy_info:
            line.append(taxonomy_info[OTU])
        output.append('\t'.join(line))
    return sort_rows(output, 1)

def output_results_paired_T_test(paired_T_results, taxonomy_info=None):
    """creates the results output using result of run_paired_T_test_OTUs
    """
    header = ['OTU', 'prob', 'T stat', 'average_diff', 'num_pairs', \
        'Bonferroni_corrected', 'FDR_corrected']
    if taxonomy_info:
            header.append('Consensus Lineage')
    output = ['\t'.join(header)]
    paired_t_results = add_bonferroni_to_results(paired_T_results)
    paired_t_results = add_fdr_correction_to_results(paired_T_results)
    for OTU in paired_T_results:
        prob = str(paired_T_results[OTU][1])
        T_stat = str(paired_T_results[OTU][0])
        average_diff = str(paired_T_results[OTU][2])
        num_pairs = str(paired_T_results[OTU][3])
        Bonf_p = str(paired_T_results[OTU][4])
        fdr_p = str(paired_T_results[OTU][5])
        line = [OTU, prob, T_stat, average_diff, num_pairs, Bonf_p, fdr_p]
        if taxonomy_info:
            line.append(taxonomy_info[OTU])
        output.append('\t'.join(line))
    return sort_rows(output, 1)

def sort_rows(output, p_val_index):
    """sorts the rows in the output by p-value
    """
    header = output[0]
    lines = output[1:]
    probs = []
    for line in lines:
        line = line.split('\t')
        p_value = line[p_val_index]
        if p_value != 'None':
            probs.append(float(line[p_val_index]))
        else:
            probs.append(1.1)
    probs = array(probs)
    x = probs.argsort()
    new_lines = [header]
    for i in x:
        new_lines.append(lines[i])
    return new_lines

def add_bonferroni_to_results(results):
    """corrects results using the false discovery rate method.

    ranks the p-values from low to high. multiplies each p-value by the #
    of comparisons divided by the rank. appends the result for each OTU to the
    supplied results dictionary.
    """
    names = []
    probs = []
    for i in results:
        names.append(i)
        probs.append(results[i][1])
    corrected_probs = []
    number_comparisons = len(names)
    for i in probs:
        if i:
            corrected_probs.append(i * number_comparisons)
        else:
            corrected_probs.append('NA')
    for index, prob in enumerate(corrected_probs):
        results[names[index]].append(prob)
    return results

def add_fdr_correction_to_results(results):
    """corrects results using the false discovery rate method.

    ranks the p-values from low to high. multiplies each p-value by the #
    of comparisons divided by the rank. appends the result for each OTU to the
    supplied results dictionary.
    """
    names = []
    probs = []
    for i in results:
        names.append(i)
        probs.append(results[i][1])
    corrected_probs = fdr_correction(probs)
    for index, prob in enumerate(corrected_probs):
        results[names[index]].append(prob)
    return results

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

def get_taxonomy_info(otu_table):
    """Returns a dict mapping OTU ids to taxonomy (if they exist)."""
    taxonomy_info = {}
    if (otu_table.ObservationMetadata is not None and
        otu_table.ObservationMetadata[0]['taxonomy'] is not None):
        for obs_id, obs_metadata in zip(otu_table.ObservationIds,
                                        otu_table.ObservationMetadata):
            curr_tax = obs_metadata['taxonomy']
            taxonomy_info[obs_id] = '; '.join(curr_tax)
    return taxonomy_info

def get_category_info(mapping_data, header, category, threshold=None):
    """parse category mapping file
    
    returns a dict mapping the sample name to the value of category
    
    also returns a list of category values

    If a numerical threshold is specified, converts to 1 or 0 depending on 
    whether the value is below that threshold.
    
    categories can have categorical or continuous data"""
    result = {}
    if threshold and threshold != 'None':
        category_values = ['0', '1']
    else:
        category_values = []
    
    category_index = header.index(category)
    for line in mapping_data:
        sample_name = line[0]
        if threshold and threshold != 'None':
            val = line[category_index]
            val = float(val)
            if val > threshold:
                result[sample_name] = '1'
            else:
                result[sample_name] = '0'
        else:
            category_val = line[category_index]
            result[sample_name] = category_val
            if category_val not in category_values:
                category_values.append(category_val)
    return result, category_values

def aggregate_multiple_results_ANOVA(all_results):
    """Finds average result values for each OTU given multiple results.
       
       all_results: is a dict of {'otuid':[results1,results2,...],...}
    """
    aggregate = {}
    for k in all_results.keys():
        # get within-category means, pvals, bonferroni
        means, pvals, bf_pvals = zip(*all_results[k])
        # get mean of each type of result
        means = array(means).mean(0)
        pval = np.mean(array(pvals))
        bf_pval = np.mean(array(bf_pvals))
        aggregate[k] = [means, pval, bf_pval]
    return aggregate

def aggregate_multiple_results_correlation(all_results):
    """Finds average result values for each OTU given multiple results.
       
       all_results: is a dict of {'otuid':[results1,results2,...],...}
    """
    aggregate = {}
    for k in all_results.keys():
        # get r-squared, pvals, bonferroni
        x = zip(*all_results[k])
        rs = x[0]
        pvals = x[1]
        otu_vals = x[2]
        cat_vals = x[3]
        # get mean of each type of result
        r = np.mean(array(rs))
        pval = np.mean(array(pvals))
        otu_vals = 'NA'
        cat_vals = 'NA'
        aggregate[k] = [r, pval, 'NA', 'NA']
    return aggregate

def aggregate_multiple_results_G_test(all_results):
    """Finds average result values for each OTU given multiple results.
       
       all_results: is a dict of {'otuid':[results1,results2,...],...}
       Note: The first contingency matrix is included instead of the average
    """

    aggregate = {}
    for k in all_results.keys():
        # get gvals, pvals, contingency matrices, bonferroni
        gs, pvals, cms, bf_pvals = zip(*all_results[k])
        # get mean of each type of result
        g = np.mean(array(gs))
        pval = np.mean(array(pvals))
        bf_pval = np.mean(array(bf_pvals))
        
        # get average of each position in contingency matrix
        cm = Dict2D(cms[0])
        for i in cm:
            for j in cm[i]:
                cm[i][j] = mean(array([cm_k[i][j] for cm_k in cms]))
        aggregate[k] = [g, pval, cm, bf_pval]
    return aggregate
    

def test_wrapper(test, otu_table, category_mapping, category, threshold, \
        filt, otu_include=None, ignore_val=None, \
        otu_table_relative_abundance=False, individual_column='individual',\
        timepoint_zero_column='timepoint_zero'):
    """runs statistical test to look for category/OTU associations"""
    if ignore_val == 'None':
        ignore_val = None
    
    if test == 'ANOVA' or test == 'correlation': 
        if not otu_table_relative_abundance:
            otu_table = otu_table.normObservationBySample()
        all_samples = False
    elif test == 'g_test':
        all_samples = True
    elif test == 'paired_T':
        pass
    else:
        raise ValueError("An invalid test statistic was given. (-s option). Valid values are ANOVA, correlation, g_test, paired_T.")

    taxonomy_info = get_taxonomy_info(otu_table)
    mapping_data, header, comments = category_mapping

    if not test == 'paired_T':
        category_info, category_values = \
            get_category_info(mapping_data, header, category, threshold)
    #do not apply the filter_OTUs to the longitudinal studies as they are
    #filtered later
    if test == 'ANOVA' or test == 'correlation' or test == 'g_test':
        OTU_list = filter_OTUs(otu_table, filt, all_samples=all_samples,
                               category_mapping_info=category_info)
    elif test == 'paired_T':
        OTU_list = otu_table.ObservationIds

    #filter OTU_list with the otu_include list
    if otu_include:
        otu_include = [line.strip() for line in otu_include]
        OTU_list = [OTU for OTU in OTU_list if OTU in otu_include]
    if len(OTU_list) == 0:
        raise ValueError("No OTUs remain after applying the filter. Try lowering the filter value (-f option)")
    if test == 'ANOVA':
        results = run_ANOVA_OTUs(OTU_list, category_info, otu_table,
                                 category_values)
        output = output_results_ANOVA(results, category_values, taxonomy_info)
    elif test == 'correlation':
        results = run_correlation_OTUs(OTU_list, category_info, \
            otu_table, ignore_val=ignore_val, filt=filt)
        output = output_results_correlation(results, taxonomy_info)
    elif test == 'g_test':
        results = run_G_test_OTUs(OTU_list, category_info, otu_table, \
                         category_values)
        output = output_results_G_test(results, taxonomy_info)
    elif test == 'paired_T':
        #category info in this case should be the timepoint_zero column.
        #The timepoint_zero column should be used as the category in the wrapper
        results = run_paired_T_test_OTUs(OTU_list, mapping_data, header,
            individual_column, timepoint_zero_column, otu_table,
            ignore_val=ignore_val, filt=filt)
        output = output_results_paired_T_test(results, taxonomy_info)
    return output

def get_common_OTUs(parsed_otu_tables, filt, category_info, \
                        filter_all_samples, otu_include):
    """Searches all OTU tables in dir, returns common OTUs and their taxonomies
       Applies filter within each OTU table."""
    OTU_list = set()
    taxonomy_all_OTUs = {}
    count = 0

    # get list of all observed OTUs and their taxonomies
    for otu_table in parsed_otu_tables:
        count += 1
        sys.stdout.flush()
        taxonomy_info = get_taxonomy_info(otu_table)
        OTU_list_i = filter_OTUs(otu_table, filt,
                                 all_samples=filter_all_samples,
                                 category_mapping_info=category_info)
        if count==1:
            OTU_list = set(OTU_list_i)
        else:
            OTU_list &= set(OTU_list_i)
            
        for OTU in taxonomy_info.keys():
            taxonomy_all_OTUs[OTU] = taxonomy_info[OTU]

    #filter OTU_list with the otu_include list
    if not otu_include is None:
        otu_include = [line.strip() for line in otu_include]
        OTU_list &= set(otu_include)
        if len(OTU_list) == 0:
            raise ValueError("No OTUs remain after applying the filter. Try lowering the filter value (-f option)")

    # remove taxonomies for OTUs not in OTU_list
    for k in taxonomy_all_OTUs.keys():
        if not k in OTU_list:
            del(taxonomy_all_OTUs[k])
    return OTU_list, taxonomy_all_OTUs

def test_wrapper_multiple(test, parsed_otu_tables, category_mapping, category, \
                threshold, filt, otu_include=None, \
                otu_table_relative_abundance=False):
    """runs statistical test to look for category/OTU associations on multiple files.
       Unlike the test_wrapper() method, this method includes all OTUs, even when 
       some have zero counts.
    """
    mapping_data, header, comments = category_mapping
    category_info, category_values = \
        get_category_info(mapping_data, header, category, threshold)
    
    # if this is the g_test, disallow otus that are present in all samples 
    filter_all_samples = test == "g_test"

    OTU_list, taxonomy_all_OTUs = get_common_OTUs(parsed_otu_tables, filt, \
                                  category_info=category_info, \
                                  filter_all_samples=filter_all_samples, \
                                  otu_include=otu_include)

    all_results = {}
    count = 0
    for otu_table in parsed_otu_tables:
        count += 1
        sys.stdout.flush()

        if test == 'ANOVA' or test == 'correlation': 
            if not otu_table_relative_abundance:
                otu_table = otu_table.normObservationBySample()
        elif not test=='g_test':
            raise ValueError("An invalid test statistic was given. (-s option). Valid values are ANOVA, correlation, and g_test.")

        taxonomy_info = get_taxonomy_info(otu_table)

        if test == 'ANOVA':
            results = run_ANOVA_OTUs(OTU_list, category_info, otu_table,
                                     category_values)
        elif test == 'correlation':
            results = run_correlation_OTUs(OTU_list, category_info, otu_table)
        elif test == 'g_test':
            results = run_G_test_OTUs(OTU_list, category_info, otu_table,
                        category_values, suppress_warnings=True)
        for OTU in results.keys():
            if not all_results.has_key(OTU):
                all_results[OTU] = []
            all_results[OTU].append(results[OTU])
    
    # aggregate multiple results and create output string
    if test == 'ANOVA':
        all_results = aggregate_multiple_results_ANOVA(all_results)
        output = output_results_ANOVA(all_results, category_values, taxonomy_all_OTUs)
    elif test == 'correlation':
        all_results = aggregate_multiple_results_correlation(all_results)
        output = output_results_correlation(all_results, taxonomy_all_OTUs)
    elif test == 'g_test':
        all_results = aggregate_multiple_results_G_test(all_results)
        output = output_results_G_test(all_results, taxonomy_all_OTUs)
    return output










