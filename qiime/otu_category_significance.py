#!/usr/bin/env python
# File created on 07 Oct 2009.
from __future__ import division

__author__ = "Catherine Lozupone"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Catherine Lozupone and Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Catherine Lozupone"
__email__ = "lozupone@colorado.edu"
__status__ = "Pre-release"


from optparse import OptionParser
from os.path import split, splitext
from numpy import argsort
from cogent.util.dict2d import Dict2D
from cogent.maths.stats.test import calc_contingency_expected, G_fit_from_Dict2D,\
    ANOVA_one_way
from cogent.maths.stats.util import Numbers
from numpy import array
import sys

usage_str = """usage: %prog [-o output_file] {-i OTU table, -m category mapping -f filter -c category, -t None}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example usage:
    Look for OTUs that are associated with a category. Currently can do:
    1) perform g-test of independence to determine whether OTU presence
    absence is associated with a category in the mapping file. Can
    be used to look for co-occurrence using qPCR data or to look for 
    associations with a categorical environmental variable.
    2) perform ANOVA to determine whether OTU abundance is significantly
    different across a category
python ~/repo/Qiime/qiime/OTU_category_significance.py -i otu_table.txt, -m category_mapping.txt -s g_test -f 10 -c category name -o output_fp -t None
"""

def parse_command_line_parameters():
    """ Parses command line arguments """
    usage = usage_str
    version = 'Version: %prog ' + __version__
    parser = OptionParser(usage=usage, version=version)

    # A binary 'verbose' flag
    parser.add_option('-v','--verbose',action='store_true',\
        dest='verbose',help='Print information during execution -- '+\
        'useful for debugging [default: %default]')

    parser.add_option('-i','--otu_table_fp', dest='otu_table_fp',\
        help='path to the otu table [REQUIRED]')

    parser.add_option('-m','--category_mapping_fp',\
        dest='category_mapping_fp',\
        help='path to category mapping file [REQUIRED]')

    parser.add_option('-c','--category', dest='category',\
        help='name of category over which to run the analysis [REQUIRED]')
    
    parser.add_option('-s','--test', dest='test', default='g_test',\
        help='the type of statistical test to run. options are: ' +\
            'g_test: g test of independence: determines whether OTU ' +\
                'presence/absence is associated with a category ' +\
            'ANOVA: determines whether OTU abundance is associated with a ' +\
                'category')
    
    parser.add_option('-o','--output_fp', dest='output_fp', \
        default= 'otu_category_G_test_results.txt',\
        help='path to output file. otu_category_G_test_results.txt by default')

    parser.add_option('-f','--filter', dest='filter',\
        default= 10, \
        help='minimum number of samples that must contain the OTU for the ' +\
        'OTU to be included in the analysis. default value=10.')
    
    parser.add_option('-t','--threshold', dest='threshold', default=None, \
        help='threshold under which to consider something absent: ' +\
        'Only used if you have numerical data that should be converted to ' +\
        'present or absent based on a threshold. Should be None for ' +\
        'categorical data. default value is None')
    # Set default values here if they should be other than None
    parser.set_defaults(verbose=False)

    required_options = ['otu_table_fp', 'category_mapping_fp', 'category']

    opts,args = parser.parse_args()
    
    for option in required_options:
        if eval('opts.%s' % option) == None:
            parser.error('Required option --%s omitted.' % option) 

    return opts,args

def parse_otu_table(otu_table):
    """parse otu_file to get dict of OTU names mapped to sample:count dicts
    
    if taxonomic information is in there parses that too
    """
    result = {}
    taxonomy_info = {}
    for line in otu_table:
        if line.startswith('#OTU ID'):
            sample_names = line.strip().split('\t')
            if not sample_names[-1] == 'Consensus Lineage':
                num_samples = len(sample_names) - 1
                taxonomy = True
            else:
                num_samples = len(sample_names) - 2
                taxonomy = False
        elif line and not line.startswith('#'):
            line = line.strip().split('\t')
            OTU_name = line[0]
            result[OTU_name] = {}
            if taxonomy:
                for index, sample in enumerate(sample_names):
                    if index != 0:
                        result[OTU_name][sample] = line[index]
            else:
                for index, sample in enumerate(sample_names):
                    if index != 0 and index != len(sample_names) -1:
                        result[OTU_name][sample] = line[index]
                taxonomy_info[OTU_name] = line[-1]
    return result, num_samples, taxonomy_info

def convert_OTU_table_relative_abundance(otu_table):
    """converts the OTU table to have relative abundances rather than raw counts
    """
    output = []
    data_lines = []
    otu_ids = []
    tax_strings = []
    taxonomy=False
    for line in otu_table:
        line = line.strip().split('\t')
        if line[0].startswith('#OTU ID'):
            output.append('\t'.join(line))
            if line[-1] == 'Consensus Lineage':
                taxonomy=True
        elif line[0].startswith('#'):
            output.append('\t'.join(line))
        else:
            if taxonomy:
                vals = [float(i) for i in line[1:-1]]
                tax_strings.append(line[-1])
            else:
                vals = [float(i) for i in line[1:]]
                tax_string = None
            data = array(vals, dtype=float)
            data_lines.append(data)
            otu_ids.append(line[0])
    data_lines = array(data_lines)
    totals = sum(data_lines)
    new_values = []
    for i in data_lines:
        new_values.append(i/totals)
    for index, i in enumerate(new_values):
        line = [otu_ids[index]]
        line.extend([str(j) for j in i])
        if taxonomy:
            line.append(tax_strings[index])
        output.append('\t'.join(line))
    return output

def parse_category_mapping(category_mapping, category, threshold=None):
    """parse category mapping file
    
    returns a dict mapping the sample name to the value of category

    the values in the category mapping must be categorical. If numerical
    threshold is specified, converts to 1 or 0 depending on whether the 
    value is below that threshold"""
    result = {}
    if threshold and threshold != 'None':
        category_values = ['0', '1']
    else:
        category_values = []
    
    for line in category_mapping:
        if line.startswith("#SampleID"):
            line = line.strip().split('\t')
            category_index = line.index(category)
        elif not line.startswith('#'):
            line = line.strip().split('\t')
            sample_name = line[0]
            if threshold and threshold != 'None':
                val = float(line[category_index])
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

def filter_OTUs(OTU_sample_info, filter, num_samples=None, all_samples=True,\
                category_mapping_info=None):
    """Get the list of OTUs found in more than <filter> samples

    optionally filters out those that are found in all_samples if True:
        G_test: set to True since can't show presence/absence patterns if
            present in all (also causes an error)
        ANOVA: set to False since can still see differences in abundance

    optionally takes a category mapping file as input and only considers
        inclusion in samples that are included therein. This is a method
        of making sure that the samples represented in the category
        mapping file and the OTU table are in sync
    """
    result = []
    if category_mapping_info:
        included_samples = category_mapping_info.keys()
        num_samples = len(included_samples)
    for OTU in OTU_sample_info:
        samples = []
        for sample in OTU_sample_info[OTU]:
            if float(OTU_sample_info[OTU][sample]) > 0:
                if category_mapping_info:
                    if sample in included_samples:
                        samples.append(sample)
                else:
                    samples.append(sample)
        if len(samples) > int(filter):
            if all_samples:
                if len(samples) < num_samples:
                    result.append(OTU)
            else:
                result.append(OTU)
    return result

def run_single_G_test(OTU_name, category_info, otu_sample_info, category_values):
    """run the G test on a single OTU
    """
    contingency_matrix = make_contingency_matrix(OTU_name, category_info, \
        otu_sample_info, category_values)
    contingency_matrix = calc_contingency_expected(contingency_matrix)
    g_val, prob = G_fit_from_Dict2D(contingency_matrix)
    return g_val, prob, contingency_matrix

def make_contingency_matrix(OTU_name, category_info, otu_sample_info, category_values):
    """make the contingency table for running the G test of independence
    
    counts OTU as present (count > 1) or absent
    makes a column in the matrix for each category value
    """
    result = {'OTU_pos':{},
                'OTU_neg':{}}
    for category in category_values:
        result['OTU_pos'][category + '_pos'] = 0
        result['OTU_neg'][category + '_pos'] = 0
    for sample in category_info:
        category = category_info[sample]
        try:
            OTU_count = int(otu_sample_info[OTU_name][sample])
            worked = True
        except KeyError: 
            print "Warning: {0} is in the  sample mapping file but not the OTU table" 
            worked = False
        if worked:
            if OTU_count == 0:
                result['OTU_neg'][category + '_pos'] += 1
            elif OTU_count > 0: 
                result['OTU_pos'][category + '_pos'] += 1
    return Dict2D(result, Default=0, Pad=True)

def run_G_test_OTUs(OTU_list, category_info, otu_sample_info, category_values):
    """run the G test on all OTUs in the OTU list
    """
    result = {}
    num_comparisons = len(OTU_list)
    for OTU in OTU_list:
        g_val, g_prob, contingency_matrix = run_single_G_test(OTU, \
            category_info, otu_sample_info, category_values)
        Bonf_p_val = g_prob * num_comparisons
        result[OTU] = [g_val, g_prob, contingency_matrix, Bonf_p_val]
    return result

def run_ANOVA_OTUs(OTU_list, category_info, otu_sample_info, category_values):
    """run ANOVA on all OTUs in the OTU list
    """
    result = {}
    num_comparisons = len(OTU_list)
    for OTU in OTU_list:
        group_means, prob = run_single_ANOVA(OTU, category_info, \
            otu_sample_info, category_values)
        Bonf_p_val = prob * num_comparisons
        result[OTU] = [group_means, prob, Bonf_p_val]
    return result

def run_single_ANOVA(OTU, category_info, otu_sample_info, category_values):
    """runs ANOVA on the designated OTU
    """
    result = {}
    #get a list of values for each category
    values = []
    for category in category_values:
        values.append(Numbers([]))
    sample_info = otu_sample_info[OTU]
    for sample in sample_info:
        count = sample_info[sample]
        if sample in category_info:
            category = category_info[sample]
            index = category_values.index(category)
            values[index].append(count)
    dfn, dfd, F, between_MS, within_MS, group_means, prob = ANOVA_one_way(values)
    return group_means, prob

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
    G_test_results = fdr_correction_G_test(G_test_results)
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
    return output

def output_results_ANOVA(ANOVA_results, category_values, taxonomy_info=None):
    """creates the results output using result of run_ANOVA_OTUs"""
    header = ['OTU', 'prob', 'Bonferroni_corrected', 'FDR_corrected']
    for i in category_values:
        header.append(i + '_mean')
    if taxonomy_info:
            header.append('Consensus Lineage')
    output = ['\t'.join(header)]
    #perform fdr correction
    ANOVA_results = fdr_correction_ANOVA(ANOVA_results)
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
    return output

def fdr_correction_G_test(G_test_results):
    """corrects G test results using the false discovery rate method.

    ranks the p-values from low to high. multiplies each p-value by the #
    of comparisons divided by the rank.
    """
    names = []
    g_probs = []
    for i in G_test_results:
        names.append(i)
        g_probs.append(G_test_results[i][1])
    corrected_probs = fdr_correction(g_probs)
    for index, prob in enumerate(corrected_probs):
        G_test_results[names[index]].append(prob)
    return G_test_results

def fdr_correction_ANOVA(ANOVA_results):
    """corrects ANOVA probs using the false discovery rate method.

    ranks the p-values from low to high. multiplies each p-value by the #
    of comparisons divided by the rank.
    """
    names = []
    probs = []
    for i in ANOVA_results:
        names.append(i)
        probs.append(ANOVA_results[i][1])
    corrected_probs = fdr_correction(probs)
    for index, prob in enumerate(corrected_probs):
        ANOVA_results[names[index]].append(prob)
    return ANOVA_results

def fdr_correction(probs):
    """corrects a list of probs using the false discovery rate method

    ranks the p-values from low to high. multiplies each p-value by the #
    of comparison divided by the rank.
    """
    corrected_probs = [None] * len(probs)
    for rank, index in enumerate(argsort(probs)):
        correction = len(probs) / float(rank + 1)
        fdr_p = probs[index] * correction
        corrected_probs[index] = fdr_p
    return corrected_probs

def G_test_wrapper(otu_table, category_mapping, category, threshold, \
                filter, output_fp):
    """runs the G test to look for category/OTU associations"""
    otu_sample_info, num_samples, taxonomy_info = parse_otu_table(otu_table)
    category_info, category_values = parse_category_mapping(category_mapping, category, threshold)
    OTU_list = filter_OTUs(otu_sample_info, filter, all_samples= True,\
        category_mapping_info=category_info)
    G_test_results = run_G_test_OTUs(OTU_list, category_info, otu_sample_info, category_values)
    output = output_results_G_test(G_test_results, taxonomy_info)
    of = open(output_fp, 'w')
    of.write('\n'.join(output))

def ANOVA_wrapper(otu_table, category_mapping, category, threshold, \
                filter, output_fp):
    """runs ANOVA to look for category/OTU associations"""

    otu_table = convert_OTU_table_relative_abundance(otu_table)
    otu_sample_info, num_samples, taxonomy_info = parse_otu_table(otu_table)
    category_info, category_values = parse_category_mapping(category_mapping, category, threshold)
    OTU_list = filter_OTUs(otu_sample_info, filter, all_samples= False, \
        category_mapping_info=category_info)
    ANOVA_results = run_ANOVA_OTUs(OTU_list, category_info, otu_sample_info, category_values)
    output = output_results_ANOVA(ANOVA_results, category_values, taxonomy_info)
    of = open(output_fp, 'w')
    of.write('\n'.join(output))

if __name__ == "__main__":
    opts,args = parse_command_line_parameters()
    verbose = opts.verbose
    
    otu_table_fp = opts.otu_table_fp
    otu_table = open(otu_table_fp)
    output_fp = opts.output_fp
    
    category_mapping_fp = opts.category_mapping_fp
    category_mapping = open(category_mapping_fp,'U')
    
    filter = opts.filter
    category = opts.category
    threshold = opts.threshold
    if threshold and threshold != 'None':
        threshold = float(threshold)
    test = opts.test

    if test == 'g_test':
        G_test_wrapper(otu_table, category_mapping, category, threshold, \
                filter, output_fp)
    elif test == 'ANOVA':
        ANOVA_wrapper(otu_table, category_mapping, category, threshold, \
                filter, output_fp)
