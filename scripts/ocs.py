#!/usr/bin/env python
# File created on 15 Aug 2013
from __future__ import division

__author__ = "Will Van Treuren, Luke Ursell"
__copyright__ = "Copyright 2013s, The QIIME project"
__credits__ = ["Will Van Treuren, Luke Ursell"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Will Van Treuren, Luke Ursell"
__email__ = "wdwvt1@gmail.com"
__status__ = "Development"

from qiime.util import parse_command_line_parameters, make_option
from qiime.ocs import (sync_biom_and_mf, get_sample_cats, get_sample_indices, 
    get_cat_sample_groups, row_generator, output_formatter, fdr_correction, 
    bonferroni_correction, sort_by_pval, run_ocs_test, two_group_tests)
from qiime.pycogent_backports.test import (G_fit, ANOVA_one_way, kruskal_wallis,
    mw_test, mw_boot)
from cogent.maths.stats.test import (t_two_sample, mc_t_two_sample)
from qiime.parse import parse_mapping_file_to_dict
from biom.parse import parse_biom_table
from numpy import array, where

test_choices = {'ANOVA': ANOVA_one_way, 'g_test': G_fit, 
    'kruskal_wallis': kruskal_wallis, 'parametric_t_test': t_two_sample,
    'nonparametric_t_test': mc_t_two_sample, 'mann_whitney_u': mw_test, 
    'bootstrap_mann_whitney_u': mw_boot}

script_info = {}
script_info['brief_description'] = """
This script is used to compare OTU frequencies in sample groups and to ascertain 
whether or not there are statistically significant differences between the OTU
abundance in the different sample groups."""
script_info['script_description'] = """
This script is used to compare OTU frequencies in sample groups and to ascertain 
whether or not there are statistically significant differences between the OTU
abundance in the different sample groups. The script will compare each OTU based
on the passed sample groupings to see if it is differentially represented. The 
sample groupings are determined by the -c option. The script will group together
samples which have the same value in the mapping file under the header passed 
with the -c option. At a basic level, the script is constructing a OTUxSample
(rowXcolumn) contingency table, and testing whether or not each OTU is 
differentially represented in certain groups of columns (determined by the 
metadata category passed). 

There are several important considerations with this script:
* The script will silenty ignore samples in the mapping file for which there are
* no values for the given category. This would cause the mapping file to fail 
* validate_mapping.py, but we check here regardless. 

* The script will only consider samples that are found in both the mapping file
* and the biom table.

* The null hypothesis for the g_test (aka goodness of fit, log-likelihood ratio
* test) is that the frequency of any given OTU is equal across all sample 
* groups.

* corrected pvalues greater than 1 will be rounded to 1

The available tests are:
ANOVA - one way analysis of variance. This test compares the within-group 
variance to the between-group variance in order to assess whether or not the 
sample groups have even frequencies of a given OTU. It generalizes the t-test to
more than two groups. This is a parametric test whose assumptions are likely 
violated by data found in most gene surveys. 

kruskal_wallis - nonparametric ANOVA. This test is functionally an expansion of
ANOVA to cases where the sample means are not equal (although other assumptions
like equal skewness remain). This is a nonparametric test. 

g_test - goodness of git or log-likelihood ratio test. This test compares the 
ratio of the OTU frequencies in the sample groups to an 'extrinsic hypothesis' 
about what their distribution should be. The extrinsic hypothesis coded in this
script is that all sample groups have equal OTU frequencies. The test compares
the ratio of the observed OTU frequencies in the sample groups to the expected
frequencies based on the extrinsic hypothesis. This is a XXXX test. 

parametric_t_test - Student's t-test. This test compares the frequencies of an
OTU in one sample group versus another sample group to see what the probability 
of drawing the samples given that each sample had an equal proportion of the OTU
in it. This is a parametric test whose assumptions are likely violated by data 
found in most gene surveys.

nonparametric_t_test - nonparametric t-test calculated using Monte Carlo 
simulation. This test performs in the same way as the t-test, but computes the 
probability based on a boot-strap procedure where the sample group values are 
permuted. The fraction of the time that a t-statistic greater than or equal to
the observed t-statistic is found is the basis of the nonparametric p-value. 
This is a nonparametric test.

mann_whitney_u - aka wilcoxon rank sum test is a nonparametric test that the
means of two populations represented by samples are even (null hypothesis). It 
is basically an extension of the t-test. This is a nonparametric test.

bootstrap_mann_whitney_u - the bootstrapped version of the mann_whitney_u test. 
Identical behavior to the nonparametric_t_test. This is a nonparametric_t_test.

This script generates a tab separated output file with the following headers.
OTU - OTU id 
Test-Statistic - the value of the test statistic for the given test
P - the raw P value returned by the given test. 
FDR_P - the P value corrected by the Benjamini-Hochberg FDR procedure for 
 multiple comparisons.
Bonferroni_P - the P value corrected by the Bonferroni procedure for multiple
 comparisons.
groupX_mean - there will be as many of these headers as there are unique values
 in the mapping file under the category passed with the -c option. Each of these
 fields will contain the mean frequency/abundance/count of the given OTU for the
 given sample group.
Taxonomy - this column will be present only if the biom table contained Taxonomy
 information. It will contain the taxonomy of the given OTU. 
"""
script_info['script_usage'] = []
script_info['script_usage'].append(("add", "me", "please"))
script_info['output_description']= """
This script generates a tab separated output file with the following headers.
OTU - OTU id 
Test-Statistic - the value of the test statistic for the given test
P - the raw P value returned by the given test. 
FDR_P - the P value corrected by the Benjamini-Hochberg FDR procedure for 
 multiple comparisons.
Bonferroni_P - the P value corrected by the Bonferroni procedure for multiple
 comparisons.
groupX_mean - there will be as many of these headers as there are unique values
 in the mapping file under the category passed with the -c option. Each of these
 fields will contain the mean frequency/abundance/count of the given OTU for the
 given sample group.
Taxonomy - this column will be present only if the biom table contained Taxonomy
 information. It will contain the taxonomy of the given OTU. 
"""
script_info['required_options']=[
    make_option('-i','--otu_table_fp',
        help='path to biom format table or to directory containing OTU tables',
        type='existing_path'),
    make_option('-m','--mapping_fp', type='existing_filepath',
        help='path to category mapping file'),
    make_option('-c', '--category', type='string',
        help='name of the category over which to run the analysis'),
    make_option('-o', '--output_fp', type='new_filepath',
        help='path to the output file or directory')]

script_info['optional_options']=[
    make_option('-s', '--test', type="choice", choices=test_choices.keys(),
        default='ANOVA', help='Test to use. Choices are:\n%s' % \
         (', '.join(test_choices.keys()))+'\n\t' + '[default: %default]'),
    make_option('-w', '--collate_results', dest='collate_results',
        action='store_true', default=False,
        help='When passing in a directory of OTU tables, '
        'this parameter gives you the option of collating those resulting '
        'values. For example, if your input directory contained multiple '
        'rarefied OTU tables at the same depth, pass the -w option '
        'in order to find the average p-value for your statistical test '
        'over all rarefied tables and collate the results into one file. '
        'If your input directory contained OTU tables '
        'that contained different taxonomic levels, filtering levels, etc '
        'then do not pass the -w option so that an individual results file '
        'is created for every input OTU table. [default: %default]')]

script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # sync the mapping file and the biom file
    bt = parse_biom_table(open(opts.otu_table_fp))
    pmf, _ = parse_mapping_file_to_dict(opts.mapping_fp)
    pmf, bt = sync_biom_and_mf(pmf, bt)
    # find group indices
    sam_cats = get_sample_cats(pmf, opts.category)
    cat_sam_groups = get_cat_sample_groups(sam_cats)
    cat_sam_indices = get_sample_indices(cat_sam_groups, bt)
    # sanity check to prevent inscrutable errors later
    if not all([len(v)>0 for k,v in cat_sam_indices.items()]):
        raise ValueError('At least one metadata group has no samples. Check '+\
            'that the mapping file has at least one sample for each value in '+\
            'the passed category.')
    
    if opts.test in two_group_tests and len(cat_sam_indices) > 2:
            raise option_parser.error('The t-test and mann_whitney_u test may'+\
                ' only be used when there'+\
                ' are two sample groups. Choose another test or another ' +\
                'metadata category.')

    data_feed = row_generator(bt, cat_sam_indices)
    test_stats, pvals, means = run_ocs_test(data_feed, opts.test, test_choices)
    # calculate corrected pvals
    fdr_pvals = array(fdr_correction(pvals))
    bon_pvals = bonferroni_correction(pvals)
    # correct for cases where values above 1.0 due to correction
    fdr_pvals = where(fdr_pvals>1.0, 1.0, fdr_pvals)
    bon_pvals = where(bon_pvals>1.0, 1.0, bon_pvals)
    # write output results after sorting
    lines = output_formatter(bt, test_stats, pvals, fdr_pvals, bon_pvals, means,
        cat_sam_indices)
    lines = sort_by_pval(lines, ind=2)
    o = open(opts.output_fp, 'w')
    o.writelines('\n'.join(lines))
    o.close()

if __name__ == "__main__":
    main()

