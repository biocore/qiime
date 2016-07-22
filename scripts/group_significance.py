#!/usr/bin/env python
# File created on 15 Aug 2013
from __future__ import division
from functools import reduce

__author__ = "Will Van Treuren, Luke Ursell"
__copyright__ = "Copyright 2013, The QIIME project"
__credits__ = ["Will Van Treuren", "Luke Ursell", "Catherine Lozupone",
               "Jesse Stombaugh", "Doug Wendel", "Dan Knights", "Greg Caporaso",
               "Jai Ram Rideout", "Daniel McDonald"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Will Van Treuren"
__email__ = "wdwvt1@gmail.com"

from qiime.util import (parse_command_line_parameters, make_option,
                        sync_biom_and_mf)
from qiime.stats import (benjamini_hochberg_step_down,
                                   bonferroni_correction)
from qiime.otu_significance import (get_sample_cats, get_sample_indices,
                                    get_cat_sample_groups, group_significance_row_generator,
                                    group_significance_output_formatter,
                                    sort_by_pval, run_group_significance_test,
                                    TWO_GROUP_TESTS, GROUP_TEST_CHOICES)
from qiime.parse import parse_mapping_file_to_dict
from biom import load_table
from numpy import array, where, seterr, allclose

# set invalid comparisons error level to ignore. when nans are compared they
# frequently trigger this error. its not informative from the user perspective
# and likely will trigger a bunch of forum posts.
seterr(invalid='ignore')

script_info = {}
script_info['brief_description'] = """Compare OTU frequencies across sample groups"""
script_info['script_description'] = """
This script is used to compare OTU frequencies in sample groups and to ascertain
whether or not there are statistically significant differences between the OTU
abundance in the different sample groups.

The script will compare each OTU based
on the passed sample groupings to see if it is differentially represented. The
sample groupings are determined by the -c option. The script will group together
samples which have the same value in the mapping file under the header passed
with the -c option. Any samples that do not contain a value under the given
header will not be included in the comparison.
At a basic level, the script is constructing a OTUxSample
(rowXcolumn) contingency table, and testing whether or not each OTU is
differentially represented in cerstain groups of columns (determined by the
metadata category passed).

There are several important considerations with this script.

Errors and behind-the-scenes processing:

- This script will ignore samples that are found in the mapping file but do not
  contain information for the passed category. This would cause the mapping
  file to fail validate_mapping_file.py, but will not cause a failure here.
- This script will silenty ignore situations where the set of samples in the
  mapping file is a superset of the samples in the biom file. If the reverse is
  true, the script will error unless --biom_samples_are_superset is passed.
- This script will round P-values greater than 1 (after correcting for multiple
  comparisons) to 1.0.
- If your results file contains nans for p values its because one or more of
  the assumptions the selected test makes about the data was not met by the
  given OTU. The inverse of this statement is not guaranteed; just because the
  test worked on the data doesn't mean all its assumptions are met, just that
  enough assumptions are met so it doesn't fail (see below).

Filtering your OTU table prior to this script is important:

- Filtering out OTUs which are found in a low percentage of samples is a good
  idea before using this script. The old otu_category_significance script
  removed OTUs that were not found in at least 25 percent of samples. This
  prevents 0 variance errors and spurious significance for really low abundance
  OTUs and focuses the hypothesis discovery process on the abundant OTUs which
  are likely playing a larger role.

Test assumptions:

- This script tests that some basic assumptions of the given statistical test
  are met by the passed data. These 'assumption tests' are *necessary* not
  *sufficient* to ensure that the given statistical test you are applying is
  appropriate for the data. For instance, the script will error if you use the
  Mann-Whitney-U test and one of your group sizes is smaller than 20. It is
  likely that assumptions about the distribution of the data, the distribution
  of the variance, etc. are not robustly met. IT IS YOUR REPSONSIBILTY TO CHECK
  THAT YOU ARE USING AN APPROPRIATE TEST. For more information on assumptions
  made by the tests, please view the following resources:

  - Biometry by Sokal and Rolhf
  - Nonparamteric Statistical Methods by Hollander and Wolfe
  - Documentation in R and Scipy packages
  - Handbook of Biological Statistics by McDonald (available at
    http://udel.edu/~mcdonald/statintro.html)

The assumptions we check for:

- Kruskal-Wallis: No assumptions are checked for Kruskal-Wallis.
- G-test: We check that all the values in the table are non-negative and we
  check that each sample grouping contains at least one non 0 value. If either
  condition is not met we return G-stat, pval = (nan, nan) but don't error.
- Mann-Whitney-U: The number of data points in the combined samples is >= 20. If
  there are fewer than 21 data points the script will error. Although R gives
  exact values for less than 50 data points,  its definition of 'exact' is
  unclear since a conditional permutation calculation with ~ 10^14 calculations
  would be required. The normal approximation is suggested for more than 16 data
  points by Mann and Whitney 1947, and more than 20 in the Scipy documentation.
  The bootstrapped version of the test does not require >20 data points. If all
  the data ranks are tied for one of the groups, the function will return U,nan.
- ANOVA: No assumptions are checked for ANOVA. However, if the within group
  variance is 0, we return nan,nan.
- T-test: No assumptions are checked for the T-test. If no variance groups are
  detected None or nan will be returned.

The assumptions we do not check for are:

- Kruskal-Wallis: The scipy documentation indicates that each group must have
  at least 5 samples to make the Chi Squared approximation for the distribution
  of the H statistic be appropriate. R has no such requirement, and we do not
  implement the requirement here. The KW test does assume that the distributions
  from which the samples come are the same (although they may be non-normal)
  except for their location parameter, and we do not check this.
- G-test: we check that the data are counts rather than relative abundance.
- Mann-Whitney-U: Equality of variance between groups. Sample 1 is IID, Sample 2
  is IID. Sample 1 and Sample 2 are mutually independent.
- ANOVA: ANOVA assumes equality of variance between groups (homoscedasticity),
  normality of the residuals, and independence of the individual observations in
  the samples. None of these conditions are checked by this script.
- T-test: the t-test assumes that the samples come from populations which are
  normally distributed. that the variances of the sampled groups are ~ equal and
  that the individual observations are independent. None of these conditions are
  checked by this script.

Null and alternate hypothesis:

- G-test: The null hypothesis for the g_test (aka goodness of fit,
  log-likelihood ratio test) is that the frequency of any given OTU is equal
  across all sample groups. The alternate hypothesis is that the frequency of
  the OTU is not the same across all sample groups.
- Kruskal-Wallis: The null hypothesis is that the location paramater of the
  groups of abundances for a given OTU is the same. The alternate hypothesis is
  that at least one of the location parameters is different.
- ANOVA: the null hypothesis is that the means of the observations in the groups
  are the same, the alternate is that at least one is not.
- Mann-Whitney-U: the null hypothesis is that the distributions of the groups
  are equal, such that there is a 50 percent chance that a value from group1 is
  greater than a value from group2. The alternate is that the distributions are
  not the same.
- T-test: the null hypothesis is that the means of the two groups are the same
  versus the alternate that they are unequal.

The available tests are:

- ANOVA: one way analysis of variance. This test compares the within-group
  variance to the between-group variance in order to assess whether or not the
  sample groups have even frequencies of a given OTU. It generalizes the t-test
  to more than two groups. This is a parametric test whose assumptions are
  likely violated by data found in most gene surveys.

- kruskal_wallis: nonparametric ANOVA. This test is functionally an expansion of
  ANOVA to cases where the sample means are unequal and the distribution is not
  normal. The assumption that the distribution from which each group (within a
  single OTU) came is the same remains. This is a nonparametric test.

- g_test: goodness of fit or log-likelihood ratio test. This test compares the
  ratio of the OTU frequencies in the sample groups to an 'extrinsic hypothesis'
  about what their distribution should be. The extrinsic hypothesis coded in this
  script is that all sample groups have equal OTU frequencies. The test compares
  the ratio of the observed OTU frequencies in the sample groups to the expected
  frequencies based on the extrinsic hypothesis. This is a parametric test.

- parametric_t_test: Student's t-test. This test compares the frequencies of an
  OTU in one sample group versus another sample group to see what the probability
  of drawing the samples given that each sample had an equal proportion of the OTU
  in it. This is a parametric test whose assumptions are likely violated by data
  found in most gene surveys.

- nonparametric_t_test: nonparametric t-test is calculated using Monte Carlo
  simulation. This test performs in the same way as the parametric t-test, but
  computes the probability based on a boot-strap procedure where the sample
  group values are permuted. The fraction of the time that a t-statistic
  greater than or equal to the observed t-statistic is found is the basis of
  the nonparametric p-value. This is a nonparametric test.

- mann_whitney_u: aka Wilcoxon rank sum test is a nonparametric test where the
  null hypothesis is that the populations from which the two samples come have
  equal means. It is basically an extension of the t-test. This is a nonparametric
  test.

- bootstrap_mann_whitney_u: the bootstrapped version of the mann_whitney_u test.
  Identical behavior to the nonparametric_t_test. This is a nonparametric test.

"""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("Find which OTUs have the highest probablilty of being differently represented depending on the sample category 'diet' using a G test:",
     "",
     "%prog -i otu_table.biom -m map_overlapping.txt -c diet -s g_test -o gtest_ocs.txt"))
script_info['script_usage'].append(
    ("Find which OTUs are differentially represented in two sample groups 'before_after' using a T-test:",
     "",
     "%prog -i otu_table.biom -m map_overlapping.txt -c before_after -s parametric_t_test -o tt_ocs.txt"))
script_info['script_usage'].append(
    ("Find which OTUs are differentially represented in the sample groups formed by 'diet' based on nonparamteric ANOVA, aka, Kruskal Wallis test. In addition, prevent the script from erroring because the biom table samples are a superset of the mapping file samples, and print the non-overlapping samples:",
     "",
     "%prog -i otu_table.biom -m map.txt -c diet -s kruskal_wallis -o kw_ocs.txt --biom_samples_are_superset --print_non_overlap"))
script_info['script_usage'].append(
    ("Find which OTUs are differentially represented in the sample groups formed by 'before_after' based on bootstrapped T-testing with 100 permutations:",
     "",
     "%prog -i otu_table.biom -m map_overlapping.txt -c before_after -s nonparametric_t_test --permutations 100 -o btt_ocs.txt"))

script_info['output_description'] = """
This script generates a tab separated output file with the following headers:

- OTU: OTU id
- Test-Statistic: the value of the test statistic for the given test
- P: the raw P value returned by the given test.
- FDR_P: the P value corrected by the Benjamini-Hochberg FDR procedure for
  multiple comparisons.
- Bonferroni_P: the P value corrected by the Bonferroni procedure for multiple
  comparisons.
- groupX_mean: there will be as many of these headers as there are unique values
  in the mapping file under the category passed with the -c option. Each of these
  fields will contain the mean frequency/abundance/count of the given OTU for the
  given sample group.
- Taxonomy: this column will be present only if the biom table contained Taxonomy
  information. It will contain the taxonomy of the given OTU.

"""
script_info['required_options'] = [
    make_option('-i', '--otu_table_fp',
                help='path to biom format table',
                type='existing_path'),
    make_option('-m', '--mapping_fp', type='existing_filepath',
                help='path to category mapping file'),
    make_option('-c', '--category', type='string',
                help='name of the category over which to run the analysis'),
    make_option('-o', '--output_fp', type='new_filepath',
                help='path to the output file')]

script_info['optional_options'] = [
    make_option(
        '-s', '--test', type="choice", choices=GROUP_TEST_CHOICES.keys(),
        default='kruskal_wallis', help='Test to use. Choices are: %s' %
        (', '.join(GROUP_TEST_CHOICES.keys())) + ' [default: %default]'),
    make_option('--metadata_key', default='taxonomy', type=str,
                help='Key to extract metadata from biom table. default: %default]'),
    make_option('--permutations', default=1000, type=int,
                help='Number of permutations to use for bootstrapped tests.' +
                '[default: %default]'),
    make_option('--biom_samples_are_superset', action='store_true',
                default=False,
                help='If this flag is passed you will be able to use a biom table ' +
                'that contains all the samples listed in the mapping file ' +
                'as well as additional samples not listed in the mapping file. ' +
                'Only their intersecting samples will be used for calculations.'),
    make_option('--print_non_overlap', action='store_true', default=False,
                help='If this flag is passed the script will display the samples that' +
                ' do not overlap between the mapping file and the biom file.')]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    # sync the mapping file and the biom file
    tmp_bt = load_table(opts.otu_table_fp)
    tmp_pmf, _ = parse_mapping_file_to_dict(opts.mapping_fp)
    pmf, bt, nonshared_samples = sync_biom_and_mf(tmp_pmf, tmp_bt)

    # test error conditions for overlapping mf and bt
    if not opts.biom_samples_are_superset:
        # user indicates biom sample should be subset of mapping file samples
        if any([i in nonshared_samples for i in tmp_bt.ids()]):
            raise ValueError('The samples in the biom table are a superset of' +
                             ' the samples in the mapping file. The script will abort in' +
                             ' this case even though the calculations wouldn\'t be' +
                             ' affected, to ensure consistency within QIIME. Pass the' +
                             ' --biom_samples_are_superset option to disable this behavior.')
    # user wants non-overlapping samples printed out
    if opts.print_non_overlap:
        print 'The following samples were not shared between the mapping file' +\
            ' and the biom file and will not be included in the analysis:\n' +\
            ' '.join(nonshared_samples)

    # find group indices
    sam_cats = get_sample_cats(pmf, opts.category)
    cat_sam_groups = get_cat_sample_groups(sam_cats)
    cat_sam_indices = get_sample_indices(cat_sam_groups, bt)

    # sanity check to prevent inscrutable errors later
    if not all([len(v) > 0 for k, v in cat_sam_indices.items()]):
        raise ValueError('At least one metadata group has no samples. Check ' +
                         'that the mapping file has at least one sample for each value in ' +
                         'the passed category.')
    if opts.test in TWO_GROUP_TESTS and len(cat_sam_indices) > 2:
        option_parser.error('The t-test and mann_whitney_u test may ' +
                            'only be used when there are two sample groups. Choose another ' +
                            'test or another metadata category.')

    # check that assumptions are met for a given test:
    if opts.test == 'mann_whitney_u':
        sams = reduce(lambda x, y: len(x) + len(y), cat_sam_indices.values())
        if sams <= 20:
            raise ValueError('The number of samples is too small to use the ' +
                             'Mann-Whitney-U normal approximation. Review the script ' +
                             'documentation.')

    # check that the G-test was not selected if the table appears to be
    # relative abundance
    if opts.test == 'g_test':
        if allclose(bt.sum(axis='sample'), 1.) or (bt.sum(axis='whole') == 1.):
            raise ValueError('It appears that the biom table you have passed '
                'is a relative abundance table where values i,j (obsevation i '
                'count in sample j) are fractional and the sum of the columns '
                'is 1.0. This will fail to work properly with the G-test. If '
                'your data sums to 1 in each column but your data is not '
                'relative abundance then the tests will fail anyway because '
                'of the reduced number of observations.')

    # run actual tests
    data_feed = group_significance_row_generator(bt, cat_sam_indices)
    test_stats, pvals, means = run_group_significance_test(
        data_feed, opts.test,
        GROUP_TEST_CHOICES, int(opts.permutations))

    # calculate corrected pvals
    fdr_pvals = array(benjamini_hochberg_step_down(pvals))
    bon_pvals = bonferroni_correction(pvals)
    # correct for cases where values above 1.0 due to correction
    fdr_pvals = where(fdr_pvals > 1.0, 1.0, fdr_pvals)
    bon_pvals = where(bon_pvals > 1.0, 1.0, bon_pvals)

    # write output results after sorting
    lines = group_significance_output_formatter(bt, test_stats, pvals,
                                                fdr_pvals, bon_pvals, means, cat_sam_indices, md_key=opts.metadata_key)
    lines = sort_by_pval(lines, ind=2)
    o = open(opts.output_fp, 'w')
    o.writelines('\n'.join(lines))
    o.close()

if __name__ == "__main__":
    main()
