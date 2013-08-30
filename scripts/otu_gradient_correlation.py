#!/usr/bin/env python
# File created on 15 Aug 2013
from __future__ import division

__author__ = "Will Van Treuren, Luke Ursell"
__copyright__ = "Copyright 2013, The QIIME project"
__credits__ = ["Will Van Treuren", "Luke Ursell", "Catherine Lozupone",
    "Jesse Stombaugh", "Doug Wendel", "Dan Knights", "Greg Caporaso", 
    "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Will Van Treuren"
__email__ = "wdwvt1@gmail.com"
__status__ = "Development"

from qiime.util import (parse_command_line_parameters, make_option,
    sync_biom_and_mf)
from qiime.pycogent_backports.test import (benjamini_hochberg_step_down, 
    bonferroni_correction)
from qiime.otu_significance import (sort_by_pval, run_correlation_test, 
    correlation_row_generator, correlation_output_formatter, 
    CORRELATION_TEST_CHOICES, paired_t_generator, run_paired_t, 
    paired_t_output_formatter, longitudinal_row_generator,
    run_longitudinal_correlation_test, longitudinal_correlation_formatter, 
    CORRELATION_TEST_CHOICES, get_sample_cats, get_cat_sample_groups, 
    get_sample_indices)
from qiime.parse import parse_mapping_file_to_dict
from biom.parse import parse_biom_table
from numpy import array, where


script_info = {}
script_info['brief_description'] = """
This script calculates the correlation between OTU values and a gradient of 
sample data. Several methods are provided to allow the user to correlate OTUs to
sample metadata values. Longitudinal correlations are also supported, where a 
single sample represents a reference point. Finally, the script allows one to 
conduct a paired t test.
"""
script_info['script_description'] = """
This script calculates the correlation between OTU values and a gradient of 
sample data. Several methods are provided to allow the user to correlate OTUs to
sample metadata values. Longitudinal correlations are also supported, where a 
single sample represents a reference point. Finally, the script allows one to 
conduct a paired t test.
The tests of OTU correlation to a metadata field are accomplished by passing a 
mapping file and a category from which to pull the data. If the data are not 
convertable to floats, the script will abort. 
The tests of OTU correlation to a metadata field with a longitudinal component 
are accomplished by passing the --individual_column which tells the script which
samples are from which individual, along with the -c option which tells the 
script which metadata field/column to use as the gradient for the longitudinal 
correlation. This will most frequently be useful when there is time series data
and you have sampes from an individual at a host of time points.  
The paired t test is accomplished by passing a paired mapping file which is just
a two column (tab separation) table with the samples that should be paired in 
each row. It should not have headers.
The available tests are Kendall's Tau, Spearmans rank correlation, and Pearson. 
This script generates a tab separated output file which differs based on which 
test you have chosen. If you have chosen simple correlation or paired_t (i.e. 
you have not passed the --individual_column option ) then you will see
the following headers:
OTU - OTU id 
Test-Statistic - the value of the test statistic for the given test
P - the raw P value returned by the given test. 
FDR_P - the P value corrected by the Benjamini-Hochberg FDR procedure for 
 multiple comparisons.
Bonferroni_P - the P value corrected by the Bonferroni procedure for multiple
 comparisons.
Taxonomy - this column will be present only if the biom table contained Taxonomy
 information. It will contain the taxonomy of the given OTU.
If you have opted for the longitudinal correlation test (with the
--individual_column) you will have:
OTU
Individual:X stat - correlation statistics from the given test. There will be as
 many of these headers as there are individuals in the individual column.
Taxonomy
"""
script_info['script_usage'] = []
script_info['script_usage'].append(("Calculate the correlation between OTUs in the table and the pH of the samples from mich they came:", "", "%prog -i otu_table.biom -m map.txt -c pH -s spearman -o spearman_otu_gradient.txt"))
script_info['script_usage'].append(("Calculate correlation between OTUs over the course of a treatment regimen assuming that 'hsid' is the column in the mapping file which identifies the individual subject from which each sample came and 'treatment_day' refers to the day after treatment was started for each individual:", "", "%prog -i otu_table.biom -m map.txt -c treatment_day -s kendall --individual_column hsid -o kendall_longitudinal_otu_gradient.txt"))
script_info['script_usage'].append(("Calculate paired t values for a before and after group of samples:", "", "%prog -i otu_table.biom --paired_t_fp=paired_samples.txt -o kendall_longitudinal_otu_gradient.txt"))

script_info['output_description']= """
This script generates a tab separated output file which differs based on which 
test you have chosen. If you have chosen simple correlation or paired_t (i.e. 
you have not passed the --individual_column option) then you will see
the following headers:
OTU - OTU id 
Test-Statistic - the value of the test statistic for the given test
P - the raw P value returned by the given test. 
FDR_P - the P value corrected by the Benjamini-Hochberg FDR procedure for 
 multiple comparisons.
Bonferroni_P - the P value corrected by the Bonferroni procedure for multiple
 comparisons.
Taxonomy - this column will be present only if the biom table contained Taxonomy
 information. It will contain the taxonomy of the given OTU.
If you have opted for the longitudinal correlation test (with the
--individual_column) you will have:
OTU
Individual:X stat - correlation statistics from the given test. There will be as
 many of these headers as there are individuals in the individual column.
Taxonomy
"""
script_info['required_options']=[
    make_option('-i','--otu_table_fp',
        help='path to biom format table or to directory containing OTU tables',
        type='existing_path'),
    make_option('-o', '--output_fp', type='new_filepath',
        help='path to the output file or directory')]

script_info['optional_options']=[
    make_option('-m','--mapping_fp', type='existing_filepath',
        help='path to category mapping file'),
    make_option('-c', '--category', type='string',
        help='name of the category over which to run the analysis'),
    make_option('-s', '--test', type="choice", 
        choices=CORRELATION_TEST_CHOICES.keys(),
        default='kendall', help='Test to use. Choices are:\n%s' % \
            (', '.join(CORRELATION_TEST_CHOICES.keys()))+'\n\t' + \
            '[default: %default]'),
    make_option('--individual_column', type='string', default=None, 
        help='Column header in mapping file that designates which sample '+\
            'is from which individual.'),
    make_option('--paired_t_fp', type='existing_filepath', default=None, 
        help='Pass a paired sample map as described in help to test with a '+\
            'paired_t_two_sample test. Overrides all other options. A '+\
            'paired sample map must be two columns without header that are '+\
            'tab separated. Each row contains samples which should be paired.'),
    make_option('--permutations', default=1000, type=int, 
        help='Number of permutations to use for bootstrapped tests.')]

script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    bt = parse_biom_table(open(opts.otu_table_fp))

    if opts.paired_t_fp is not None: #user wants to conduct paired t_test
        o = open(opts.paired_t_fp, 'U')
        lines = o.readlines()
        o.close()
        b_samples = []
        a_samples = []
        for i in lines:
            a,b = i.strip().split('\t')
            a_samples.append(a)
            b_samples.append(b)
        data_feed = paired_t_generator(bt, b_samples, a_samples)
        test_stats, pvals = run_paired_t(data_feed)
        # calculate corrected pvals
        fdr_pvals = array(benjamini_hochberg_step_down(pvals))
        bon_pvals = bonferroni_correction(pvals)
        # write output results after sorting
        lines = paired_t_output_formatter(bt, test_stats, pvals, fdr_pvals, 
            bon_pvals)
        lines = sort_by_pval(lines, ind=2)
        o = open(opts.output_fp, 'w')
        o.writelines('\n'.join(lines))
        o.close()
    elif opts.individual_column is not None: #user wants longitudinal corr
        pmf, _ = parse_mapping_file_to_dict(opts.mapping_fp)
        pmf, bt = sync_biom_and_mf(pmf, bt)
        sample_to_hsid = get_sample_cats(pmf, opts.individual_column)
        hsid_to_samples = get_cat_sample_groups(sample_to_hsid)
        hsid_to_sample_indices = get_sample_indices(hsid_to_samples, bt)
        data_feed = longitudinal_row_generator(bt, pmf, opts.category, 
            hsid_to_samples, hsid_to_sample_indices)
        rs, combo_pvals, combo_rhos, homogenous = \
            run_longitudinal_correlation_test(data_feed, opts.test, 
                CORRELATION_TEST_CHOICES)
        # make corrections
        fdr_ps = array(benjamini_hochberg_step_down(combo_pvals))
        fdr_ps = where(fdr_ps>1.0, 1.0, fdr_ps)
        bon_ps = array(bonferroni_correction(combo_pvals))
        bon_ps = where(bon_ps>1.0, 1.0, bon_ps)
        lines = longitudinal_correlation_formatter(bt, combo_rhos, combo_pvals, 
            homogenous, fdr_ps, bon_ps, rs, hsid_to_samples)
        # arange by fdr_ps
        lines = sort_by_pval(lines, 4)
        o = open(opts.output_fp, 'w')
        o.writelines('\n'.join(lines))
        o.close()
    else: #simple correlation analysis requested
        pmf, _ = parse_mapping_file_to_dict(opts.mapping_fp)
        pmf, bt = sync_biom_and_mf(pmf, bt)
        data_feed = correlation_row_generator(bt, pmf, opts.category, 
            opts.individual_column)
        corr_coefs, p_pvals, np_pvals, ci_highs, ci_lows = \
            run_correlation_test(data_feed, opts.test, CORRELATION_TEST_CHOICES)
        # calculate corrected pvals for both parametric and non-parametric 
        p_pvals_fdr = array(benjamini_hochberg_step_down(p_pvals))
        p_pvals_bon = bonferroni_correction(p_pvals)
        np_pvals_fdr = array(benjamini_hochberg_step_down(np_pvals))
        np_pvals_bon = bonferroni_correction(np_pvals)
        # correct for cases where values above 1.0 due to correction
        p_pvals_fdr = where(p_pvals_fdr>1.0, 1.0, p_pvals_fdr)
        p_pvals_bon = where(p_pvals_bon>1.0, 1.0, p_pvals_bon)
        np_pvals_fdr = where(np_pvals_fdr>1.0, 1.0, np_pvals_fdr)
        np_pvals_bon = where(np_pvals_bon>1.0, 1.0, np_pvals_bon)
        # write output results after sorting
        lines = correlation_output_formatter(bt, corr_coefs, p_pvals,
            p_pvals_fdr, p_pvals_bon, np_pvals, np_pvals_fdr, np_pvals_bon, 
            ci_highs, ci_lows)
        lines = sort_by_pval(lines, ind=2)
        o = open(opts.output_fp, 'w')
        o.writelines('\n'.join(lines))
        o.close()

if __name__ == "__main__":
    main()

