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
    paired_t_output_formatter, get_sample_cats, get_cat_sample_groups, 
    get_sample_indices, CORRELATION_PVALUE_CHOICES, 
    grouped_correlation_row_generator, run_grouped_correlation, 
    grouped_correlation_formatter)
from qiime.parse import parse_mapping_file_to_dict
from biom.parse import parse_biom_table
from numpy import array, where, logical_xor


script_info = {}
script_info['brief_description'] = """
This script allows the calculation of:
    1. Correlations between OTU abundances (relative or absolute) and numeric 
       metadata.
    2. Correlations between OTU abundances and numeric metadata taking into 
       account samples groupings based on passed categorical metadata. 
    3. Paired t-tests between two groups of samples.
The available methods for calculating correlations (1 and 2 above) are Spearmans
Rho, Pearson, Kendall's Tau, and the C or checkerboard score. The available 
methods for assigning p-values to the calculated correlation scores are 
bootstrapping, Fisher's Z transformation, a parametric t-distribution, and 
a Kendall's Tau specific p-value calculation.
"""
script_info['script_description'] = """
This script calculates the correlation between OTU values and a gradient of 
sample data. Several methods are provided to allow the user to correlate OTUs to
sample metadata values. Grouped correlations are also supported. 
The tests of OTU correlation to a metadata field are accomplished by passing a 
mapping file and a category from which to pull the data. If the data are not 
convertable to floats, the script will abort. 
The tests of OTU correlation to a metadata field with a grouped component 
are accomplished by passing the --grouping_category which tells the script which
samples should be formed in groups, along with the -c option which tells the 
script which metadata field/column to use as the gradient for the correlation. 
This will most frequently be useful when there is time series data
and you have samples from an individual at a number of time points.  
The paired t test is accomplished by passing a paired mapping file which is just
a two column (tab separation) table with the samples that should be paired in 
each row. It should not have a header.
The available tests are Kendall's Tau, Spearmans rank correlation, Pearsons 
product moment correlation, and the C-score (or checkerboard score). 
This script generates a tab separated output file which differs based on which 
test you have chosen. If you have chosen simple correlation or paired_t (i.e. 
you have not passed the --grouping_category option ) then you will see
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
--grouping_category) you will have:
OTU
Individual:X stat - correlation statistics from the given test. There will be as
 many of these headers as there are individuals in the individual column.
Taxonomy

Warnings:
The only supported metric for P-value assignment with the C-score is 
bootstrapping. For more information on the C-score, read Stone and Roberts 1990
Oecologea paper 85: 74-79. If you fail to pass 
pval_assignment_method='bootstrapped' while you have -s cscore, the script will 
error. 

Assigning pvalues to Kendall's Tau scores with the bootstrapping method is 
very slow.


"""
script_info['script_usage'] = []
script_info['script_usage'].append(("Calculate the correlation between OTUs in the table and the pH of the samples from mich they came:", "", "%prog -i otu_table.biom -m map.txt -c pH -s spearman -o spearman_otu_gradient.txt"))
script_info['script_usage'].append(("Calculate correlation between OTUs over the course of a treatment regimen assuming that 'hsid' is the column in the mapping file which identifies the individual subject from which each sample came and 'treatment_day' refers to the day after treatment was started for each individual:", "", "%prog -i otu_table.biom -m map.txt -c treatment_day -s kendall --grouping_category hsid -o kendall_longitudinal_otu_gradient.txt"))
script_info['script_usage'].append(("Calculate paired t values for a before and after group of samples:", "", "%prog -i otu_table.biom --paired_t_fp=paired_samples.txt -o kendall_longitudinal_otu_gradient.txt"))

script_info['output_description']= """
This script generates a tab separated output file which differs based on which 
test you have chosen. If you have chosen simple correlation or paired_t (i.e. 
you have not passed the --grouping_category option) then you will see
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
--grouping_category) you will have:
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
        default='spearman', help='Test to use. Choices are:\n%s' % \
            (', '.join(CORRELATION_TEST_CHOICES.keys()))+'\n\t' + \
            '[default: %default]'),
    make_option('--pval_assignment_method', type="choice", 
        choices=CORRELATION_PVALUE_CHOICES,
        default='fisher_z_transform', help='Test to use. Choices are:\n%s' % \
            (', '.join(CORRELATION_PVALUE_CHOICES))+'\n\t' + \
            '[default: %default]'),
    make_option('--metadata_key', default='taxonomy', type=str, 
        help='Key to extract metadata from biom table. default: %default]'),
    make_option('--grouping_category', type='string', default=None, 
        help='Column header in mapping file that designates which category '+\
            'to group samples by.'),
    make_option('--paired_t_fp', type='existing_filepath', default=None, 
        help='Pass a paired sample map as described in help to test with a '+\
            'paired_t_two_sample test. Overrides all other options. A '+\
            'paired sample map must be two columns without header that are '+\
            'tab separated. Each row contains samples which should be paired.'),
    make_option('--permutations', default=1000, type=int, 
        help='Number of permutations to use for bootstrapped tests.'+\
            '[default: %default]'),
    make_option('--biom_samples_are_superset', action='store_true', 
        default=False, 
        help='If this flag is passed you will be able to use a biom table '+\
            'that contains all the samples listed in the mapping file '+\
            'as well as additional samples not listed in the mapping file. '+\
            'Only their intersecting samples will be used for calculations.'),
    make_option('--print_non_overlap', action='store_true', default=False, 
        help='If this flag is passed the script will display the samples that'+\
            ' do not overlap between the mapping file and the biom file.')]

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
            bon_pvals, md_key=opts.metadata_key)
        lines = sort_by_pval(lines, ind=2)
        o = open(opts.output_fp, 'w')
        o.writelines('\n'.join(lines))
        o.close()
    elif opts.grouping_category is not None: #user wants grouped corr
        tmp_bt = parse_biom_table(open(opts.otu_table_fp))
        tmp_pmf, _ = parse_mapping_file_to_dict(opts.mapping_fp)
        pmf, bt, nonshared_samples = sync_biom_and_mf(tmp_pmf, tmp_bt)
        if not opts.biom_samples_are_superset: 
            # user indicates biom sample should be subset of mapping file samples
            if any([i in nonshared_samples for i in tmp_bt.SampleIds]):
                raise ValueError('The samples in the biom table are a' +\
                    ' superset of the samples in the mapping file. The' +\
                    ' script will abort in this case even though the' +\
                    ' calculations wouldn\'t be affected, to ensure' +\
                    ' consistency within QIIME. Pass the' +\
                    ' --biom_samples_are_superset option to disable this'+\
                    ' behavior.')
        if opts.print_non_overlap: #user wants non-overlapping samples printed out
            print 'The following samples were not shared between the mapping '+\
                'file and the biom file and will not be included in the '+\
                'analysis:\n'+' '.join(nonshared_samples)
        sample_to_group = get_sample_cats(pmf, opts.grouping_category)
        group_to_samples = get_cat_sample_groups(sample_to_group)
        category_values, metadata_values, otu_values = \
            grouped_correlation_row_generator(bt, pmf, opts.category,
                group_to_samples)
        rhos, pvals, f_pvals, f_rhos, f_hs = \
            run_grouped_correlation(metadata_values, otu_values, opts.test, 
                CORRELATION_TEST_CHOICES, opts.pval_assignment_method)
        # make lines
        lines = grouped_correlation_formatter(bt, rhos, pvals, f_rhos, f_pvals,
            f_hs, opts.grouping_category, category_values, opts.metadata_key)
        # arange by population correlation coefficient 
        lines = sort_by_pval(lines, ind=-4)
        o = open(opts.output_fp, 'w')
        o.writelines('\n'.join(lines))
        o.close()
    else: #simple correlation analysis requested
        if opts.category is None:
            option_parser.error('If simple correlation analysis is requested '+\
                ' a category must be passed with the -c option.')
        if logical_xor(opts.test=='kendall', 
            opts.pval_assignment_method=='kendall'):
            print ('Warning: you are using %s to calculate correlation and '+
                '%s to calculate the pvalue associated with that correlation '+
                'score. Using both -s kendall and --pval_assignment_method '+
                'kendall together or neither is advised.') % (opts.test, 
                opts.pval_assignment_method)
        if opts.test=='cscore' and opts.pval_assignment_method!='bootstrapped':
            option_parser.error('The only p-value assignment method supported '+
                'for the checkerboard score is bootstrapping.')
        pmf, _ = parse_mapping_file_to_dict(opts.mapping_fp)
        pmf, bt, nss = sync_biom_and_mf(pmf, bt)
        data_feed = correlation_row_generator(bt, pmf, opts.category)
        corr_coefs, pvals = run_correlation_test(data_feed, opts.test, 
            CORRELATION_TEST_CHOICES, opts.pval_assignment_method, 
            permutations=opts.permutations)
        # calculate corrected pvals for both parametric and non-parametric 
        pvals_fdr = array(benjamini_hochberg_step_down(pvals))
        pvals_bon = bonferroni_correction(pvals)
        # correct for cases where values above 1.0 due to correction
        pvals_fdr = where(pvals_fdr>1.0, 1.0, pvals_fdr)
        pvals_bon = where(pvals_bon>1.0, 1.0, pvals_bon)
        # write output results after sorting
        lines = correlation_output_formatter(bt, corr_coefs, pvals,
            pvals_fdr, pvals_bon, md_key=opts.metadata_key)
        lines = sort_by_pval(lines, ind=2)
        o = open(opts.output_fp, 'w')
        o.writelines('\n'.join(lines))
        o.close()

if __name__ == "__main__":
    main()

