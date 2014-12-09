#!/usr/bin/env python

__author__ = "Will Van Treuren"
__copyright__ = "Copyright 2014, The QIIME project"
__credits__ = ["Will Van Treuren", "Luke Ursell", "Catherine Lozupone",
               "Jesse Stombaugh", "Doug Wendel", "Dan Knights", "Greg Caporaso", 
               "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Will Van Treuren"
__email__ = "wdwvt1@gmail.com"

from qiime.util import (parse_command_line_parameters, make_option,
                        sync_biom_and_mf)
from qiime.stats import (benjamini_hochberg_step_down, bonferroni_correction,
                         assign_correlation_pval, correlate, pearson,
                         spearman, kendall, cscore)
from qiime.otu_significance import (correlate_output_formatter, sort_by_pval,
                                    run_paired_t, is_computable_float)
from qiime.parse import parse_mapping_file_to_dict
from biom import load_table
from numpy import array, where

filtration_error_text = '''The biom table does not have enough samples after
filtration to perform correlations (it has 3 or fewer). The filtration steps
that have occurred are:
1. Removing samples that were not found in the mapping file. 
2. Removing samples that were found in the mapping file but whose metadata for
the given category was not convertable to float.'''

cscore_error_text = '''The only supported metric for P-value assignment with the
C-score is bootstrapping. For more information on the C-score, read Stone and
Roberts 1990 Oecologea paper 85: 74-79.'''

correlation_assignment_choices = ['spearman', 'pearson', 'kendall', 'cscore']
pvalue_assignment_choices = ['fisher_z_transform', 'parametric_t_distribution',
                            'bootstrapped', 'kendall']
bootstrap_functions = {'spearman': spearman, 'pearson': pearson,
                       'kendall': kendall, 'cscore': cscore}

script_info = {}
script_info['brief_description'] = """This script calculates correlations between feature abundances and continuous-valued metadata."""
script_info['script_description'] = """
This script allows the calculation of:
    1. Correlations between feature abundances (relative or absolute) and 
       numeric metadata.
    2. Paired t-tests between two groups of samples.

Several methods are provided to allow the user to correlate 
features to sample metadata values including Spearmans Rho, Pearson, Kendall's
Tau, and the C or checkerboard score. 

The available methods for assigning p-values to the calculated correlation
scores are bootstrapping, Fisher's Z transformation, a parametric
t-distribution, and a Kendall's Tau specific p-value calculation.

This script also allows paired t testing. The paired t test is accomplished by
passing a paired mapping file which is just a two column (tab separation) table
with the samples that should be paired in each row. It should not have a header.

Notes:

The only supported metric for P-value assignment with the C-score is 
bootstrapping. For more information on the C-score, read Stone and Roberts 1990
Oecologea paper 85: 74-79. If you fail to pass 
pval_assignment_method='bootstrapped' while you have -s cscore, the script will 
error. 

Assigning pvalues to Kendall's Tau scores with the bootstrapping method is 
very slow.

"""

script_info['script_usage'] = []
script_info['script_usage'].append(
    ("Calculate the correlation between OTUs in the table and the pH of the samples from mich they came:",
     "",
     "%prog -i otu_table.biom -m map.txt -c pH -s spearman -o spearman_otu_gradient.txt"))
script_info['script_usage'].append(
    ("Calculate paired t values for a before and after group of samples:",
     "",
     "%prog -i otu_table.biom --paired_t_fp=paired_samples.txt -o paired.txt"))
script_info['script_usage'].append(
    ("Calculate the correlation between OTUs in the table and the pH of the samples from mich they came using bootstrapping and pearon correlation:",
     "",
     "%prog -i otu_table.biom -m map.txt -c pH -s pearson --pval_assignment_method bootstrapped --permutations 100 -o pearson_bootstrapped.txt"))

script_info['output_description']= """
The output will be a tab delimited file with the following headers. Each row
will record the values calculated for a fiven featue:
- Feature ID: ID of the features being correlated. If these are OTUs, then they
  will take the form of GreenGenes identifiers or de-novo identifiers. 
- Test-Statistic: the value of the test statistic for the given test.
- P: the raw P value returned by the given test. 
- FDR_P: the P value corrected by the Benjamini-Hochberg FDR procedure for 
  multiple comparisons.
- Bonferroni_P: the P value corrected by the Bonferroni procedure for multiple
  comparisons.
- Metadata - this column will be present only if the biom table contained
  metadata information for your features. If these are OTUs, and taxonomy is
  present in the biom table, this category will contain that taxonomy (or other
  metadata).
"""
script_info['required_options']=[
    make_option('-i','--otu_table_fp',
        help='path to input biom format table',
        type='existing_path'),
    make_option('-o', '--output_fp', type='new_filepath',
        help='path to the output file to be created')]

script_info['optional_options']=[
    make_option('-m','--mapping_fp', type='existing_filepath',
        help='path to category mapping file'),
    make_option('-c', '--category', type='string',
        help='name of the category over which to run the analysis'),
    make_option('-s', '--test', type="choice", 
        choices=correlation_assignment_choices,
        default='spearman', help='Correlation method to use. Choices are: %s' %
        (', '.join(correlation_assignment_choices)) + ' [default: %default]'),
    make_option('--pval_assignment_method', type="choice", 
        choices=pvalue_assignment_choices,
        default='fisher_z_transform', help='Pvalue method to use. Choices are: '
        '%s' % (', '.join(pvalue_assignment_choices)) + ' [default: %default]'),
    make_option('--metadata_key', default='taxonomy', type=str, 
        help='Key to extract metadata from biom table. [default: %default]'),
    make_option('--paired_t_fp', type='existing_filepath', default=None, 
        help='Pass a paired sample map as described in help to test with a '
            'paired_t_two_sample test. Overrides all other options. A '
            'paired sample map must be two columns without header that are '
            'tab separated. Each row contains samples which should be paired.'
            ' [default: %default]'),
    make_option('--permutations', default=1000, type=int, 
        help='Number of permutations to use for bootstrapped tests.'
            ' [default: %default]')]

script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if opts.test == 'cscore' and opts.pval_assignment_method != 'bootstrapped':
        raise ValueError(cscore_error_text)
    
    bt = load_table(opts.otu_table_fp)

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
        test_stats, pvals = run_paired_t(bt, a_samples, b_samples)
        # calculate corrected pvals
        fdr_pvals = array(benjamini_hochberg_step_down(pvals))
        bon_pvals = bonferroni_correction(pvals)
        # correct for cases where values above 1.0 due to correction
        fdr_pvals = where(array(fdr_pvals) > 1.0, 1.0, fdr_pvals)
        bon_pvals = where(array(bon_pvals) > 1.0, 1.0, bon_pvals)
        # write output results after sorting
        lines = correlate_output_formatter(bt, test_stats, pvals, fdr_pvals, 
                                           bon_pvals, md_key=opts.metadata_key)
        lines = sort_by_pval(lines, ind=2)
        o = open(opts.output_fp, 'w')
        o.writelines('\n'.join(lines))
        o.close()

    else:  # user wants normal correlation analysis
        pmf, _ = parse_mapping_file_to_dict(opts.mapping_fp)
        category = opts.category

        samples_to_correlate = []
        md_values_to_correlate = []
        bt_sample_ids = bt.ids(axis='sample')

        for sample_id, sample_md in pmf.items():
            if sample_id in bt_sample_ids:
                try:
                    v = is_computable_float(sample_md[category])
                    samples_to_correlate.append(sample_id)
                    md_values_to_correlate.append(v)
                except KeyError:
                    raise ValueError('The category (%s)' % opts.category +
                        ' was not found in the mapping file.')
            else:
                pass  # sample in mf, but not bt

        # remove samples which are not found in the mapping file or do not have
        # metadata that converts to float
        bt.filter(ids_to_keep = samples_to_correlate, axis='sample')

        # sort the biom table so that feature values are retrieved in the same 
        # order as the metadata in the samples they correspond to
        bt.sort(sort_f = lambda _: samples_to_correlate, axis='sample')

        if bt.shape[1] <= 3:
            raise ValueError(filtration_error_text)

        rhos = []
        pvals = []
        for feature_vector in bt.iter_data(axis='observation'):
            rho = correlate(feature_vector, md_values_to_correlate,
                            method=opts.test)
            pval = assign_correlation_pval(rho, len(feature_vector),
                                           method=opts.pval_assignment_method,
                                           permutations=opts.permutations,
                                           perm_test_fn=\
                                                bootstrap_functions[opts.test],
                                           v1=feature_vector,
                                           v2=md_values_to_correlate)
            rhos.append(rho)
            pvals.append(pval)

        fdr_pvals = benjamini_hochberg_step_down(pvals)
        bon_pvals = bonferroni_correction(pvals)
        # correct for cases where values above 1.0 due to correction
        fdr_pvals = where(array(fdr_pvals) > 1.0, 1.0, fdr_pvals)
        bon_pvals = where(array(bon_pvals) > 1.0, 1.0, bon_pvals)

        lines = correlate_output_formatter(bt, rhos, pvals, fdr_pvals,
                                           bon_pvals, opts.metadata_key)
        lines = sort_by_pval(lines, ind=2)

        o = open(opts.output_fp, 'w')
        o.writelines('\n'.join(lines))
        o.close()

if __name__ == "__main__":
    main()


