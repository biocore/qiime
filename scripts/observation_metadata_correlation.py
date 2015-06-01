#!/usr/bin/env python

__author__ = "Will Van Treuren"
__copyright__ = "Copyright 2014, The QIIME project"
__credits__ = ["Will Van Treuren", "Luke Ursell", "Catherine Lozupone",
               "Jesse Stombaugh", "Doug Wendel", "Dan Knights", "Greg Caporaso",
               "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
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

cscore_error_text = '''The only supported metric for p-value assignment with the
C-score is bootstrapping. For more information on the C-score, read Stone and
Roberts 1990 Oecologea paper 85: 74-79.'''

correlation_assignment_choices = ['spearman', 'pearson', 'kendall', 'cscore']
pvalue_assignment_choices = ['fisher_z_transform', 'parametric_t_distribution',
                            'bootstrapped', 'kendall']
bootstrap_functions = {'spearman': spearman, 'pearson': pearson,
                       'kendall': kendall, 'cscore': cscore}

script_info = {}
script_info['brief_description'] = """Correlation between observation abundances and continuous-valued metadata"""
script_info['script_description'] = """
This script calculates correlations between feature (aka observation) abundances
(relative or absolute) and numeric metadata. Several methods are provided to
allow the user to correlate features to sample metadata values including
Spearman's Rho, Pearson, Kendall's Tau, and the C or checkerboard score.
References for these methods are numerous, but good descriptions may be found in
'Biometry' by Sokal and Rolhf. A brief description of the available tests
follows:

- Pearson score: The Pearson score, aka Pearson's Product Moment correlation, is
  a scaled measure of the degree to which two sequences of numbers co-vary. For
  'correlated' sequences, Pearson > 0, and for 'anticorrelated' sequences
  Pearson < 0 (uncorrelated implies Pearson = 0). Pearson is a paramateric
  and linear measure of correlation.

- Spearman's Rho: The Spearman correlation is a non-paramateric measure of
  correlation between two sequences of numbers. Spearman correlation is
  appropriate for data where the values of the observations are not necessarily
  accurate, but for which their relative magnitudes are (see Biometry for more
  details).

- Kendall's Tau: Kendall's Tau is an alternative method of calculating
  correlation between two sequences of numbers. It is slower and less widely
  utilized than Spearman or Pearson scores.

- Cscore: The c-score or 'checkerboard score' is a measure of covariation
  between two sequences that is derived from traditional ecology (Stone and
  Roberts. 1990, Oecologia 85:74-79).

Raw correlation statistics alone reflect only the degree of association between
two sequences of numbers or vectors. Assigning a likelihood to these score via
a p-value can be done with several methods depending on the particular
assumptions that are used. This script allows four methods for calculating
p-values:

- Bootrapping: Bootstrapping is the most robust, but slowest procedure for
  calculating the p-value of a given correlation score. Bootstrapping takes the
  input sequences, shuffles the order of one, and then recomputes the
  correlation score. The p-value is then the number of times (out of the given
  number of permutations) that the score of the permuted sequence pair was more
  extreme than the observed pair. Bootstrapping is good when the underlying
  properties of the distributions are unknown.

- Parametric t distribution: the traditional method for calculating the
  significance of a correlation score, this method assumes that the scores are
  normally distributed and computes a t statistic for each correlation score in
  conjunction with the length of the sequences being correlated.

- Fisher Z transform: Fisher's Z transform is a way to make the distribution of
  correlation scores (especially when there are many highly correlated scores)
  look more normal. It is not to be confused with Fisher's transformation for
  the combination of p-values.

- Kendall's Tau: for the Kendall's Tau calculation, the specific Kendall's Tau
  p-value is provided.

Notes:

- The only supported metric for p-value assignment with the C-score is
  bootstrapping. For more information on the C-score, read Stone and Roberts
  1990 Oecologea paper 85: 74-79. If you don't pass
  pval_assignment_method='bootstrapped' while you have -s cscore, the script
  will raise an error.

- Assigning p-values to Kendall's Tau scores with the bootstrapping method is
  very slow.

"""

script_info['script_usage'] = []
script_info['script_usage'].append(
    ("Example 1:",
     "Calculate the correlation between OTUs in the table and the pH of the samples from whence they came:",
     "%prog -i otu_table.biom -m map.txt -c pH -s spearman -o spearman_otu_gradient.txt"))
script_info['script_usage'].append(
    ("Example 2:",
     "Calculate the correlation between OTUs in the table and the pH of the samples from whence they came using bootstrapping and pearson correlation:",
     "%prog -i otu_table.biom -m map.txt -c pH -s pearson --pval_assignment_method bootstrapped --permutations 100 -o pearson_bootstrapped.txt"))

script_info['output_description']= """
The output will be a tab-delimited file with the following headers. Each row
will record the values calculated for a given feature:

- Feature ID: ID of the features being correlated. These are the observation IDs
  in the BIOM table.
- Test stat.: the value of the test statistic for the given test.
- pval: the raw p-value returned by the given test.
- pval_fdr: the p-value corrected by the Benjamini-Hochberg FDR procedure for
  multiple comparisons.
- pval_bon: the p-value corrected by the Bonferroni procedure for multiple
  comparisons.
- [metadata]: this column will be present only if the BIOM table contained
  metadata information for your features. For example, if these are OTUs, and
  taxonomy is present in the BIOM table, this column will contain OTU
  taxonomy and will be named 'taxonomy'.

"""
script_info['required_options']=[
    make_option('-i','--otu_table_fp',
        help='path to input BIOM table',
        type='existing_path'),
    make_option('-o', '--output_fp', type='new_filepath',
        help='path to the output file to be created'),
    make_option('-m','--mapping_fp', type='existing_filepath',
        help='path to metadata mapping file'),
    make_option('-c', '--category',
        help='name of the category in the metadata mapping file over which to '
        'run the analysis')]

script_info['optional_options']=[
    make_option('-s', '--test', type="choice",
        choices=correlation_assignment_choices,
        default='spearman', help='Correlation method to use. Choices are: %s' %
        (', '.join(correlation_assignment_choices)) + ' [default: %default]'),
    make_option('--pval_assignment_method', type="choice",
        choices=pvalue_assignment_choices,
        default='fisher_z_transform', help='p-value method to use. Choices are: '
        '%s' % (', '.join(pvalue_assignment_choices)) + ' [default: %default]'),
    make_option('--metadata_key', default='taxonomy',
        help='Key to extract metadata from BIOM table. [default: %default]'),
    make_option('--permutations', default=1000, type=int,
        help='Number of permutations to use for bootstrapped p-value '
        'calculations. [default: %default]')]

script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if opts.test == 'cscore' and opts.pval_assignment_method != 'bootstrapped':
        option_parser.error(cscore_error_text)

    bt = load_table(opts.otu_table_fp)
    pmf, _ = parse_mapping_file_to_dict(opts.mapping_fp)

    samples_to_correlate = []
    md_values_to_correlate = []
    bt_sample_ids = bt.ids(axis='sample')

    for sample_id, sample_md in pmf.items():
        if sample_id in bt_sample_ids:
            try:
                v = is_computable_float(sample_md[opts.category])
                samples_to_correlate.append(sample_id)
                md_values_to_correlate.append(v)
            except KeyError:
                option_parser.error('The category (%s)' % opts.category +
                    ' was not found in the mapping file.')
            except ValueError:
                pass  # value couldn't be converted to float, ignore this sample
        else:
            pass  # sample in mf, but not bt

    # remove samples which are not found in the mapping file or do not have
    # metadata that converts to float
    bt.filter(ids_to_keep = samples_to_correlate, axis='sample')

    # sort the biom table so that feature values are retrieved in the same
    # order as the metadata in the samples they correspond to
    bt = bt.sort(sort_f = lambda _: samples_to_correlate, axis='sample')

    if bt.shape[1] <= 3:
        option_parser.error(filtration_error_text)

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
    o.write('\n')
    o.close()

if __name__ == "__main__":
    main()
