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
from qiime.ocs import (sync_biom_and_mf, fdr_correction, bonferroni_correction, 
    sort_by_pval, run_correlation_test, correlation_row_generator)
from qiime.pycogent_backports.test import (pearson, spearman, 
    kendall_correlation)
from qiime.parse import parse_mapping_file_to_dict
from biom.parse import parse_biom_table
from numpy import array, where

correlation_test_choices = {'pearson': pearson, 'spearman': spearman,
    'kendall': kendall_correlation}

script_info = {}
script_info['brief_description'] = """
"""
script_info['script_description'] = """
"""
script_info['script_usage'] = []
script_info['script_usage'].append(("add", "me", "please"))
script_info['output_description']= """
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

    data_feed = correlation_row_generator(bt, pmf, opts.category)
    corr_coefs, p_pvals, np_pvals, ci_highs, ci_lows = \
        run_correlation_test(data_feed, opts.test, test_choices)
    # calculate corrected pvals
    fdr_pvals = array(fdr_correction(pvals))
    bon_pvals = bonferroni_correction(pvals)
    # correct for cases where values above 1.0 due to correction
    fdr_pvals = where(fdr_pvals>1.0, 1.0, fdr_pvals)
    bon_pvals = where(bon_pvals>1.0, 1.0, bon_pvals)
    # write output results after sorting
    lines = output_formatter(bt, test_stats, pvals, fdr_pvals, bon_pvals, means,
        cat_sam_indices)
    lines = sort_by_pval(lines)
    o = open(opts.output_fp, 'w')
    o.writelines('\n'.join(lines))
    o.close()

if __name__ == "__main__":
    main()

