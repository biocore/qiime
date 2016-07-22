#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from os.path import basename, join
from skbio.util import create_dir
from qiime.parse import parse_sample_id_map, parse_taxa_summary_table
from qiime.util import (add_filename_suffix, parse_command_line_parameters,
                        get_options_lookup, make_option)

from qiime.compare_taxa_summaries import (comparison_modes,
                                          compare_taxa_summaries, correlation_types, tail_types)

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = "Compares taxa summary files"
script_info['script_description'] = """
This script compares two taxa summary files by computing the correlation
coefficient between pairs of samples. This is useful, for example, if you want
to compare the taxonomic composition of mock communities that were assigned
using different taxonomy assigners in order to see if they are correlated or
not. Another example use-case is to compare the taxonomic composition of
several mock community replicate samples to a single expected, or known, sample
community.

This script is also useful for sorting and filling taxa summary files so that
each sample has the same taxa listed in the same order (with missing taxa
reporting an abundance of zero). The sorted and filled taxa summary files can
then be passed to a script, such as plot_taxa_summary.py, to visually compare
the differences using the same taxa coloring scheme.

For more information and examples pertaining to this script, please refer to
the accompanying tutorial, which can be found at
http://qiime.org/tutorials/taxa_summary_comparison.html.
"""

script_info['script_usage'] = []
script_info['script_usage'].append(("Paired sample comparison",
                                    "Compare all samples that have matching sample IDs between the two input taxa "
                                    "summary files using the pearson correlation coefficient. The first input "
                                    "taxa summary file is from the overview tutorial, using the RDP classifier "
                                    "with a confidence level of 0.60 and the gg_otus_4feb2011 97% representative "
                                    "set. The second input taxa summary file was generated the same way, except "
                                    "for using a confidence level of 0.80.",
                                    "%prog -i ts_rdp_0.60.txt,ts_rdp_0.80.txt -m paired -o taxa_comp"))

script_info['script_usage'].append(("Paired sample comparison with sample ID "
                                    "map", "Compare samples based on the mappings in the sample ID map using the "
                                    "spearman correlation coefficient. The second input taxa summary file is "
                                    "simply the original ts_rdp_0.60.txt file with all sample IDs containing "
                                    "'PC.' renamed to 'S.'.",
                                    "%prog -i ts_rdp_0.80.txt,ts_rdp_0.60_renamed.txt -m paired -o "
                                    "taxa_comp_using_sample_id_map -s sample_id_map.txt -c spearman"))

script_info['script_usage'].append(("Detailed paired sample comparison",
                                    "Compare all samples that have matching sample IDs between the two input taxa "
                                    "summary files using the pearson correlation coefficient. Additionally, "
                                    "compute the correlation coefficient between each pair of samples "
                                    "individually.",
                                    "%prog -i ts_rdp_0.60.txt,ts_rdp_0.80.txt -m paired -o taxa_comp_detailed "
                                    "--perform_detailed_comparisons"))

script_info['script_usage'].append(("One-tailed test",
                                    "Compare all samples that have matching sample IDs between the two input taxa "
                                    "summary files using the pearson correlation coefficient. Perform a "
                                    "one-tailed (negative association) test of significance for both parametric "
                                    "and nonparametric tests. Additionally, compute a 90% confidence interval for "
                                    "the correlation coefficient. Note that the confidence interval will still be "
                                    "two-sided.",
                                    "%prog -i ts_rdp_0.60.txt,ts_rdp_0.80.txt -m paired -o taxa_comp_one_tailed "
                                    "-t low -l 0.90"))

script_info['output_description'] = """
The script will always output at least three files to the specified output
directory. Two files will be the sorted and filled versions of the input taxa
summary files, which can then be used in plot_taxa_summary.py to visualize the
differences in taxonomic composition. These files will be named based on the
basename of the input files. If the input files' basenames are the same, the
output files will have '0' and '1' appended to their names to keep the
filenames unique. The first input taxa summary file will have '0' in its
filename and the second input taxa summary file will have '1' in its filename.

The third output file will contain the results of the overall comparison of the
input taxa summary files using the specified sample pairings. The correlation
coefficient, parametric p-value, nonparametric p-value, and a confidence
interval for the correlation coefficient will be included.

If --perform_detailed_comparisons is specified, the fourth output file is a
tab-separated file containing the correlation coefficients that were computed
between each of the paired samples. Each line will contain the sample IDs of
the samples that were compared, followed by the correlation coefficient that
was computed, followed by the parametric and nonparametric p-values
(uncorrrected and Bonferroni-corrected) and a confidence interval for the
correlation coefficient.

The output files will contain comments at the top explaining the types of tests
that were performed.
"""

script_info['required_options'] = [
    make_option('-i', '--taxa_summary_fps', type='existing_filepaths',
                help='the two input taxa summary filepaths, comma-separated. These '
                'will usually be the files that are output by summarize_taxa.py. '
                'These taxa summary files do not need to have the same taxa in the '
                'same order, as the script will make them compatible before comparing '
                'them'),
    options_lookup['output_dir'],
    make_option('-m', '--comparison_mode', type='choice',
                choices=comparison_modes, help='the type of comparison to '
                'perform. Valid choices: ' + ' or '.join(comparison_modes) +
                '. "paired" will compare each sample in the taxa summary '
                'files that match based on sample ID, or that match given a sample ID '
                'map (see the --sample_id_map_fp option for more information). '
                '"expected" will compare each sample in the first taxa summary file '
                'to an expected sample (contained in the second taxa summary file). '
                'If "expected", the second taxa summary file must contain only a '
                'single sample that all other samples will be compared to (unless the '
                '--expected_sample_id option is provided)')
]
script_info['optional_options'] = [
    make_option('-c', '--correlation_type', type='choice',
                choices=correlation_types, help='the type of correlation coefficient '
                'to compute. Valid choices: ' + ' or '.join(correlation_types) +
                ' [default: %default]', default='pearson'),
    make_option('-t', '--tail_type', type='choice',
                choices=tail_types, help='the type of tail test to compute when '
                'calculating the p-values. "high" specifies a one-tailed test for '
                'values greater than the observed correlation coefficient (positive '
                'association), while "low" specifies a one-tailed test for values '
                'less than the observed correlation coefficient (negative '
                'association). "two-sided" specifies a two-tailed test for values '
                'greater in magnitude than the observed correlation coefficient. '
                'Valid choices: ' +
                ' or '.join(tail_types) + ' [default: %default]',
                default='two-sided'),
    make_option('-n', '--num_permutations', type='int',
                help='the number of permutations to perform when calculating the '
                'nonparametric p-value. Must be an integer greater than or equal to '
                'zero. If zero, the nonparametric test of significance will not be '
                'performed and the nonparametric p-value will be reported as "N/A" '
                '[default: %default]', default=999),
    make_option('-l', '--confidence_level', type='float',
                help='the confidence level of the correlation coefficient confidence '
                'interval. Must be a value between 0 and 1 (exclusive). For example, '
                'a 95% confidence interval would be 0.95 [default: %default]',
                default=0.95),
    make_option('-s', '--sample_id_map_fp', type='existing_filepath',
                help='map of original sample IDs to new sample IDs. Use this to match '
                'up sample IDs that should be compared between the two taxa summary '
                'files. Each line should contain an original sample ID, a tab, and '
                'the new sample ID. All original sample IDs from the two input taxa '
                'summary files must be mapped. This option only applies if the '
                'comparison mode is "paired". If not provided, only sample IDs that '
                'exist in both taxa summary files will be compared '
                '[default: %default]', default=None),
    make_option('-e', '--expected_sample_id', type='string',
                help='the sample ID in the second "expected" taxa summary file to '
                'compare all samples to. This option only applies if the comparison '
                'mode is "expected". If not provided, the second taxa summary file '
                'must have only one sample [default: %default]', default=None),
    make_option('--perform_detailed_comparisons', action='store_true',
                help='Perform a comparison for each sample pair in addition to the '
                'single overall comparison. The results will include the '
                'Bonferroni-corrected p-values in addition to the original p-values '
                '[default: %default]', default=False)
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if len(opts.taxa_summary_fps) != 2:
        option_parser.error("Exactly two taxa summary files are required. You "
                            "provided %d." % len(opts.taxa_summary_fps))

    # Create the output dir if it doesn't already exist.
    try:
        create_dir(opts.output_dir)
    except:
        option_parser.error("Could not create or access output directory "
                            "specified with the -o option.")

    sample_id_map = None
    if opts.sample_id_map_fp:
        sample_id_map = parse_sample_id_map(open(opts.sample_id_map_fp, 'U'))

    results = compare_taxa_summaries(
        parse_taxa_summary_table(open(opts.taxa_summary_fps[0], 'U')),
        parse_taxa_summary_table(open(opts.taxa_summary_fps[1], 'U')),
        opts.comparison_mode, correlation_type=opts.correlation_type,
        tail_type=opts.tail_type, num_permutations=opts.num_permutations,
        confidence_level=opts.confidence_level,
        perform_detailed_comparisons=opts.perform_detailed_comparisons,
        sample_id_map=sample_id_map,
        expected_sample_id=opts.expected_sample_id)

    # Write out the sorted and filled taxa summaries, basing their
    # filenames on the original input filenames. If the filenames are the same,
    # append a number to each filename.
    same_filenames = False
    if basename(opts.taxa_summary_fps[0]) == \
       basename(opts.taxa_summary_fps[1]):
        same_filenames = True

    for orig_ts_fp, filled_ts_lines, file_num in zip(opts.taxa_summary_fps,
                                                     results[:2], range(0, 2)):
        filename_suffix = '_sorted_and_filled'
        if same_filenames:
            filename_suffix += '_%d' % file_num
        filled_ts_fp = add_filename_suffix(orig_ts_fp, filename_suffix)
        filled_ts_f = open(join(opts.output_dir, filled_ts_fp), 'w')
        filled_ts_f.write(filled_ts_lines)
        filled_ts_f.close()

    # Write the overall comparison result.
    overall_comp_f = open(join(opts.output_dir, 'overall_comparison.txt'), 'w')
    overall_comp_f.write(results[2])
    overall_comp_f.close()

    # Write the correlation vector containing the pairwise sample comparisons.
    if opts.perform_detailed_comparisons:
        corr_vec_f = open(join(opts.output_dir,
                               'detailed_comparisons.txt'), 'w')
        corr_vec_f.write(results[3])
        corr_vec_f.close()


if __name__ == "__main__":
    main()
