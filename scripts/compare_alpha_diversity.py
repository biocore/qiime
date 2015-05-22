#!/usr/bin/env python
# File created on 06 Jun 2011
from __future__ import division

__author__ = "William Van Treuren"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["William Van Treuren", "Greg Caparaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "William Van Treuren"
__email__ = "vantreur@colorado.edu"

import os
from os.path import join
from qiime.util import (parse_command_line_parameters,
                        make_option,
                        create_dir)
from qiime.compare_alpha_diversity import (compare_alpha_diversities,
                                           _correct_compare_alpha_results,
                                           test_types,
                                           correction_types,
                                           generate_alpha_diversity_boxplots)

script_info = {}
script_info[
    'brief_description'] = """This script compares alpha diversities based on a two-sample t-test using either parametric or non-parametric (Monte Carlo) methods."""

script_info['script_description'] = """
This script compares the alpha diversity of samples found in a collated alpha
diversity file. The comparison is done not between samples, but between groups
of samples. The groupings are created via the input category passed via
-c/--category. Any samples which have the same value under the catgory will be
grouped.

For example, if your mapping file had a category called 'Treatment' that
separated your samples into three groups (Treatment='Control', Treatment='Drug',
Treatment='2xDose'), passing 'Treatment' to this script would cause it to
compare (Control,Drug), (Control,2xDose), (2xDose, Drug) alpha diversity
values. By default the two-sample t-test will be nonparametric (i.e. using
Monte Carlo permutations to calculate the p-value), though the user has the
option to make the test a parametric t-test.

The script creates an output file in tab-separated format where each row is a
different group comparison. The columns in each row denote which two groups of
samples are being compared, as well as the mean and standard deviation of each
group's alpha diversity. Finally, the t-statistic and p-value are reported for
the comparison. This file can be most easily viewed in a spreadsheet program
such as Excel.

Note: Any iterations of a rarefaction at a given depth will be averaged. For
instance, if your collated_alpha file had 10 iterations of the rarefaction at
depth 480, the scores for the alpha diversity metrics of those 10 iterations
would be averaged (within sample). The iterations are not controlled by this
script; when multiple_rarefactions.py is called, the -n option specifies the
number of iterations that have occurred. The multiple comparison correction
takes into account the number of between group comparisons. If you do not know
the rarefaction depth available or you want to use the deepest rarefaction
level available then do not pass -d/--depth and it will default to using the
deepest available.

If t-statistics and/or p-values are None for any of your comparisons, there are
three possible reasons. The first is that there were undefined values in your
collated alpha diversity input file. This occurs if there were too few
sequences in one or more of the samples in the groups involved in those
comparisons to compute alpha diversity at that depth. You can either rerun
%prog passing a lower value for --depth, or you can re-run alpha diversity
after filtering samples with too few sequences. The second is that you had some
comparison where each treatment was represented by only a single sample. It is
not possible to perform a two-sample t-test on two samples each of length 1, so
None will be reported instead. The third possibility occurs when using the
nonparamteric t-test with small datasets where the Monte Carlo permutations
don't return a p-value because the distribution of the data has no variance.
The multiple comparisons correction will not penalize you for comparisons that
return as None regardless of origin.

If the means/standard deviations are None for any treatment group, the likely
cause is that there is an \'n\\a\' value in the collated_alpha file that was
passed.
"""

script_info['script_usage'] = []

script_info['script_usage'].append(("Comparing alpha diversities",
                                    "The following command takes the following input: a mapping file (which "
                                    "associaties each sample with a number of characteristics), alpha diversity "
                                    "metric (the results of collate_alpha for an alpha diverity metric, like "
                                    "PD_whole_tree), depth (the rarefaction depth to use for comparison), "
                                    "category (the category in the mapping file to determine which samples to "
                                    "compare to each other), and output filepath (a path to the output file to be created). A "
                                    "nonparametric two sample t-test is run to compare the alpha diversities "
                                    "using the default number of Monte Carlo permutations (999).",
                                    "%prog -i PD_whole_tree.txt -m mapping.txt -c Treatment -d 100 -o Treatment_PD100"))

script_info['script_usage'].append(("Comparing alpha diversities",
                                    "Similar to above, but performs comparisons for two categories.",
                                    "%prog -i PD_whole_tree.txt -m mapping.txt -c Treatment,DOB -d 100 -o Treatment_DOB_PD100"))

script_info['script_usage'].append(("Parametric t-test",
                                    "The following command runs a parametric two sample t-test using the "
                                    "t-distribution instead of Monte Carlo permutations at rarefaction depth 100.",
                                    "%prog -i PD_whole_tree.txt -m mapping.txt -c Treatment -d 100 -o "
                                    "PD_d100_parametric -t parametric"))

script_info['script_usage'].append(("Parametric t-test",
                                    "The following command runs a parametric two sample t-test using the "
                                    "t-distribution instead of Monte Carlo permutations at the greatest depth available.",
                                    "%prog -i PD_whole_tree.txt -m mapping.txt -c Treatment -o "
                                    "PD_dmax_parametric -t parametric"))

script_info['output_description'] = """
Generates a tsv stats file and pdf of boxplots for each input category.
Each row in the tsv file corresponds to a comparison between two groups of treatment values,
and includes the means and standard deviations of the two groups' alpha diversities,
along with the results of the two-sample t-test.
"""

script_info[
    'script_usage_output_to_remove'] = [
    '$PWD/PD_dmax_parametric.txt',
    '$PWD/PD_d100_parametric.txt',
    '$PWD/PD_d100.txt']

script_info['required_options'] = [
    make_option('-i',
                '--alpha_diversity_fp',
                action='store',
                type='existing_filepath',
                help='path to collated alpha diversity file (as generated by '
                'collate_alpha.py) [REQUIRED]'),
    make_option('-m',
                '--mapping_fp',
                action='store',
                type='existing_filepath',
                help='path to the mapping file [REQUIRED]'),
    make_option('-c',
                '--categories',
                action='store',
                type='string',
                help='comma-separated list of categories for comparison [REQUIRED]'),
    make_option('-o',
                '--output_dir',
                action='store',
                type='new_dirpath',
                help='directory where output files should be stored [REQUIRED]')]

script_info['optional_options'] = [
    make_option('-t', '--test_type', type='choice', choices=test_types,
                help='the type of test to perform when calculating the p-values. Valid '
                'choices: ' + ', '.join(test_types) + '. If test_type is '
                'nonparametric, Monte Carlo permutations will be used to determine the '
                'p-value. If test_type is parametric, the num_permutations option will '
                'be ignored and the t-distribution will be used instead [default: '
                '%default]', default='nonparametric'),
    make_option('-n', '--num_permutations', type='int', default=999,
                help='the number of permutations to perform when calculating the '
                'p-value. Must be greater than 10. Only applies if test_type is '
                'nonparametric [default: %default]'),
    make_option('-p', '--correction_method', type='choice',
                choices=correction_types, help='method to use for correcting multiple '
                'comparisons. Available methods are bonferroni, fdr, or none. '
                '[default: %default]', default='bonferroni'),
    make_option('-d', '--depth', type='int', default=None,
                help='depth of rarefaction file to use [default: greatest depth]')]


script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    mapping_fp = opts.mapping_fp
    alpha_diversity_fp = opts.alpha_diversity_fp
    categories = opts.categories.split(',')
    depth = opts.depth
    output_dir = opts.output_dir
    correction_method = opts.correction_method
    test_type = opts.test_type
    num_permutations = opts.num_permutations

    if num_permutations < 10:
        option_parser.error('Number of permuations must be greater than or '
                            'equal to 10.')

    create_dir(output_dir)
    for category in categories:
        stat_output_fp = join(output_dir, '%s_stats.txt' % category)
        boxplot_output_fp = join(output_dir, '%s_boxplots.pdf' % category)

        alpha_diversity_f = open(alpha_diversity_fp, 'U')
        mapping_f = open(mapping_fp, 'U')
        ttest_result, alphadiv_avgs = \
            compare_alpha_diversities(alpha_diversity_f,
                                      mapping_f,
                                      category,
                                      depth,
                                      test_type,
                                      num_permutations)
        alpha_diversity_f.close()
        mapping_f.close()

        corrected_result = _correct_compare_alpha_results(ttest_result,
                                                          correction_method)

        # write stats results
        stat_output_f = open(stat_output_fp, 'w')
        header = ('Group1\tGroup2\tGroup1 mean\tGroup1 std\tGroup2 mean\t'
                  'Group2 std\tt stat\tp-value')
        lines = [header]
        for (t0, t1), v in corrected_result.items():
            lines.append('\t'.join(map(str, [t0,
                                             t1,
                                             alphadiv_avgs[t0][0],
                                             alphadiv_avgs[t0][1],
                                             alphadiv_avgs[t1][0],
                                             alphadiv_avgs[t1][1],
                                             v[0],
                                             v[1]])))
        stat_output_f.write('\n'.join(lines) + '\n')
        stat_output_f.close()

        # write box plots
        alpha_diversity_f = open(alpha_diversity_fp, 'U')
        mapping_f = open(mapping_fp, 'U')
        boxplot = generate_alpha_diversity_boxplots(alpha_diversity_f,
                                                    mapping_f,
                                                    category,
                                                    depth)
        alpha_diversity_f.close()
        mapping_f.close()
        boxplot.savefig(boxplot_output_fp)

if __name__ == "__main__":
    main()
