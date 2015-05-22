#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from os.path import join
from string import strip

from skbio.util import create_dir

from qiime.make_distance_boxplots import make_distance_boxplots, SORT_TYPES
from qiime.stats import all_pairs_t_test, tail_types
from qiime.util import (get_options_lookup, make_option,
                        parse_command_line_parameters)

script_info = {}
script_info[
    'brief_description'] = "Creates boxplots to compare distances between categories"
script_info['script_description'] = """
This script creates boxplots that allow for the comparison between different
categories found within the mapping file. The boxplots that are created compare
distances within all samples of a field value, as well as between different
field values. Individual within and between distances are also plotted.

The script also performs two-sample t-tests for all pairs of boxplots to help
determine which boxplots (distributions) are significantly different.

Tip: the script tries its best to fit everything into the plot, but there are
cases where plot elements may get cut off (e.g. if axis labels are extremely
long), or things may appear squashed, cluttered, or too small (e.g. if
there are many boxplots in one plot). Increasing the width and/or height of the
plot (using --width and --height) usually fixes these problems.

For more information and examples pertaining to this script, please refer to
the accompanying tutorial, which can be found at
http://qiime.org/tutorials/creating_distance_comparison_plots.html.
"""

script_info['script_usage'] = []
script_info['script_usage'].append((
    "Compare distances between Fast and Control samples",
    "This example will generate an image with boxplots for all within and all "
    "between distances for the field Treatment, and will also include plots for "
    "individual within (e.g. Control vs. Control, Fast vs. Fast) and individual "
    "between (e.g. Control vs. Fast). The generated plot PDF and signifiance "
    "testing results will be written to the output directory 'out1'.",
    "%prog -d unweighted_unifrac_dm.txt -m Fasting_Map.txt -f \"Treatment\" -o "
    "out1"))

script_info['script_usage'].append((
    "Only plot individual field value distances",
    "This example will generate a PNG of all individual field value distances "
    "(within and between) for the Treatment field.",
    "%prog -d unweighted_unifrac_dm.txt -m Fasting_Map.txt -f \"Treatment\" -o "
    "out2 -g png --suppress_all_within --suppress_all_between"))

script_info['script_usage'].append((
    "Save raw data",
    "This example will generate an SVG image of the boxplots and also output the "
    "plotting data to a tab-delimited file.",
    "%prog -d unweighted_unifrac_dm.txt -m Fasting_Map.txt -f \"Treatment\" -o "
    "out3 -g svg --save_raw_data"))

script_info['script_usage'].append((
    "Suppress significance tests",
    "This example will only generate a plot and skip the significance testing "
    "step. This can be useful if you are operating on a large dataset and are not "
    "interested in performing the statistical tests (or at least not initially).",
    "%prog -d unweighted_unifrac_dm.txt -m Fasting_Map.txt -f \"Treatment\" -o "
    "out4 --suppress_significance_tests"))

script_info['script_usage'].append((
    "Sort boxplots",
    "To sort the boxplots by increasing median, supply the --sort option.",
    "%prog -d unweighted_unifrac_dm.txt -m Fasting_Map.txt -f \"Treatment\" "
    "-o out5 --sort median"))

script_info['output_description'] = """
Images of the plots are written to the specified output directory (one image
per field). The raw data used in the plots and the results of significance
tests can optionally be written into tab-delimited files (one file per field)
that are most easily viewed in a spreadsheet program such as Microsoft Excel.
"""

options = get_options_lookup()

script_info['required_options'] = [
    options['mapping_fp'],
    options['output_dir'],
    make_option('-d', '--distance_matrix_fp',
                help='input distance matrix filepath (i.e. the result of '
                'beta_diversity.py). WARNING: Only symmetric, hollow distance '
                'matrices may be used as input. Asymmetric distance matrices, such as '
                'those obtained by the UniFrac Gain metric (i.e. beta_diversity.py '
                '-m unifrac_g), should not be used as input',
                type='existing_filepath'),
    make_option('-f', '--fields', type='string',
                help='comma-separated list of fields to compare, where the list of '
                'fields should be in quotes (e.g. "Field1,Field2,Field3")')]

script_info['optional_options'] = [
    make_option('-g', '--imagetype',
                help='type of image to produce (i.e. png, svg, pdf) '
                '[default: %default]', default='pdf', type="choice",
                choices=['pdf', 'png', 'svg']),
    make_option('--save_raw_data', action='store_true',
                help='store raw data used to create boxplots in tab-delimited files '
                '[default: %default]',
                default=False),
    make_option('--suppress_all_within', action='store_true',
                help='suppress plotting of "all within" boxplot [default: %default]',
                default=False),
    make_option('--suppress_all_between', action='store_true',
                help='suppress plotting of "all between" boxplot [default: %default]',
                default=False),
    make_option('--suppress_individual_within', action='store_true',
                help='suppress plotting of individual "within" boxplot(s) '
                '[default: %default]',
                default=False),
    make_option('--suppress_individual_between', action='store_true',
                help='suppress plotting of individual "between" boxplot(s) '
                '[default: %default]',
                default=False),
    make_option('--suppress_significance_tests', action='store_true',
                help='suppress performing signifance tests between each pair of '
                'boxplots [default: %default]', default=False),
    make_option('-n', '--num_permutations', type='int',
                help='the number of Monte Carlo permutations to perform when '
                'calculating the nonparametric p-value in the significance tests. '
                'Must be an integer greater than or equal to zero. If zero, the '
                'nonparametric p-value will not be calculated and will instead be '
                'reported as "N/A". This option has no effect if '
                '--suppress_significance_tests is supplied [default: %default]',
                default=0),
    make_option('-t', '--tail_type', type='choice',
                choices=tail_types, help='the type of tail test to compute when '
                'calculating the p-values in the significance tests. "high" specifies '
                'a one-tailed test for values greater than the observed t statistic, '
                'while "low" specifies a one-tailed test for values less than the '
                'observed t statistic. "two-sided" specifies a two-tailed test for '
                'values greater in magnitude than the observed t statistic. This '
                'option has no effect if --suppress_significance_tests is supplied. '
                'Valid choices: ' +
                ' or '.join(tail_types) + ' [default: %default]',
                default='two-sided'),
    make_option('--y_min',
                help='the minimum y-axis value in the resulting plot. If "auto", '
                'it is automatically calculated [default: %default]',
                default=0, type='string'),
    make_option('--y_max',
                help='the maximum y-axis value in the resulting plot. If "auto", '
                'it is automatically calculated [default: %default]',
                default=1, type='string'),
    make_option('--width',
                help='width of the output image in inches. If not provided, '
                'a "best guess" width will be used [default: auto]',
                default=None, type='float'),
    make_option('--height',
                help='height of the output image in inches [default: %default]',
                default=6, type='float'),
    make_option('--transparent', action='store_true',
                help='make output images transparent (useful for overlaying an image '
                'on top of a colored background) [default: %default]',
                default=False),
    make_option('--whisker_length',
                help='length of the whiskers as a function of the IQR. For example, '
                'if 1.5, the whiskers extend to 1.5 * IQR. Anything outside of '
                'that range is seen as an outlier [default: %default]',
                default='1.5', type='float'),
    make_option('--box_width',
                help='width of each box in plot units [default: %default]',
                default='0.5', type='float'),
    make_option('--box_color',
                help='the color of the boxes. Can be any valid matplotlib color '
                'string, such as "black", "magenta", "blue", etc. See '
                'http://matplotlib.sourceforge.net/api/colors_api.html for more '
                'examples of valid color strings that may be used. Will be ignored if '
                '--color_individual_within_by_field is supplied [default: '
                'same as plot background, which is white unless --transparent is '
                'enabled]',
                default=None, type='string'),
    make_option('--color_individual_within_by_field',
                help='field in the the mapping file to color the individual '
                '"within" boxes by. A legend will be provided to match boxplot colors '
                'to field states. A one-to-one mapping must exist between the field '
                'to be colored and the field to color by, otherwise the coloring will '
                'be ambiguous. If this option is supplied, --box_color will be '
                'ignored. If --suppress_individual_within is supplied, this option '
                'will be ignored [default: %default]',
                default=None, type='string'),
    make_option('--sort', type='choice', choices=SORT_TYPES,
                help='If "median", sort boxplots by increasing median. If '
                '"alphabetical", sort boxplots alphabetically by their '
                'labels. If this option is not supplied (the default), '
                'boxplots will be grouped logically as follows: all within, '
                'all between, individual within, and individual between '
                '[default: %default]', default=None)]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # Create the output dir if it doesn't already exist.
    out_dir = opts.output_dir
    try:
        create_dir(out_dir)
    except:
        option_parser.error("Could not create or access output directory "
                            "specified with the -o option.")

    map_f = open(opts.mapping_fp, 'U')
    dm_f = open(opts.distance_matrix_fp, 'U')

    fields = map(strip, opts.fields.split(','))
    fields = [field.strip('"').strip("'") for field in fields]

    color_individual_within_by_field = opts.color_individual_within_by_field
    results = make_distance_boxplots(dm_f, map_f, fields, width=opts.width,
                                     height=opts.height, suppress_all_within=opts.suppress_all_within,
                                     suppress_all_between=opts.suppress_all_between,
                                     suppress_individual_within=opts.suppress_individual_within,
                                     suppress_individual_between=opts.suppress_individual_between,
                                     y_min=opts.y_min, y_max=opts.y_max,
                                     whisker_length=opts.whisker_length, box_width=opts.box_width,
                                     box_color=opts.box_color,
                                     color_individual_within_by_field=color_individual_within_by_field,
                                     sort=opts.sort)

    for field, plot_figure, plot_data, plot_labels, plot_colors in results:
        output_plot_fp = join(out_dir, "%s_Distances.%s" %
                              (field, opts.imagetype))
        plot_figure.savefig(output_plot_fp, format=opts.imagetype,
                            transparent=opts.transparent)

        if not opts.suppress_significance_tests:
            sig_tests_f = open(join(out_dir, "%s_Stats.txt" % field), 'w')
            sig_tests_results = all_pairs_t_test(plot_labels, plot_data,
                                                 tail_type=opts.tail_type,
                                                 num_permutations=opts.num_permutations)
            sig_tests_f.write(sig_tests_results)
            sig_tests_f.close()

        if opts.save_raw_data:
            # Write the raw plot data into a tab-delimited file.
            assert(len(plot_labels) == len(plot_data))
            raw_data_fp = join(out_dir, "%s_Distances.txt" % field)
            raw_data_f = open(raw_data_fp, 'w')

            for label, data in zip(plot_labels, plot_data):
                raw_data_f.write(label.replace(" ", "_") + "\t")
                raw_data_f.write("\t".join(map(str, data)))
                raw_data_f.write("\n")
            raw_data_f.close()


if __name__ == "__main__":
    main()
