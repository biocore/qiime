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
from skbio.draw import grouped_distributions
from qiime.colors import data_colors, data_color_order
from qiime.group import get_field_state_comparisons
from qiime.colors import matplotlib_rgb_color
from qiime.parse import (group_by_field, parse_distmat, parse_mapping_file,
                         QiimeParseError)
from qiime.stats import all_pairs_t_test, tail_types
from qiime.util import (get_options_lookup, make_option,
                        parse_command_line_parameters)

script_info = {}
script_info[
    'brief_description'] = "Creates plots comparing distances between sample groupings"

script_info['script_description'] = """
This script creates plots (bar charts, scatter plots, or box plots) that
allow for the comparison between samples grouped at different field states
of a mapping file field.

This script can work with any field in the mapping file, and it can compare
any number of field states to all other field states within that field.
This script may be especially useful for fields that represent a time series,
because a plot can be generated showing the distances between samples at
certain timepoints against all other timepoints.

For example, a time field might contain the values 1, 2, 3, 4, and 5, which
label samples that are from day 1, day 2, day 3, and so on. This time field
can be specified when the script is run, as well as the timepoint(s) to
compare to every other timepoint. For example, two comparison groups
might be timepoints 1 and 2. The resulting plot would contain timepoints for
days 3, 4, and 5 along the x-axis, and at each of those timepoints, the
distances between day 1 and that timepoint would be plotted, as well as the
distances between day 2 and the timepoint.

The script also performs two-sample t-tests for all pairs of distributions to
help determine which distributions are significantly different from each other.

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
    "Compare distances between Native and Input samples for each timepoint in the "
    "Time field",
    "This example will generate a PDF containing a bar chart with the distances "
    "between Native samples and every other timepoint, as well as the distances "
    "between Input samples and every other timepoint. The output image will be "
    "put in the 'out1' directory. For more details about this example input data, "
    "please refer to the accompanying tutorial.",
    "%prog -d forearm_only_unweighted_unifrac_dm.txt -m "
    "costello_timeseries_map.txt -f TIME_SINCE_TRANSPLANT -c \"Native,Input\" -o "
    "out1"))

script_info['output_description'] = """
An image of the plot is written to the specified output directory. The raw data
used in the plots and the results of significance tests can optionally be
written into tab-delimited files that are most easily viewed in a spreadsheet
program such as Microsoft Excel.
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
    make_option('-f', '--field', type='string',
                help='field in the mapping file to make comparisons on'),
    make_option('-c', '--comparison_groups', type='string',
                help='comma-separated list of field states to compare to every other '
                'field state, where the list of field states should be in quotes '
                '(e.g. "FieldState1,FieldState2,FieldState3")')]

script_info['optional_options'] = [
    make_option('-t', '--plot_type',
                help='type of plot to produce ("bar" is bar chart, "scatter" is '
                'scatter plot, and "box" is box plot) [default: %default]',
                default='bar', type='choice', choices=['bar', 'scatter', 'box']),
    make_option('-g', '--imagetype',
                help='type of image to produce (i.e. png, svg, pdf) '
                '[default: %default]', default='pdf', type="choice",
                choices=['pdf', 'png', 'svg']),
    make_option('--save_raw_data', action='store_true',
                help='store raw data used to create plot in a tab-delimited file '
                '[default: %default]',
                default=False),
    make_option('--suppress_significance_tests', action='store_true',
                help='suppress performing signifance tests between each pair of '
                'distributions [default: %default]', default=False),
    make_option('-n', '--num_permutations', type='int',
                help='the number of Monte Carlo permutations to perform when '
                'calculating the nonparametric p-value in the significance tests. '
                'Must be an integer greater than or equal to zero. If zero, the '
                'nonparametric p-value will not be calculated and will instead be '
                'reported as "N/A". This option has no effect if '
                '--suppress_significance_tests is supplied [default: %default]',
                default=0),
    make_option('--tail_type', type='choice',
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
    make_option('--width',
                help='width of the output image in inches [default: %default]',
                default=12, type='float'),
    make_option('--height',
                help='height of the output image in inches [default: %default]',
                default=6, type='float'),
    make_option('--x_tick_labels_orientation',
                help='type of orientation for x-axis tick labels [default: %default]',
                default='vertical', type='choice', choices=['vertical', 'horizontal']),
    make_option('-a', '--label_type',
                help='Label type ("numeric" or "categorical"). '
                'If the label type is defined as numeric, the x-axis will be '
                'scaled accordingly. Otherwise the x-values will treated '
                'categorically and will be evenly spaced [default: %default].',
                default='categorical',
                type='choice', choices=['categorical', 'numeric']),
    make_option('--y_min',
                help='the minimum y-axis value in the resulting plot. If "auto", '
                'it is automatically calculated [default: %default]',
                default=0, type='string'),
    make_option('--y_max',
                help='the maximum y-axis value in the resulting plot. If "auto", '
                'it is automatically calculated [default: %default]',
                default=1, type='string'),
    make_option('--transparent', action='store_true',
                help='make output images transparent (useful for overlaying an image '
                'on top of a colored background ) [default: %default]',
                default=False),
    make_option('--whisker_length',
                help='if --plot_type is "box", determines the length of the whiskers '
                'as a function of the IQR. For example, if 1.5, the whiskers '
                'extend to 1.5 * IQR. Anything outside of that range is seen as '
                'an outlier. If --plot_type is not "box", this option is ignored '
                '[default: %default]',
                default='1.5', type='float'),
    make_option('--error_bar_type',
                help='if --plot_type is "bar", determines the type of error bars to '
                'use. "stdv" is standard deviation and "sem" is the standard '
                'error of the mean. If --plot_type is not "bar", this option is '
                'ignored [default: %default]',
                default='stdv', type='choice', choices=['stdv', 'sem']),
    make_option('--distribution_width',
                help='width (in plot units) of each individual distribution (e.g. each '
                'bar if the plot type is a bar chart, or the width of each box '
                'if the plot type is a boxplot) [default: auto]',
                default=None, type='float')]

script_info[
    'option_label'] = {'mapping_fp': 'QIIME-formatted mapping filepath',
                       'output_dir': 'output directory',
                       'distance_matrix_fp':
                       'distance matrix filepath',
                       'field': 'field in mapping file',
                       'comparison_groups': 'field states to compare',
                       'imagetype': 'output image format',
                       'save_raw_data': 'save raw data used in plot',
                       'plot_type': 'output plot type',
                       'width': 'image width',
                       'height': 'image height',
                       'x_tick_labels_orientation':
                       'x-axis tick label '
                       'orientation',
                       'label_type': 'x-axis label type',
                       'y_min': 'y-axis min',
                       'y_max': 'y-axis max',
                       'transparent': 'make images transparent',
                       'whisker_length': 'whisker length as function '
                       'of IQR',
                       'error_bar_type': 'type of error bars to use ',
                       'distribution_width': 'width of each '
                       'distribution'}

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # Create the output dir if it doesn't already exist.
    try:
        create_dir(opts.output_dir)
    except:
        option_parser.error("Could not create or access output directory "
                            "specified with the -o option.")

    # Parse the distance matrix and mapping file.
    try:
        dist_matrix_header, dist_matrix = parse_distmat(
            open(opts.distance_matrix_fp, 'U'))
    except:
        option_parser.error("This does not look like a valid distance matrix "
                            "file. Please supply a valid distance matrix file using the -d "
                            "option.")
    try:
        mapping, mapping_header, mapping_comments = parse_mapping_file(
            open(opts.mapping_fp, 'U'))
    except QiimeParseError:
        option_parser.error("This does not look like a valid metadata mapping "
                            "file. Please supply a valid mapping file using the -m option.")

    # Make sure the y_min and y_max options make sense, as they can be either
    # 'auto' or a number.
    y_min = opts.y_min
    y_max = opts.y_max
    try:
        y_min = float(y_min)
    except ValueError:
        if y_min == 'auto':
            y_min = None
        else:
            option_parser.error("The --y_min option must be either a number "
                                "or 'auto'.")
    try:
        y_max = float(y_max)
    except ValueError:
        if y_max == 'auto':
            y_max = None
        else:
            option_parser.error("The --y_max option must be either a number "
                                "or 'auto'.")

    # Parse the field states that will be compared to every other field state.
    comparison_field_states = opts.comparison_groups
    comparison_field_states = map(strip, comparison_field_states.split(','))
    comparison_field_states = [field_state.strip('"').strip("'")
                               for field_state in comparison_field_states]
    if comparison_field_states is None:
        option_parser.error("You must provide at least one field state to "
                            "compare (using the -c option).")

    # Get distance comparisons between each field state and each of the
    # comparison field states.
    field = opts.field
    comparison_groupings = get_field_state_comparisons(dist_matrix_header,
                                                       dist_matrix, mapping_header, mapping, field,
                                                       comparison_field_states)

    # Grab a list of all field states that had the comparison field states
    # compared against them. These will be plotted along the x-axis.
    field_states = comparison_groupings.keys()

    def custom_comparator(x, y):
        try:
            num_x = float(x)
            num_y = float(y)
            return int(num_x - num_y)
        except:
            if x < y:
                return -1
            elif x > y:
                return 1
            else:
                return 0

    # Sort the field states as numbers if the elements are numbers, else sort
    # them lexically.
    field_states.sort(custom_comparator)

    # If the label type is numeric, get a list of all field states in sorted
    # numeric order. These will be used to determine the spacing of the
    # field state 'points' along the x-axis.
    x_spacing = None
    if opts.label_type == "numeric":
        try:
            x_spacing = sorted(map(float, field_states))
        except:
            option_parser.error("The 'numeric' label type is invalid because "
                                "not all field states could be converted into "
                                "numbers. Please specify a different label "
                                "type.")

    # Accumulate the data for each field state 'point' along the x-axis.
    plot_data = []
    plot_x_axis_labels = []
    for field_state in field_states:
        field_state_data = []
        for comp_field_state in comparison_field_states:
            field_state_data.append(
                comparison_groupings[field_state][comp_field_state])
        plot_data.append(field_state_data)
        plot_x_axis_labels.append(field_state)

    # Plot the data and labels.
    plot_title = "Distance Comparisons"
    plot_x_label = field
    plot_y_label = "Distance"

    # If we are creating a bar chart or box plot, grab a list of good data
    # colors to use.
    plot_type = opts.plot_type
    plot_colors = None
    if plot_type == "bar" or plot_type == "box":
        plot_colors = [matplotlib_rgb_color(data_colors[color].toRGB())
                       for color in data_color_order]

    assert plot_data, "Error: there is no data to plot!"

    width = opts.width
    height = opts.height
    if width <= 0 or height <= 0:
        option_parser.error("The specified width and height of the image must "
                            "be greater than zero.")

    plot_figure = grouped_distributions(
        opts.plot_type, plot_data, x_values=x_spacing,
        data_point_labels=plot_x_axis_labels,
        distribution_labels=comparison_field_states,
        distribution_markers=plot_colors, x_label=plot_x_label,
        y_label=plot_y_label, title=plot_title,
        x_tick_labels_orientation=opts.x_tick_labels_orientation, y_min=y_min,
        y_max=y_max, whisker_length=opts.whisker_length,
        error_bar_type=opts.error_bar_type,
        distribution_width=opts.distribution_width, figure_width=width,
        figure_height=height)

    # Save the plot in the specified format.
    output_plot_fp = join(opts.output_dir, "%s_Distance_Comparisons.%s" %
                          (field, opts.imagetype))
    plot_figure.savefig(output_plot_fp, format=opts.imagetype,
                        transparent=opts.transparent)

    if not opts.suppress_significance_tests:
        sig_tests_f = open(join(opts.output_dir, "%s_Stats.txt" % field), 'w')

        # Rearrange the plot data into a format suitable for all_pairs_t_test.
        sig_tests_labels = []
        sig_tests_data = []
        for data_point, data_point_label in zip(plot_data, plot_x_axis_labels):
            for dist, comp_field in zip(data_point, comparison_field_states):
                sig_tests_labels.append('%s vs %s' % (data_point_label,
                                                      comp_field))
                sig_tests_data.append(dist)

        sig_tests_results = all_pairs_t_test(sig_tests_labels, sig_tests_data,
                                             tail_type=opts.tail_type,
                                             num_permutations=opts.num_permutations)
        sig_tests_f.write(sig_tests_results)
        sig_tests_f.close()

    if opts.save_raw_data:
        # Write the raw plot data into a tab-delimited file, where each line
        # has the distances between a comparison group and another field state
        # 'point' along the x-axis.
        assert (len(plot_x_axis_labels) == len(plot_data)), "The number of " +\
            "labels do not match the number of points along the x-axis."
        raw_data_fp = join(opts.output_dir,
                           "%s_Distance_Comparisons.txt" % field)
        raw_data_f = open(raw_data_fp, 'w')

        raw_data_f.write("#ComparisonGroup\tFieldState\tDistances\n")
        for label, data in zip(plot_x_axis_labels, plot_data):
            assert (len(comparison_field_states) == len(data)), "The " +\
                "number of specified comparison groups does not match " +\
                "the number of groups found at the current point along " +\
                "the x-axis."
            for comp_field_state, comp_grp_data in zip(comparison_field_states, data):
                raw_data_f.write(comp_field_state + "\t" + label + "\t" +
                                 "\t".join(map(str, comp_grp_data)) + "\n")
        raw_data_f.close()

if __name__ == "__main__":
    main()
