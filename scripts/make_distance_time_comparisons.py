#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Rideout"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Jai Rideout"]
__license__ = "GPL"
__version__ = "1.3.0-dev"
__maintainer__ = "Jai Rideout"
__email__ = "jr378@nau.edu"
__status__ = "Development"

from os import path
from string import strip
from cogent.util.misc import create_dir
from qiime.colors import data_colors, data_color_order
from qiime.group import get_grouped_distances
from qiime.make_distance_histograms import matplotlib_rgb_color
from qiime.parse import group_by_field, parse_distmat, parse_mapping_file, \
                        QiimeParseError
from qiime.pycogent_backports.distribution_plots import \
        generate_comparative_plots
from qiime.util import get_options_lookup, make_option, \
                       parse_command_line_parameters

script_info = {}
script_info['brief_description'] = "Creates plots to compare distances \
                                    between timepoints"

script_info['script_description'] = """
This script creates plots (bar charts or scatter plots) that allow for the \
comparison between field values at various timepoints found within the \
mapping file.

A special time field is required to be in the mapping file. This field should \
contain timepoint values. For example, a time field might contain the values \
1, 2, 3, 4, and 5, which label samples that are from day 1, day 2, day 3, and \
so on. This time field can be specified when the script is run, as well as \
the groups to compare at each timepoint. For example, two comparison groups \
might be 1 and 2. The resulting plot would contain timepoints for days 3, 4, \
and 5, and at each timepoint, the distances between day 1 and that timepoint \
would be plotted, as well as the distances between day 2 and the timepoint.
"""

script_info['script_usage'] = [("Compare distances between Native and Input "
                                "samples for each timepoint in the Time field",
                                "This example will generate a bar chart "
                                "with the distances between Native samples "
                                "and each timepoint, as well as the distances "
                                "between Input samples and each timepoint. "
                                "The output image will be put in the "
                                "'out_files' directory.",
                                "%prog -d dist_matrix.txt -m map.txt -f "
                                "Time -c \"Native,Input\" -o out_files")]

script_info['output_description'] = "An image of the plot is written to \
                                    the specified output directory."

# Get a dictionary of common QIIME options.
options = get_options_lookup()

script_info['required_options'] = [
    options['mapping_fp'],
    options['output_dir'],
    make_option('-d', '--distance_matrix_fp',
        help='input distance matrix filepath (i.e. the result of '
             'beta_diversity.py)',
        type='existing_filepath'),
    make_option('-f', '--time_field',
        help='time field in the mapping file to compare distances on'),
    make_option('-c', '--comparison_groups',
        help='comma-separated list of time field values to compare at each '
             'timepoint, where the list of field values should be in quotes '
             '(e.g. "FieldVal1,FieldVal2,FieldVal3")')]

script_info['optional_options'] = [
    make_option('-g', '--imagetype',
        help='type of image to produce (i.e. png, svg, pdf) '
             '[default: %default]', default='pdf', type="choice",
        choices=['pdf','png','svg']),
    make_option('--save_raw_data', action='store_true',
        help='store raw data used to create plot in a tab-delimited file '
             '[default: %default]',
        default=False),
    make_option('-t', '--plot_type',
        help='type of plot to produce (bar chart, scatter plot, or box plot) '
             '[default: %default]',
        default='bar', type='choice', choices=['bar', 'scatter', 'box']),
    make_option('--width',
        help='width of the output image in inches [default: %default]',
        default=12, type='float'),
    make_option('--height',
        help='height of the output image in inches [default: %default]',
        default=6, type='float'),
    make_option('--x_tick_labels_orientation',
        help='type of orientation for x-axis tick labels [default: %default]',
        default='vertical', type='choice', choices=['vertical', 'horizontal']),
    make_option('-a', '--label_type', type='choice',
        help='Label type ("numeric" or "categorical"). '
        'If the label type is defined as numeric, the x-axis will be '
        'scaled accordingly. Otherwise the x-values will treated '
        'categorically and be evenly spaced [default: %default].',
        choices=['categorical','numeric'], default='categorical'),
    make_option('--y_min', type='string',
        help='the minimum y-axis value in the resulting plot. If "auto", '
             'it is automatically calculated [default: %default]',
        default=0),
    make_option('--y_max', type='string',
        help='the maximum y-axis value in the resulting plot. If "auto", '
             'it is automatically calculated [default: %default]',
        default=1),
    make_option('--transparent', action='store_true',
        help='make output images transparent (useful for overlaying an image '
             'on top of a colored background ) [default: %default]',
        default=False),
    make_option('--whisker_length',
        help='if plot_type is "box", determines the length of the whiskers as '
             'a function of the IQR. For example, if 1.5, the whiskers extend '
             'to 1.5 * IQR. Anything outside of that range is seen as an '
             'outlier. If plot_type is not "box", this option is ignored '
             '[default: %default]',
        default='1.5', type='float'),
    make_option('--distribution_width',
        help='width (in plot units) of each individual distribution (e.g. each '
             'bar if the plot type is a bar chart, or the width of each box '
             'if the plot type is a boxplot) [default: %default]',
        default='0.4', type='float'),
    make_option('--group_spacing',
        help='width (in plot units) of the gap between each timepoint (i.e. '
             'the width between each group of distributions) '
             '[default: %default]',
        default='0.5', type='float')]

script_info['option_label'] = {'mapping_fp':'QIIME-formatted mapping filepath',
                               'output_dir':'output directory',
                               'distance_matrix_fp':'distance matrix filepath',
                               'time_field':'time field in mapping file',
                               'comparison_groups':'time field values to '
                                   'compare',
                               'imagetype':'output image format',
                               'save_raw_data':'save raw data used in plot',
                               'plot_type':'output plot type',
                               'width':'image width',
                               'height':'image height',
                               'x_tick_labels_orientation':'x-axis tick label '
                                   'orientation',
                               'label_type':'x-axis label type',
                               'y_min':'y-axis min',
                               'y_max':'y-axis max',
                               'transparent':'make images transparent',
                               'whisker_length':'whisker length as function '
                                   'of IQR',
                               'distribution_width':'width of each '
                                   'distribution',
                               'group_spacing':'width of gap between '
                                   'timepoints'}

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

    # Get the time field and make sure it is in the mapping file.
    time_field = opts.time_field
    if time_field not in mapping_header:
        option_parser.error("The field '%s' is not in the provided "
            "mapping file. Please supply a time field (using the -f "
            "option) corresponding to the time field in the mapping file."
            % time_field)

    # Parse the time field values that will be compared at each timepoint.
    comp_groups = opts.comparison_groups
    comp_groups = map(strip, comp_groups.split(','))
    comp_groups = [group.strip('"').strip("'") for group in comp_groups]

    if comp_groups is None:
        option_parser.error("You must provide at least one group to compare "
                            "using the -c option.")

    # Make sure each comparison group is in the mapping file's time field.
    mapping_data = [mapping_header]
    mapping_data.extend(mapping)
    groups = group_by_field(mapping_data, time_field)
    for group in comp_groups:
        if group not in groups:
            option_parser.error("The comparison group '%s' is not in the "
                                "provided mapping file's time field '%s'. "
                                "Please supply correct comparison groups "
                                "(using the -c option) corresponding to valid "
                                "values in the time field."
                                % (group, time_field))

    # Grab a list of the timepoints that will be plotted. This list of
    # timepoints should not include the values of the time field that are being
    # compared at each timepoint.
    timepoints = [group for group in groups.keys() if group not in comp_groups]

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

    # Sort the timepoints as numbers if the elements are numbers, else sort
    # them lexically.
    timepoints.sort(custom_comparator)

    # If the label type is numeric, get a list of all timepoints in sorted
    # numeric order. These will be used to determine the spacing of the
    # timepoints along the x-axis.
    x_spacing = None
    if opts.label_type == "numeric":
        try:
            x_spacing = map(float, timepoints)
            x_spacing.sort()
        except:
            option_parser.error("The 'numeric' label type is invalid because "
                                "not all timepoints could be converted into "
                                "numbers. Please specify a different label "
                                "type.")

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

    # Get all between distance groupings.
    groupings = get_grouped_distances(dist_matrix_header, dist_matrix,
            mapping_header, mapping, time_field, within=False)

    # Accumulate the data for each timepoint.
    plot_data = []
    plot_timepoint_labels = []
    for timepoint in timepoints:
        timepoint_data = []
        for comp_group in comp_groups:
            matching_data = []
            for group in groupings:
                if ((group[0] == timepoint or group[1] == timepoint)
                    and (group[0] == comp_group or group[1] == comp_group)):
                    matching_data = group[2]
            timepoint_data.append(matching_data)
        plot_data.append(timepoint_data)
        plot_timepoint_labels.append(timepoint)

    # Plot the data and labels.
    plot_title = "Timepoint Distances"
    plot_x_label = time_field
    plot_y_label = "Distance"

    # If we are creating a bar chart or box plot, grab a list of good data
    # colors to use.
    plot_type = opts.plot_type
    plot_colors = None
    if plot_type == "bar" or plot_type == "box":
        plot_colors = [matplotlib_rgb_color(data_colors[color].toRGB()) \
                       for color in data_color_order]

    assert(plot_data)
    plot_figure = generate_comparative_plots(opts.plot_type, plot_data,
            x_values=x_spacing, data_point_labels=plot_timepoint_labels,
            distribution_labels=comp_groups, distribution_markers=plot_colors,
            x_label=plot_x_label, y_label=plot_y_label, title=plot_title,
            x_tick_labels_orientation=opts.x_tick_labels_orientation,
            y_min=y_min, y_max=y_max, whisker_length=opts.whisker_length,
            distribution_width=opts.distribution_width,
            group_spacing=opts.group_spacing)

    # Save the plot in the specified format and size.
    width = opts.width
    height = opts.height
    if width > 0 and height > 0:
        plot_figure.set_size_inches(width, height)
    else:
        option_parser.error("The specified width and height of the image must "
                            "be greater than zero.")
    output_plot_fp = path.join(opts.output_dir, "%s_Timepoint_Distances.%s"
                               % (time_field, opts.imagetype))
    plot_figure.savefig(output_plot_fp, format=opts.imagetype,
            transparent=opts.transparent)

    if opts.save_raw_data:
        # Write the raw plot data into a tab-delimited file, where each line is
        # a comparison group's data at a specific timepoint.
        assert(len(plot_timepoint_labels) == len(plot_data))
        raw_data_fp = path.join(opts.output_dir, "%s_Timepoint_Distances.xls"
                                % time_field)
        raw_data_f = open(raw_data_fp, 'w')

        raw_data_f.write("#ComparisonGroup\tTimePoint\tDistances\n")
        for label, data in zip(plot_timepoint_labels, plot_data):
            assert(len(comp_groups) == len(data))
            for comp_grp, comp_grp_data in zip(comp_groups, data):
                raw_data_f.write(comp_grp + "\t" + label + "\t" +
                        "\t".join(map(str, comp_grp_data)) + "\n")
        raw_data_f.close()


if __name__ == "__main__":
    main()
