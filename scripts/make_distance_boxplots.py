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

from operator import itemgetter
from os import path
from string import strip
from cogent.util.misc import create_dir
from numpy import median
from qiime.group import get_all_grouped_distances, get_grouped_distances
from qiime.parse import parse_distmat, parse_mapping_file, QiimeParseError
from qiime.pycogent_backports.distribution_plots import generate_box_plots
from qiime.util import get_options_lookup, make_option, \
                       parse_command_line_parameters

script_info = {}
script_info['brief_description'] = "Creates boxplots to compare distances \
                                    between categories"

script_info['script_description'] = """
This script creates boxplots that allow for the comparison between different \
categories found within the mapping file. The boxplots that are created show \
the distances within all samples of a field value, as well as between \
different field values. Individual within and individual between distances \
are also plotted.
"""

script_info['script_usage'] = [("Compare distances between Fast and Control "
                                "for Treatment field",
                                "This example will generate an image "
                                "with boxplots for all within and all between "
                                "distances for the field Treatment, and will "
                                "also include plots for individual within "
                                "(e.g. Control vs. Control, Fast vs. Fast) "
                                "and individual between (e.g. Control vs. "
                                "Fast). The generated plot image PDF will be "
                                "written to the output directory 'out_files'.",
                                "%prog -d dist_matrix.txt -m map.txt -f "
                                "\"Treatment\" -o out_files"),
                               ("Only plot individual distances",
                                "This example will generate a PNG of all of "
                                "individual distances (within and between) "
                                "for the Treatment field.",
                                "%prog -d dist_matrix.txt -m map.txt -f "
                                "\"Treatment\" -o out_files -g png "
                                "--suppress_all_within "
                                "--suppress_all_between"),
                               ("Save raw data",
                                "This example will generate an SVG image of "
                                "the boxplots and also output the plotting "
                                "data to a tab-delimited file.",
                                "%prog -d dist_matrix.txt -m map.txt -f "
                                "\"Treatment\" -o out_files -g svg "
                                "--save_raw_data")]

script_info['output_description'] = "Images of the plots are written to \
                                    the specified output directory (one image \
                                    per field). The raw data used in the \
                                    plots can optionally be written into \
                                    tab-delimited files (one file per field)."

# Get a dictionary of common QIIME options.
options = get_options_lookup()

script_info['required_options'] = [
    options['mapping_fp'],
    options['output_dir'],
    make_option('-d', '--distance_matrix_fp',
        help='input distance matrix filepath (i.e. the result of '
             'beta_diversity.py)',
        type='existing_filepath'),
    make_option('-f', '--fields',
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
    make_option('--y_min', type='string',
        help='the minimum y-axis value in the resulting plot. If "auto", '
             'it is automatically calculated [default: %default]',
        default=0),
    make_option('--y_max', type='string',
        help='the maximum y-axis value in the resulting plot. If "auto", '
             'it is automatically calculated [default: %default]',
        default=1),
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
    make_option('--sort', action='store_true',
        help='sort boxplots by increasing median. If no sorting is applied, '
             'boxplots will be grouped logically as follows: all within, all '
             'between, individual within, and individual between '
             '[default: %default]', default=False)]

script_info['option_label'] = {'mapping_fp':'QIIME-formatted mapping filepath',
                               'output_dir':'output directory',
                               'distance_matrix_fp':'distance matrix filepath',
                               'fields':'categories to compare',
                               'imagetype':'output image format',
                               'save_raw_data':'save raw data used in plots',
                               'suppress_all_within':'suppress all within plot',
                               'suppress_all_between':'suppress all between '
                                   'plot',
                               'suppress_individual_within':'suppress '
                                   'individual within plot(s)',
                               'suppress_individual_between':'suppress '
                                   'individual between plot(s)',
                               'y_min':'y-axis min',
                               'y_max':'y-axis max',
                               'width':'image width',
                               'height':'image height',
                               'transparent':'make images transparent',
                               'whisker_length':'whisker length as function '
                                   'of IQR',
                               'box_width':'width of boxes',
                               'sort':'sort boxplots by ascending median'}

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

    fields = opts.fields
    fields = map(strip, fields.split(','))
    fields = [field.strip('"').strip("'") for field in fields]

    if fields is None:
        option_parser.error("You must provide at least one field using the -f "
                            "option.")

    # Make sure each field is in the mapping file.
    for field in fields:
        if field not in mapping_header:
            option_parser.error("The field '%s' is not in the provided "
                "mapping file. Please supply correct fields (using the -f "
                "option) corresponding to fields in the mapping file."
                % field)

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

    # Generate the various boxplots, depending on what the user wanted
    # suppressed. Add them all to one encompassing plot.
    for field in fields:
        plot_data = []
        plot_labels = []

        if not opts.suppress_all_within:
            plot_data.append(get_all_grouped_distances(dist_matrix_header,
                    dist_matrix, mapping_header, mapping, field, within=True))
            plot_labels.append("All within %s" % field)
        if not opts.suppress_all_between:
            plot_data.append(get_all_grouped_distances(dist_matrix_header,
                    dist_matrix, mapping_header, mapping, field, within=False))
            plot_labels.append("All between %s" % field)
        if not opts.suppress_individual_within:
            within_dists = get_grouped_distances(dist_matrix_header,
                    dist_matrix, mapping_header, mapping, field, within=True)
            for grouping in within_dists:
                plot_data.append(grouping[2])
                plot_labels.append("%s vs. %s" % (grouping[0], grouping[1]))
        if not opts.suppress_individual_between:
            between_dists = get_grouped_distances(dist_matrix_header,
                    dist_matrix, mapping_header, mapping, field, within=False)
            for grouping in between_dists:
                plot_data.append(grouping[2])
                plot_labels.append("%s vs. %s" % (grouping[0], grouping[1]))

        # We now have our data and labels ready, so plot them!
        assert (len(plot_data) == len(plot_labels)), "The number " +\
                "of boxplot labels does not match the number of " +\
                "boxplots."
        if plot_data:
            if opts.sort:
                # Sort our plot data in order of increasing median.
                sorted_data = []
                for label, distribution in zip(plot_labels, plot_data):
                    sorted_data.append((label, distribution,
                        median(distribution)))
                sorted_data.sort(key=itemgetter(2))
                plot_labels = []
                plot_data = []
                for label, distribution, median_value in sorted_data:
                    plot_labels.append(label)
                    plot_data.append(distribution)

            plot_figure = generate_box_plots(plot_data,
                    x_tick_labels=plot_labels, title="%s Distances" % field,
                    x_label="Grouping", y_label="Distance",
                    x_tick_labels_orientation='vertical', y_min=y_min,
                    y_max=y_max, whisker_length=opts.whisker_length,
                    box_width=opts.box_width)
            width = opts.width
            height = opts.height
            if width is None:
                width = len(plot_data) * opts.box_width + 2
            if width > 0 and height > 0:
                plot_figure.set_size_inches(width, height)
            else:
                option_parser.error("The specified width and height of the "
                                    "image must be greater than zero.")
            output_plot_fp = path.join(opts.output_dir, "%s_Distances.%s"
                                       % (field, opts.imagetype))
            plot_figure.savefig(output_plot_fp, format=opts.imagetype,
                    transparent=opts.transparent)
        else:
            option_parser.error("You have chosen to suppress all plots. At "
                                "least one type of plot must be unsuppressed.")

        if opts.save_raw_data:
            # Write the raw plot data into a tab-delimited file.
            assert(len(plot_labels) == len(plot_data))
            raw_data_fp = path.join(opts.output_dir, "%s_Distances.xls"
                                    % field)
            raw_data_f = open(raw_data_fp, 'w')

            for label, data in zip(plot_labels, plot_data):
                raw_data_f.write(label.replace(" ", "_") + "\t")
                raw_data_f.write("\t".join(map(str, data)))
                raw_data_f.write("\n")
            raw_data_f.close()


if __name__ == "__main__":
    main()
