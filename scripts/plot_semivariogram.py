#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Antonio Gonzalez Pena", "Kyle Patnode", "Yoshiki Vazquez-Baeza"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"

from qiime.plot_taxa_summary import make_legend
from qiime.colors import get_qiime_hex_string_color
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from qiime.plot_semivariogram import fit_semivariogram, FitModel
from qiime.parse import parse_distmat, parse_mapping_file
from qiime.filter import (filter_samples_from_distance_matrix,
                          sample_ids_from_metadata_description)
from matplotlib import use
use('Agg', warn=False)
from matplotlib.pyplot import (plot, xlabel, ylabel, title, savefig, ylim,
                               xlim, legend, show, figure)
from numpy import asarray
import os
from os.path import splitext
from StringIO import StringIO
from copy import deepcopy

options_lookup = get_options_lookup()

script_info = {}

script_info['brief_description'] = "Fits a model between two distance matrices " +\
    "and plots the result"
script_info['script_description'] = "Fits a spatial autocorrelation model " +\
    "between two matrices and plots the result. This script will work with " +\
    "two distance matrices but will ignore the 0s at the diagonal and the " +\
    "values that go to N/A. See distance_matrix_from_mapping.py."
script_info['script_usage'] = []
script_info['script_usage'].append(("Fitting", "For this script, the user "
                                    "supplies two distance matrices (i.e. resulting file from "
                                    "beta_diversity.py), along with the output filename (e.g. semivariogram), "
                                    "and the model to fit, as follows:", "%prog -x distance.txt -y unifrac.txt "
                                    "-o semivariogram_exponential.png"))
script_info['script_usage'].append(("", "Modify the the default method to "
                                    "gaussian", "%prog -x distance.txt -y unifrac.txt --model gaussian -o "
                                    "semivariogram_gaussian.png"))
script_info['script_usage'].append(("Color semivariograms by a category in"
                                    " the metadata mapping file", "Using a header name in the mapping file"
                                    " (Time), create two separate semivariograms in the same plot, an "
                                    "accompanying file with the color coding will be created"
                                    "(categories_legend.eps), both the legends and the plot will be in eps "
                                    "format.", "%prog -y unweighted_unifrac_dm.txt -x time_dm.txt --model "
                                    "gaussian -m Fasting_Map.txt -o categories.eps -c Treatment"))
script_info['output_description'] = "The resulting output file consists of a " +\
    "pdf image containing the plot between the two distances matrices and the" +\
    " fitted model"

script_info['required_options'] = [
    make_option('-x', '--input_path_x', type='existing_filepath',
                help='path to distance matrix to be displayed in the x axis'),
    make_option('-y', '--input_path_y', type='existing_filepath',
                help='path to distance matrix to be displayed in the y axis'),
    make_option('-o', '--output_path', type='new_path',
                help='output path. directory for batch processing, ' +
                'filename for single file operation'),
]
script_info['optional_options'] = [
    make_option('-b', '--binning', type='string',
                default=None, help='binning ranges. Format: [increment,top_limit], when ' +
                'top_limit is -1=infinitum; you can specify several ranges using the same ' +
                'format, i.e. [2.5,10][50,-1] will set two bins, one from 0-10 using 2.5 ' +
                'size steps and from 10-inf using 50 size steps. Note that the binning is ' +
                'used to clean the plots (reduce number of points) but ignored to fit the ' +
                'model. [default: %default]'),
    make_option('--ignore_missing_samples', help='This will overpass the error raised ' +
                'when the matrices have different sizes/samples', action='store_true', default=False),
    make_option(
        '--x_max',
        type='float',
        help='x axis max limit [default: auto]',
        default=None),

    make_option(
        '--x_min',
        type='float',
        help='x axis min limit [default: auto]',
        default=None),

    make_option(
        '--y_max',
        type='float',
        help='y axis max limit [default: auto]',
        default=None),

    make_option(
        '--y_min',
        type='float',
        help='y axis min limit [default: auto]',
        default=None),
    make_option(
        '-X', '--x_label', default='Distance Dissimilarity (m)', type='string',
        help='Label for the x axis [default: %default]'),
    make_option(
        '-Y', '--y_label', default='Community Dissimilarity', type='string',
        help='Label for the y axis [default: %default]'),
    make_option('-t', '--fig_title', default='Semivariogram', type='string',
                help='Title of the plot [default: %default]'),
    make_option('--dot_color', type='string', help='dot color for plot, more info:' +
                ' http://matplotlib.sourceforge.net/api/pyplot_api.html' +
                ' [default: %default]', default="white"),
    make_option('--dot_marker', type='string', help='dot color for plot, more info:' +
                ' http://matplotlib.sourceforge.net/api/pyplot_api.html' +
                ' [default: %default]', default="o"),
    make_option('--line_color', type='string', help='line color for plot, more info:' +
                ' http://matplotlib.sourceforge.net/api/pyplot_api.html' +
                ' [default: %default]', default="blue"),
    make_option('--dot_alpha', type='float', help='alpha for dots, more info:' +
                ' http://matplotlib.sourceforge.net/api/pyplot_api.html' +
                ' [default: %default]', default=1),
    make_option('--line_alpha', type='float', help='alpha for dots, more info:' +
                ' http://matplotlib.sourceforge.net/api/pyplot_api.html' +
                ' [default: %default]', default=1),
    make_option('--model', type='choice',
                choices=FitModel.options, default='exponential',
                help='model to be fitted to the data. Valid ' +
                'choices are:' + ', '.join(FitModel.options) + '. [default: %default]'),
    make_option('-p', '--print_model', action='store_true',
                help='Print in the title of the plot the function of the fit. ' +
                '[default: %default]', default=False),
    make_option('-c', '--category', type='string', help='category to color each of'
                ' the trajectories when you have multiple treatments [default: %default]',
                default=None),
    make_option('-m', '--mapping_fp', type='existing_filepath', help='metadata '
                'mapping file, only used when coloring by a category, a file with the '
                'legends and color coding will be created with the suffix legend '
                '[default: %default]',
                default=None)
]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    category = opts.category
    mapping_fp = opts.mapping_fp

    colors_used = []

    if (category and mapping_fp is None) or (category is None and mapping_fp):
        option_parser.error('If coloring by a metadata category, both the '
                            'category and the mapping file must be supplied.')
    elif mapping_fp and category:
        mapping_data, mapping_headers, _ = parse_mapping_file(open(mapping_fp,
                                                                   'U'))
        if category not in mapping_headers:
            option_parser.error("The category supplied must exist in the "
                                "metadata mapping file, '%s' does not exist." % category)
        index = mapping_headers.index(category)
        categories = list(set([line[index] for line in mapping_data]))
    list_of_plots = []

    if opts.binning is None:
        ranges = []
    else:
        # simple ranges format validation
        if opts.binning.count('[') != opts.binning.count(']') or\
                opts.binning.count('[') != opts.binning.count(','):
            raise ValueError("The binning input has an error: '%s'; " % +
                             "\nthe format should be [increment1,top_limit1][increment2,top_limit2]")
        # spliting in ranges
        rgn_txt = opts.binning.split('][')
        # removing left [ and right ]
        rgn_txt[0] = rgn_txt[0][1:]
        rgn_txt[-1] = rgn_txt[-1][:-1]
        # converting into int
        ranges = []
        max = 0

        for i, r in enumerate(rgn_txt):
            try:
                values = map(float, r.split(','))
            except ValueError:
                raise ValueError(
                    "Not a valid format for binning %s" %
                    opts.binning)
            if len(values) != 2:
                raise ValueError(
                    "All ranges must have only 2 values: [%s]" %
                    r)
            elif i + 1 != len(rgn_txt):
                if values[0] > values[1]:
                    raise ValueError(
                        "The bin value can't be greater than the max value: [%s]" %
                        r)
                elif values < 0:
                    raise ValueError(
                        "This value can not be negative: [%s]" %
                        r)
                elif max > values[1]:
                    raise ValueError(
                        "This value can not smaller than the previous one: [%s]" %
                        r)
                else:
                    max = values[1]

            ranges.append(values)

    x_samples, x_distmtx = parse_distmat(open(opts.input_path_x, 'U'))
    y_samples, y_distmtx = parse_distmat(open(opts.input_path_y, 'U'))

    if opts.ignore_missing_samples:
        ignoring_from_x = list(set(x_samples) - set(y_samples))
        ignoring_from_y = list(set(y_samples) - set(x_samples))

        if opts.verbose:
            print '\nFrom %s we are ignoring: %s\n' % (opts.input_path_x, ignoring_from_x)
            print '\nFrom %s we are ignoring: %s\n' % (opts.input_path_y, ignoring_from_y)
            print '\nOnly using: %s\n' % (list(set(x_samples) & set(y_samples)))

        x_file = StringIO(
            filter_samples_from_distance_matrix((x_samples, x_distmtx), ignoring_from_x))
        x_samples, x_distmtx = parse_distmat(x_file)

        y_file = StringIO(
            filter_samples_from_distance_matrix((y_samples, y_distmtx), ignoring_from_y))
        y_samples, y_distmtx = parse_distmat(y_file)
    else:
        if x_distmtx.shape != y_distmtx.shape:
            raise ValueError('The distance matrices have different sizes. ' +
                             'You can cancel this error by passing --ignore_missing_samples')

    figure()
    if category is None:
        x_val, y_val, x_fit, y_fit, func_text = fit_semivariogram(
            (x_samples, x_distmtx), (y_samples, y_distmtx), opts.model, ranges)

        plot(
            x_val,
            y_val,
            color=opts.dot_color,
            marker=opts.dot_marker,
            linestyle="None",
            alpha=opts.dot_alpha)
        plot(
            x_fit,
            y_fit,
            linewidth=2.0,
            color=opts.line_color,
            alpha=opts.line_alpha)
    else:
        # not all the categories that are going to be enumerated are found in
        # the distance matrices i.e. the mapping file is a superset that can
        # contain more samples than the distance matrices
        used_categories = deepcopy(categories)

        for index, single_category in enumerate(categories):
            good_sample_ids = sample_ids_from_metadata_description(
                open(mapping_fp), '%s:%s' % (category, single_category))

            try:
                _y_samples, _y_distmtx = parse_distmat(StringIO(
                    filter_samples_from_distance_matrix((y_samples, y_distmtx),
                                                        good_sample_ids, negate=True)))
                _x_samples, _x_distmtx = parse_distmat(StringIO(
                    filter_samples_from_distance_matrix((x_samples, x_distmtx),
                                                        good_sample_ids, negate=True)))
            except ValueError:
                # no samples found for this category
                used_categories.remove(single_category)
                continue

            x_val, y_val, x_fit, y_fit, func_text = fit_semivariogram(
                (_x_samples, _x_distmtx), (_y_samples, _y_distmtx),
                opts.model, ranges)

            # retrieve one of the colors the "QIIME" colors and add it to the
            # list of used colors for the creation of the legends in the plot
            color_only = get_qiime_hex_string_color(index)
            colors_used.append(color_only)

            plot(x_val, y_val, color=color_only, marker=opts.dot_marker,
                 linestyle="None", alpha=opts.dot_alpha)
            plot(x_fit, y_fit, linewidth=2.0, color=color_only,
                 alpha=opts.line_alpha, label=single_category)

    # set plot limits if requested
    x_lb, x_ub = xlim()
    y_lb, y_ub = ylim()
    if opts.x_min is not None:
        x_lb = opts.x_min
    if opts.x_max is not None:
        x_ub = opts.x_max
    if opts.y_min is not None:
        y_lb = opts.y_min
    if opts.y_max is not None:
        y_ub = opts.y_max
    xlim(x_lb, x_ub)
    ylim(y_lb, y_ub)


    x_label = opts.x_label
    y_label = opts.y_label
    fig_title = '%s (%s)' % (opts.fig_title, opts.model)

    xlabel(x_label)
    ylabel(y_label)
    if opts.print_model:
        title(fig_title + ' ' + func_text)
    else:
        title(fig_title)

    savefig(opts.output_path)

    # print the legends after the figure is exported to avoid conflicts
    if category:
        # if there's a desired format, use that, else default it to png
        _, extension = splitext(opts.output_path)

        # remove the dot, else, make_legend will add it to the filename
        extension = extension.replace('.', '')

        if extension == '':
            extension = 'png'
        make_legend(used_categories, colors_used, 0, 0, 'black', 'white',
                    opts.output_path, extension, 80)

if __name__ == "__main__":
    main()
