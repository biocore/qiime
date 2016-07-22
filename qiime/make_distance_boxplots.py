#!/usr/bin/env python
from __future__ import division

__author__ = "Jai Ram Rideout"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

"""Contains functions used in the make_distance_boxplots.py script."""

from collections import defaultdict

import numpy as np

from qiime.colors import data_colors, data_color_order
from qiime.format import format_mapping_file
from qiime.group import get_all_grouped_distances, get_grouped_distances
from qiime.colors import matplotlib_rgb_color
from qiime.parse import parse_distmat, parse_mapping_file
from skbio.draw import boxplots
from qiime.util import MetadataMap

SORT_TYPES = ('median', 'alphabetical')

def make_distance_boxplots(dm_f,
                           map_f,
                           fields,
                           width=None,
                           height=6.0,
                           suppress_all_within=False,
                           suppress_all_between=False,
                           suppress_individual_within=False,
                           suppress_individual_between=False,
                           y_min=0.0,
                           y_max=1.0,
                           whisker_length=1.5,
                           box_width=0.5,
                           box_color=None,
                           color_individual_within_by_field=None,
                           sort=None):
    """Generates various types of boxplots for distance comparisons.

    Returns a list of tuples, one for each field. Each tuple contains the
    following:
        1) the name of the field (string)
        2) a matplotlib.figure.Figure object containing the boxplots
        3) a list of lists containing the raw plot data that was passed to mpl
        4) a list of labels for each of the boxplots (string)
        5) a list of mpl-compatible colors (one for each boxplot)

    The Figure can be saved, and the raw data and labels can be useful (for
    example) performing statistical tests or writing the raw data to disk.

    The input arguments are exactly derived from the make_distance_boxplots.py
    script (see the script options for details). To avoid duplicated effort,
    their descriptions are not reproduced here.
    """
    # Parse data files and do some preliminary error checking.
    dm_header, dm_data = parse_distmat(dm_f)
    map_data, map_header, map_comments = parse_mapping_file(map_f)

    if fields is None or len(fields) < 1:
        raise ValueError("You must provide at least one field to analyze.")

    for field in fields:
        if field not in map_header:
            raise ValueError("The field '%s' is not in the provided mapping "
                             "file. Please supply correct fields "
                             "corresponding to fields in the mapping file." %
                             field)

    # Make sure the y_min and y_max options make sense, as they can be either
    # 'auto' or a number.
    y_min = _cast_y_axis_extrema(y_min)
    y_max = _cast_y_axis_extrema(y_max)

    # Collate the distributions of distances that will comprise each boxplot.
    # Suppress the generation of the indicated types of boxplots.
    results = []
    for field in fields:
        plot_data = []
        plot_labels = []
        plot_colors = []
        legend = None

        # Little bit of duplicate code here... not sure it's worth the effort
        # to clean up though.
        if not suppress_all_within:
            plot_data.append(get_all_grouped_distances(dm_header, dm_data,
                             map_header, map_data, field, within=True))
            plot_labels.append("All within %s" % field)

            if color_individual_within_by_field is not None:
                plot_colors.append(None)
            else:
                plot_colors.append(box_color)

        if not suppress_all_between:
            plot_data.append(get_all_grouped_distances(dm_header, dm_data,
                             map_header, map_data, field, within=False))
            plot_labels.append("All between %s" % field)

            if color_individual_within_by_field is not None:
                plot_colors.append(None)
            else:
                plot_colors.append(box_color)

        if not suppress_individual_within:
            within_dists = get_grouped_distances(dm_header, dm_data,
                                                 map_header, map_data, field,
                                                 within=True)
            field_states = []
            for grouping in within_dists:
                plot_data.append(grouping[2])
                plot_labels.append("%s vs. %s" % (grouping[0], grouping[1]))
                field_states.append(grouping[0])

            # If we need to color these boxplots by a field, build up a
            # list of colors and a legend.
            if color_individual_within_by_field is not None:
                colors, color_mapping = _color_field_states(
                    format_mapping_file(map_header, map_data).split('\n'),
                    dm_header, field, field_states,
                    color_individual_within_by_field)
                plot_colors.extend(colors)
                legend = (color_mapping.values(), color_mapping.keys())
            else:
                plot_colors.extend([box_color] * len(field_states))

        if not suppress_individual_between:
            between_dists = get_grouped_distances(dm_header, dm_data,
                                                  map_header, map_data, field, within=False)

            for grouping in between_dists:
                plot_data.append(grouping[2])
                plot_labels.append("%s vs. %s" % (grouping[0], grouping[1]))

                if color_individual_within_by_field is not None:
                    plot_colors.append(None)
                else:
                    plot_colors.append(box_color)

        assert (len(plot_data) == len(plot_labels) and
                len(plot_labels) == len(plot_colors)), "The number " +\
            "of boxplot labels and colors do not match the number of " +\
            "boxplots."

        # We now have our data and labels ready, so plot them!
        if plot_data:
            if sort is not None:
                plot_data, plot_labels, plot_colors = _sort_distributions(
                    plot_data, plot_labels, plot_colors, sort)

            if width is None:
                width = len(plot_data) * box_width + 2
            if width <= 0 or height <= 0:
                raise ValueError("The specified width and height of the plot "
                                 "must be greater than zero.")

            plot_figure = boxplots(plot_data, x_tick_labels=plot_labels,
                                   title="%s Distances" % field,
                                   x_label="Grouping", y_label="Distance",
                                   x_tick_labels_orientation='vertical',
                                   y_min=y_min, y_max=y_max,
                                   whisker_length=whisker_length,
                                   box_width=box_width, box_colors=plot_colors,
                                   figure_width=width, figure_height=height,
                                   legend=legend)

            results.append((field, plot_figure, plot_data, plot_labels,
                            plot_colors))
        else:
            raise ValueError("The generation of all plots was suppressed. At "
                             "least one type of plot must be unsuppressed.")

    return results


def _cast_y_axis_extrema(y_axis_extrema):
    """Casts to float; handles 'auto' value as well."""
    try:
        y_axis_extrema = float(y_axis_extrema)
    except ValueError:
        if y_axis_extrema == 'auto':
            y_axis_extrema = None
        else:
            raise ValueError("The min and max y-axis values must be numbers "
                             "or 'auto'. Couldn't handle the value %r." %
                             y_axis_extrema)
    return y_axis_extrema


def _sort_distributions(plot_data, plot_labels, plot_colors, sort_type):
    """Sort plot data, labels, and colors.

    Can sort in order of increasing median or alphabetically (based on plot
    labels).

    """
    if sort_type == 'median':
        sort_key = lambda item: np.median(item[0])
    elif sort_type == 'alphabetical':
        sort_key = lambda item: item[1]
    else:
        raise ValueError("Invalid sort type '%s'." % sort_type)

    # Taken from http://stackoverflow.com/a/9764364
    return zip(*sorted(zip(plot_data, plot_labels, plot_colors), key=sort_key))

def _color_field_states(map_f, samp_ids, field, field_states, color_by_field):
    """Colors one field by another.

    Returns a list of matplotlib-compatible colors, one for each of the input
    field_states. Also returns a dictionary mapping color_by_field states to
    colors (useful for building a legend, for example).

    If there are not enough colors available (they are drawn from
    qiime.colors.data_colors), an error will be raised as the color mapping
    (and legend) will be ambiguous.

    A one-to-one mapping must exist between each field_state and its
    corresponding color_by field state (otherwise it is unclear which
    corresponding color_by field state should be used to color it by). An error
    will be raised if this one-to-one mapping does not exist.

    Arguments:
        map_f - the mapping file (file-like object)
        samp_ids - a list of sample IDs to consider in the mapping file. Only
            these sample IDs will be used when coloring field states
        field - the field in the mapping file to color
        field_states - the field states in field to color
        color_by_field - the field in the mapping file to color field_states by
    """
    colors = []
    color_pool = [matplotlib_rgb_color(data_colors[color].toRGB())
                  for color in data_color_order]
    metadata_map = MetadataMap.parseMetadataMap(map_f)

    for field_to_check in field, color_by_field:
        if field_to_check not in metadata_map.CategoryNames:
            raise ValueError("The field '%s' is not in the metadata mapping "
                             "file's column headers." % field_to_check)

    all_field_states = metadata_map.getCategoryValues(samp_ids, field)
    all_color_by_states = metadata_map.getCategoryValues(samp_ids,
                                                         color_by_field)

    if len(set(field_states) - set(all_field_states)) != 0:
        raise ValueError("Encountered unrecognizable field state(s) in %r "
                         "for field '%s'." % (field_states, field))

    # Build mapping from one field to the other.
    field_mapping = defaultdict(list)
    for field_state, color_by_state in zip(all_field_states,
                                           all_color_by_states):
        if field_state in field_states:
            field_mapping[field_state].append(color_by_state)

    # For each of the specified input field states, find its corresponding
    # "color by" field state and give it a color if it hasn't been assigned one
    # yet. Make sure we have enough colors and there is a one-to-one mapping.
    color_mapping = {}
    for field_state in field_states:
        color_by_states = set(field_mapping[field_state])

        if len(color_by_states) != 1:
            raise ValueError("The field '%s' to color by does not have a "
                             "one-to-one mapping with field '%s'. Coloring "
                             "would be ambiguous." % (color_by_field, field))

        color_by_state = list(color_by_states)[0]
        if color_by_state not in color_mapping:
            if len(color_pool) > 0:
                color_mapping[color_by_state] = color_pool.pop(0)
            else:
                raise ValueError("There are not enough available QIIME colors "
                                 "to color each of the field states in field "
                                 "'%s'. Coloring would be ambiguous." %
                                 color_by_field)

        colors.append(color_mapping[color_by_state])

    return colors, color_mapping
