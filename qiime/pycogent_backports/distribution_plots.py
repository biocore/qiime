#!/usr/bin/env python 

"""This module contains functions for plotting distributions in various ways.

There are two different types of plotting functions:

generate_box_plots() plots several boxplots next to each other for easy
comparison.

generate_comparative_plots() plots groupings of distributions at data
points along the x-axis.
"""

from math import isnan
from matplotlib import use
use('Agg', warn=False)
from matplotlib.colors import colorConverter
from matplotlib.patches import Polygon, Rectangle
from matplotlib.pyplot import boxplot, figure
from matplotlib.transforms import Bbox
from numpy import array, mean, random, std

__author__ = "Jai Rideout"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jai Rideout"]
__license__ = "GPL"
__version__ = "1.3.0dev"
__maintainer__ = "Jai Rideout"
__email__ = "jr378@nau.edu"
__status__ = "Development"

def generate_box_plots(distributions, x_values=None, x_tick_labels=None,
                       title=None, x_label=None, y_label=None,
                       x_tick_labels_orientation='vertical', y_min=None,
                       y_max=None):
    """Returns a matplotlib.figure.Figure object containing a boxplot for each
    distribution.

    Arguments:
        - distributions: A list of lists containing each distribution.
        - x_values: A list indicating where each boxplot should be placed. Must
            be the same length as distributions if provided.
        - x_tick_labels: A list of labels to be used to label x-axis ticks.
        - title: A string containing the title of the plot.
        - x_label: A string containing the x-axis label.
        - y_label: A string containing the y-axis label.
        - x_tick_labels_orientation: A string specifying the orientation of the
            x-axis labels (either "vertical" or "horizontal").
        - y_min: The minimum value of the y-axis. If None, uses matplotlib's
            autoscale.
        - y_max: The maximum value of the y-axis. If None, uses matplotlib's
            autoscale.
    """
    # Make sure our input makes sense.
    for distribution in distributions:
        if len(distribution) == 0:
            raise ValueError("Some of the provided distributions are empty.")
        try:
            map(float, distribution)
        except:
            raise ValueError("Each value in each distribution must be a "
                             "number.")

    _validate_x_values(x_values, x_tick_labels, len(distributions));

    # Create a new figure to plot our data on, and then plot the distributions.
    result, plot_axes = _create_plot()
    boxplot(distributions, positions=x_values)

    # Set up the various plotting options, such as x- and y-axis labels, plot
    # title, and x-axis values if they have been supplied.
    _set_axes_options(plot_axes, title, x_label, y_label,
                      x_tick_labels=x_tick_labels,
                      x_tick_labels_orientation=x_tick_labels_orientation,
                      y_min=y_min, y_max=y_max)

    # We need to explicitly draw the figure before returning it because
    # certain rendering backends don't trigger the "draw_event" for
    # certain file formats (e.g. svg and pdf). This event needs to be
    # fired off so that the plots can automatically resize themselves
    # so that tick labels don't get cut off.
    result.canvas.draw()
    return result

def generate_comparative_plots(plot_type, data, x_values=None,
        data_point_labels=None, distribution_labels=None,
        distribution_markers=None, x_label=None, y_label=None, title=None,
        x_tick_labels_orientation='vertical', y_min=None, y_max=None):
    """Returns a matplotlib.figure.Figure object containing plots of the
    specified type grouped at points along the x-axis.

    Arguments:
        - plot_type: A string indicating what type of plot should be created.
            Can be one of "bar", "scatter", or "box", where "bar" is a bar
            chart, "scatter" is a scatter plot, and "box" is a box plot.
        - data: A list of lists that represent each data point along the
            x-axis. Each data point contains lists of data for each
            distribution in the group at that point. This nesting allows for
            the grouping of distributions at each data point.
        - x_values: A list indicating the spacing along the x-axis. Must
            be the same length as the number of data points if provided. If not
            provided, plots will be spaced evenly.
        - data_point_labels: A list of strings containing the label for each
            data point.
        - distribution_labels: A list of strings containing the label for each
            distribution in a data point grouping.
        - distribution_markers: A list of matplotlib-compatible strings or
            tuples that indicate the color or symbol to be used to distinguish
            each distribution in a data point grouping. Colors will be used for
            bar charts or box plots, while markers are used for scatter plots.
        - x_label: A string containing the x-axis label.
        - y_label: A string containing the y-axis label.
        - title: A string containing the title of the plot.
        - x_tick_labels_orientation: A string specifying the orientation of the
            x-axis labels (either "vertical" or "horizontal").
        - y_min: The minimum value of the y-axis. If None, uses matplotlib's
            autoscale.
        - y_max: The maximum value of the y-axis. If None, uses matplotlib's
            autoscale.
    """
    # Set up different behavior based on the plot type.
    if plot_type == 'bar':
        plotting_function = _plot_bar_data
        legend_function = _create_standard_legend
        distribution_centered = False
        colors = True
    elif plot_type == 'scatter':
        plotting_function = _plot_scatter_data
        legend_function = _create_standard_legend
        distribution_centered = True
        colors = False
    elif plot_type == 'box':
        plotting_function = _plot_box_data
        legend_function = _create_box_plot_legend
        distribution_centered = True
        colors = True
    else:
        raise ValueError("Invalid plot type '%s'. Supported plot types are "
                "'bar', 'scatter', or 'box'." % plot_type)

    num_points, num_distributions = _validate_input(data, x_values,
            data_point_labels, distribution_labels)

    # Create a list of matplotlib markers (colors or symbols) that can be used
    # to distinguish each of the distributions.
    if distribution_markers is None:
        distribution_markers = _get_enumerated_values(colors)

    if len(distribution_markers) < num_distributions:
        # We don't have enough markers to represent each distribution uniquely,
        # so let the user know that they need to provide enough markers for
        # their large group size.
        raise ValueError("There are not enough markers to uniquely represent "
                         "each distribution in your dataset. Please provide a "
                         "list of markers that is at least as large as the "
                         "number of distributions in your dataset.")

    # Keep track of how wide each distribution is and how much spacing is put
    # between each grouping of distributions on the horizontal axis.
    dist_width = 0.4
    group_spacing = 0.5

    # Now calculate where each of the data points will start on the x-axis.
    x_locations = _calc_data_point_locations(x_values, num_points,
            num_distributions, dist_width, group_spacing)
    assert(len(x_locations) == num_points)

    # Create the figure to put the plots on, as well as a list to store each of
    # the plots (needed for the legend).
    result, plot_axes = _create_plot()
    plots = []

    # Iterate over each data point, and plot each of the distributions at that
    # data point. Increase the offset after each distribution is plotted,
    # so that the grouped distributions don't overlap.
    for point, x_pos in zip(data, x_locations):
        dist_offset = 0
        for dist, dist_marker in zip(point, distribution_markers):
            dist_location = x_pos + dist_offset
            plots.append(plotting_function(plot_axes, dist, dist_marker,
                                           dist_width, dist_location))
            dist_offset += dist_width

    # Set up various plot options that are best set after the plotting is done.
    # The x-axis tick marks (one per data point) are centered on each group of
    # distributions.
    plot_axes.set_xticks(_calc_data_point_ticks(x_locations,
            num_distributions, dist_width, distribution_centered))
    _set_axes_options(plot_axes, title, x_label, y_label, x_values,
                      data_point_labels, x_tick_labels_orientation, y_min,
                      y_max)

    # Add a legend for the different distribution markers.
    if distribution_labels is not None:
        legend_function(plots, plot_axes, distribution_markers,
                        num_distributions, distribution_labels)
    result.canvas.draw()

    # matplotlib seems to sometimes plot points on the rightmost edge of the
    # plot without adding padding, so we need to add our own to both sides of
    # the plot. For some reason this has to go after the call to draw(),
    # otherwise matplotlib throws an exception saying it doesn't have a
    # renderer.
    plot_axes.set_xlim(plot_axes.get_xlim()[0] - group_spacing,
                       plot_axes.get_xlim()[1] + group_spacing)
    return result

def _validate_input(data, x_values, data_point_labels, distribution_labels):
    """Returns a tuple containing the number of data points and distributions
    in the data.
    
    Validates plotting options to make sure they are valid with the supplied
    data.
    """
    if data is None or not data or isinstance(data, basestring):
        raise ValueError("The data must be a list type, and it cannot be "
                         "None or empty.")

    num_points = len(data)
    num_distributions = len(data[0])

    empty_data_error_msg = "The data must contain at least one data " + \
                           "point, and each data point must contain at " + \
                           "least one distribution to plot."
    if num_points == 0 or num_distributions == 0:
        raise ValueError(empty_data_error_msg)

    for point in data:
        if len(point) == 0:
            raise ValueError(empty_data_error_msg)
        if len(point) != num_distributions:
            raise ValueError("The number of distributions in each data point "
                             "grouping must be the same for all data points.")

    # Make sure we have the right number of x values (one for each data point),
    # and make sure they are numbers.
    _validate_x_values(x_values, data_point_labels, num_points)

    if (distribution_labels is not None and
        len(distribution_labels) != num_distributions):
        raise ValueError("The number of distribution labels must be equal "
                         "to the number of distributions.")
    return num_points, num_distributions

def _validate_x_values(x_values, x_tick_labels, num_expected_values):
    """Validates the x values provided by the user, making sure they are the
    correct length and are all numbers.
    
    Also validates the number of x-axis tick labels.
    
    Raises a ValueError if these conditions are not met.
    """
    if x_values is not None:
        if len(x_values) != num_expected_values:
            raise ValueError("The number of x values must match the number "
                             "of data points.")
        try:
            map(float, x_values)
        except:
            raise ValueError("Each x value must be a number.")

    if x_tick_labels is not None:
        if len(x_tick_labels) != num_expected_values:
            raise ValueError("The number of x-axis tick labels must match the "
                             "number of data points.")

def _get_enumerated_values(colors=True):
    """Returns a list of valid matplotlib colors or symbols.

    If colors is False, a list of symbols is returned.
    """
    return ['b', 'g', 'r', 'c', 'm', 'y', 'w'] if colors else\
           ['s', 'o', '^', '>', 'v', '<', 'd', 'p', 'h', '8', '+', 'x']

def _calc_data_point_locations(x_values, num_points, num_distributions,
                               dist_width, group_spacing):
    """Returns a numpy array of x-axis locations for each of the data points to
    start at.

    Note: A numpy array is returned so that the overloaded "+" operator can be
    used on the array.
    
    The x locations are spaced according to the spacing between points, and the
    width of each distribution grouping at each point. The x locations are also
    scaled by the x_values that may have been supplied by the user. If none are
    supplied, the x locations are evenly spaced.
    """
    if dist_width <= 0 or group_spacing < 0:
        raise ValueError("The width of a distribution cannot be zero or "
                         "negative. The width of the spacing between groups "
                         "of distributions cannot be negative.")
    if x_values is None:
        # Evenly space the x locations.
        x_values = range(1, num_points + 1)

    assert(len(x_values) == num_points)

    # Calculate the width of each grouping of distributions at a data point.
    # This is multiplied by the current x value to give us our final
    # absolute horizontal position for the current point.
    return array([(dist_width * num_distributions + group_spacing) * x_val\
                     for x_val in x_values])

def _calc_data_point_ticks(x_locations, num_distributions, distribution_width,
                           distribution_centered):
    """Returns a 1D numpy array of x-axis tick positions.

    These positions will be centered on each data point.
    
    Set distribution_centered to True for scatter and box plots because their
    plot types naturally center over a given horizontal position. Bar charts 
    should use distribution_centered = False because the leftmost edge of a bar
    starts at a given horizontal position and extends to the right for the
    width of the bar.
    """
    dist_size = num_distributions - 1 if distribution_centered else\
            num_distributions
    return x_locations + ((dist_size * distribution_width) / 2)

def _create_plot():
    """Creates a plot and returns the matplotlib Figure and Axes objects
    associated with it.
    
    The _on_draw() event handler is registered with the figure's canvas, as
    this is needed to auto-resize the plot to encompass text labels that might
    otherwise get cut off when the plot is rendered.
    """
    fig = figure()
    ax = fig.add_subplot(111)
    fig.canvas.mpl_connect('draw_event', _on_draw)
    return fig, ax

def _plot_bar_data(plot_axes, distribution, distribution_color,
                   distribution_width, x_position):
    """Returns the result of plotting a single bar in matplotlib."""
    result = None
    avg = mean(distribution)
    std_dev = std(distribution)
    if not isnan(avg) and not isnan(std_dev):
        # numpy's mean() and std() functions will return NaN for empty
        # lists of data, and we do not want to plot these because
        # matplotlib will not be able to render them as PDFs.
        result = plot_axes.bar(x_position, avg, distribution_width,
                               yerr=std_dev, ecolor='black',
                               facecolor=distribution_color)
    return result

def _plot_scatter_data(plot_axes, distribution, distribution_symbol,
                       distribution_width, x_position):
    """Returns the result of plotting a single scatterplot in matplotlib."""
    result = None
    x_vals = [x_position] * len(distribution)
    # matplotlib's scatter function doesn't like plotting empty data.
    if len(x_vals) > 0 and len(distribution) > 0:
        result = plot_axes.scatter(x_vals, distribution,
                                   marker=distribution_symbol, c='k')
    return result

def _plot_box_data(plot_axes, distribution, distribution_color,
                   distribution_width, x_position):
    """Returns the result of plotting a single boxplot in matplotlib."""
    box_plot = plot_axes.boxplot([distribution], positions=[x_position],
                                 widths=distribution_width)

    # Fill the distribution's box with the desired color.
    # Note: the following boxplot-coloring code is largely taken from a
    # matplotlib boxplot example:
    # http://matplotlib.sourceforge.net/examples/pylab_examples/
    #     boxplot_demo2.html
    num_boxes = len(box_plot['boxes'])
    for box_num in range(num_boxes):
        box = box_plot['boxes'][box_num]
        boxX = []
        boxY = []
        # There are five points in the box. The first is the same as
        # the last.
        for j in range(5):
            boxX.append(box.get_xdata()[j])
            boxY.append(box.get_ydata()[j])
        boxCoords = zip(boxX,boxY)
        boxPolygon = Polygon(boxCoords, facecolor=distribution_color)
        plot_axes.add_patch(boxPolygon)

        # Draw the median lines back over what we just filled in with
        # color.
        median = box_plot['medians'][box_num]
        medianX = []
        medianY = []
        for j in range(2):
            medianX.append(median.get_xdata()[j])
            medianY.append(median.get_ydata()[j])
            plot_axes.plot(medianX, medianY, 'black')
    return box_plot

def _set_axes_options(plot_axes, title=None, x_label=None, y_label=None,
                      x_values=None, x_tick_labels=None,
                      x_tick_labels_orientation='vertical', y_min=None,
                      y_max=None):
    """Applies various labelling options to the plot axes."""
    if title is not None:
        plot_axes.set_title(title)
    if x_label is not None:
        plot_axes.set_xlabel(x_label)
    if y_label is not None:
        plot_axes.set_ylabel(y_label)

    if (x_tick_labels_orientation != 'vertical' and
        x_tick_labels_orientation != 'horizontal'):
        raise ValueError("Invalid orientation for x-axis tick labels: %s. "
                         "Valid orientations are 'vertical' or 'horizontal'."
                         % x_tick_labels_rotation)

    # If labels are provided, always use them. If they aren't, use the x_values
    # that denote the spacing between data points as labels. If that isn't
    # available, simply label the data points in an incremental fashion,
    # i.e. 1, 2, 3,...,n, where n is the number of data points on the plot.
    if x_tick_labels is not None:
        labels = plot_axes.set_xticklabels(x_tick_labels,
                                           rotation=x_tick_labels_orientation)
    elif x_tick_labels is None and x_values is not None:
        labels = plot_axes.set_xticklabels(x_values,
                                           rotation=x_tick_labels_orientation)
    else:
        labels = plot_axes.set_xticklabels(
                    range(1, len(plot_axes.get_xticklabels()) + 1),
                    rotation=x_tick_labels_orientation)

    # Set the y-axis range if specified.
    if y_min is not None:
        plot_axes.set_ylim(bottom=float(y_min))
    if y_max is not None:
        plot_axes.set_ylim(top=float(y_max))

def _create_standard_legend(plots, plot_axes, distribution_markers,
                            num_distributions, distribution_labels):
    """Creates a default matplotlib legend on the supplied axes."""
    # We need to let matplotlib know of the first group of distributions that
    # we plotted, as well as their labels.
    plot_axes.legend(plots[:num_distributions], distribution_labels,
                     loc='best')

def _create_box_plot_legend(plots, plot_axes, distribution_colors,
                            num_distributions, distribution_labels):
    """Creates a custom matplotlib legend on the supplied axes.
    
    This function is useful for boxplots or other plots that do not have a
    useful default legend in matplotlib. It creates rectangles of different
    colors and pairs them up with the supplied labels.
    """
    # We have to fake the legend because box plots currently don't have a very
    # useful legend in matplotlib. Note: This code was taken from a matplotlib
    # example:
    # http://matplotlib.sourceforge.net/users/legend_guide.html
    legend_rects = [Rectangle((0, 0), 1, 1, fc=color) for color in
            distribution_colors[:num_distributions]]
    assert(len(legend_rects) == len(distribution_labels))
    plot_axes.legend(legend_rects, distribution_labels, loc='best')

def _on_draw(event):
    """Makes room for vertical x-axis tick labels on the plot (if necessary) so
    that they don't get cut off when the plot is rendered.
    
    This function must be connected to matplotlib's draw_event in order to be
    called.
    """
    # Note: the following code is largely taken from two separate code
    # postings. The first reference dealt with text labels being cut off when
    # the plot was rendered, and the second showed how to prevent infinite
    # recursion when handling a draw event that draws.
    # (1) http://matplotlib.sourceforge.net/faq/howto_faq.html#automatically-
    #       make-room-for-tick-labels
    # (2) http://stackoverflow.com/questions/4018860/text-box-in-matplotlib/
    #       4056853#4056853
    fig = event.canvas.figure
    for ax in fig.axes:
        x_tick_labels = ax.get_xticklabels()
        bboxes = []

        for label in x_tick_labels:
            bbox = label.get_window_extent()
            # The figure transform goes from relative coords->pixels and we
            # want the inverse of that.
            bboxi = bbox.inverse_transformed(fig.transFigure)
            bboxes.append(bboxi)

        # This is the bbox that bounds all the bboxes, again in relative
        # figure coords.
        bbox = Bbox.union(bboxes)

        if fig.subplotpars.bottom < bbox.height:
            # We need to move it up to add some padding for the labels.
            fig.subplots_adjust(bottom=1.2*bbox.height)

            # Temporarily disconnect any callbacks to the draw event to avoid
            # recursion.
            func_handles = fig.canvas.callbacks.callbacks[event.name]
            fig.canvas.callbacks.callbacks[event.name] = {}

            # Re-draw the figure and reset the draw event callbacks.
            fig.canvas.draw()
            fig.canvas.callbacks.callbacks[event.name] = func_handles
    return False
