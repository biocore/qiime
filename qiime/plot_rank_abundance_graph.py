#!/usr/bin/env python
# File created on 17 Aug 2010
from __future__ import division

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jens Reeder", "Emily TerAvest"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

from matplotlib import use
use('Agg', warn=False)

from numpy import arange, array, sort
from itertools import cycle
from matplotlib.pyplot import plot, gca,  ylim, xlim, show, legend, \
    savefig
from os.path import join
from qiime.colors import data_color_order, data_colors
from biom.table import UnknownIDError


def make_sorted_frequencies(counts, absolute=False):
    """transform and sort a vector of count.

    counts: a column of an OTU table
    absolute: if True return absolute values instead of frequencies.
    """

    c = sort(counts)
    c = c[c.nonzero()]
    c = c[::-1]
    if absolute:
        return c
    else:
        f = c / float(c.sum())
        return f


def plot_rank_abundance_graph(
        otu_count_vector, color='red', absolute=False, label=None):
    """Plots rank-abundance curve.

    otu_count_vector: a vector of otu counts for a single sample
    color: color of the series to plot
    absolute: if True plot absolute counts instead of  freqs
    label: text for the legend of this series
    """

    f = make_sorted_frequencies(otu_count_vector, absolute)
    x = arange(1, len(f) + 1)
    plot(x, f, color=color, alpha=0.8, label=label)
    ax = gca()
    return ax


def plot_rank_abundance_graphs(result_fp, sample_names, otu_table,
                               file_type='pdf',
                               absolute_counts=False,
                               x_linear_scale=False,
                               y_linear_scale=False,
                               no_legend=False,
                               log_fh=None):
    """plot rank-abundance curves for sample specified in sample_name.

    result_fp: filename of output figure
    sample_names: comma separated string of sample names
    otu_table_fh: open file handle to otu table
    file_type: valid matplotlib file type
    x_linear_scale: if True draw x axis in linear scale, otherwise use log
    y_linear_scale: if True draw y axis in linear scale, otherwise use log
    no_legend: if True don't draw legend
    log_fh: open file handle to log file, if not None used to log
"""
    # figure out which samples to draw
    if sample_names == '*':
        user_sample_names = otu_table.ids()
    else:
        user_sample_names = sample_names.split(',')
        if len(user_sample_names) < 1:
            raise ValueError("sample IDs must be comma separated list of "
                             + "sample names - found %s" % sample_names)

    # do the actual drawing
    ax = None
    for sample_name, color in zip(user_sample_names, cycle(data_color_order)):
        color = data_colors[color].toHex()
        try:
            otu_count_vector = otu_table.data(sample_name, 'sample')
        except UnknownIDError:
            if log_fh:
                log_fh.write(
                    "UnknownIDError: Sample name %s not in OTU table - skipping." %
                    sample_name)
            continue

        ax = plot_rank_abundance_graph(otu_count_vector,
                                       color=color,
                                       absolute=absolute_counts,
                                       label=sample_name)
        ax.set_label(sample_name)

    if ax is None:
        # ax should be defined if at least one series has been drawn
        raise ValueError(
            "No data series drawn. Check your OTU table and sample names")

    # settings for all series
    ax.grid()
    ax.set_xlabel('Species rank')
    if absolute_counts:
        ax.set_ylabel('Absolute abundance')
    else:
        ax.set_ylabel('Relative abundance')

    if not x_linear_scale:
        ax.set_xscale('log')
    if not y_linear_scale:
        ax.set_yscale('log')

    if not no_legend:
        legend()

    if not result_fp.endswith(file_type):
        result_fp += '.' + file_type
    savefig(result_fp)
