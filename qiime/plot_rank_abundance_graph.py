#!/usr/bin/env python
# File created on 17 Aug 2010
from __future__ import division

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Jens Reeder"]
__license__ = "GPL"
__version__ = "1.2.0"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Release"

from matplotlib import use
use('Agg')

from numpy import arange,  array
from itertools import cycle
from matplotlib.pyplot import plot, gca,  ylim, xlim, show, legend, \
    savefig

from qiime.parse import parse_otu_table
from qiime.colors import data_color_order

def make_sorted_frequencies(counts, absolute=False):
    """transform and sort a vector of count.

    counts: a column of an OTU table
    absolute: if True return absolute values instead of frequencies.
    """

    c = counts
    c.sort()
    c = c[c.nonzero()]
    c = c[::-1]
    if absolute:
        return c
    else:
        f = c/float(c.sum())
        return f

def plot_rank_abundance_graph(counts, color='red', absolute=False, label=None):
    """Plots rank-abundance curve.

    counts: a column of an OTU table
    color: color of the series to plot
    absolute: if True plot absolute counts instead of  freqs
    label: text for the legend of this series
    """

    f= make_sorted_frequencies(counts, absolute)
    x = arange(1,len(f)+1) 
    plot(x, f, color=color, alpha= 0.8, label=label)
    ax = gca()
    return ax


def plot_rank_abundance_graphs(sample_names, otu_table_fh,
                               output_dir, file_type='pdf',
                               absolute_counts=False,
                               x_linear_scale=False,
                               y_linear_scale=False,
                               no_legend=False,
                               log_fh=None):

    """plot rank-abundance curves for sample specified in sample_name.
    
    sample_names: comma separated string of sample names
    otu_table_fh: open file handle to otu table
    output_dir: existing directory to which files are written
    file_type: valid matplotlib file type
    x_linear_scale: if True draw x axis in linear scale, otherwise use log
    y_linear_scale: if True draw y axis in linear scale, otherwise use log
    no_legend: if True don't draw legend
    log_fh: open file handle to log file, if not None used to log 
"""
    sample_ids, otu_ids, otu_table, lineages = parse_otu_table(otu_table_fh)

    #figure out which samples to draw
    if sample_names=='*':
        user_sample_names = sample_ids
    else:
        user_sample_names = sample_names.split(',')
        if len(user_sample_names)<1:
            raise ValueError, "sample IDs must be comma separated list of "\
            +"sample names - found %s" % sample_names 

    # do the actual drawing
    ax=None
    for sample_name,color in zip(user_sample_names, cycle(data_color_order)):
        try:
            index = sample_ids.index(sample_name)
        except ValueError:
            if log_fh:
                log_fh.write("Warning: Sample name %s not in OTU table - skipping." % sample_name)
            continue     
        ax = plot_rank_abundance_graph(otu_table[:,index], color=color,
                                       absolute=absolute_counts,
                                       label=sample_name)
        ax.set_label(sample_name)
    
    if ax==None:
        #ax should be defined if at least one series has been drawn
        raise ValueError("No data series drawn. Check your OTU table and sample names")

    #settings for all series
    ax.grid()      
    ax.set_xlabel('Species rank')
    ax.set_ylabel('Relative abundance')

    if not x_linear_scale:
        ax.set_xscale('log')
    if not y_linear_scale:
        ax.set_yscale('log')
  
    if not no_legend:
        legend()

    #build output fp, if less than MAX_SAMPLES_... append the samples names    
    output_fp = output_dir+ "/rank_abundance"
    MAX_SAMPLES_TO_SHOW_IN_FILENAME = 6
    if len(user_sample_names) < MAX_SAMPLES_TO_SHOW_IN_FILENAME:
        output_fp += '_'.join(user_sample_names) 
    output_fp += ".%s" % file_type

    savefig(output_fp, format=file_type)
