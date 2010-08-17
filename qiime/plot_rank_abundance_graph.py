#!/usr/bin/env python
# File created on 17 Aug 2010
from __future__ import division

__author__ = "Jens Reeder"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Jens Reeder"]
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Jens Reeder"
__email__ = "jens.reeder@gmail.com"
__status__ = "Development"

from numpy import arange,  array
from matplotlib.pyplot import plot, gca

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
