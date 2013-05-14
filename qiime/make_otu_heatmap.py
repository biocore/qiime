#!/usr/bin/env python
#file make_otu_heatmap.py

from __future__ import division

__author__ = "Dan Knights"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Dan Knights"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Dan Knights"
__email__ = "daniel.knights@colorado.edu"
__status__ = "Development"


from numpy import array,concatenate,asarray,transpose,log,invert,asarray,\
    float32,float64, unique, fliplr
from cogent.parse.table import SeparatorFormatParser
from optparse import OptionParser
from qiime.util import MissingFileError
import os
from matplotlib import use
use('Agg',warn=False)
import matplotlib
from matplotlib.pylab import *
from qiime.beta_diversity import get_nonphylogenetic_metric
from cogent.core.tree import PhyloNode
from cogent.cluster.UPGMA import UPGMA_cluster
from qiime.parse import parse_newick, PhyloNode
from qiime.filter import filter_samples_from_otu_table

def get_overlapping_samples(map_rows, otu_table):
    """Extracts only samples contained in otu table and mapping file.
    
       Returns: new_map_rows, new_otu_table
    """
    map_sample_ids = zip(*map_rows)[0]
    shared_ids = set(map_sample_ids) & set(otu_table.SampleIds)

    otu_table = filter_samples_from_otu_table(otu_table, shared_ids, 0, inf)

    new_map = []
    for sam_id in map_sample_ids:
        if sam_id in shared_ids:
            ix = map_sample_ids.index(sam_id)
            new_map.append(map_rows[ix])
    
    return new_map, otu_table

def extract_metadata_column(sample_ids, metadata, category):
    """Extracts values from the given metadata column"""
    col_ix = metadata[1].index(category)
    map_sample_ids = zip(*metadata[0])[0]
    category_labels = []
    
    for i,sample_id in enumerate(sample_ids):
        if sample_id in map_sample_ids:
            row_ix = map_sample_ids.index(sample_id)
            entry = metadata[0][row_ix][col_ix]
            category_labels.append(entry)
    return category_labels

def get_order_from_categories(otu_table, category_labels):
    """Groups samples by category values; clusters within each group"""
    category_labels = array(category_labels)
    sample_order = []

    for label in unique(category_labels):
        label_ix = category_labels==label
        selected = [s for (i,s) in zip(label_ix, otu_table.SampleIds) if i]
        sub_otu_table = filter_samples_from_otu_table(otu_table, selected, 0, inf)
        data = asarray([val for val in sub_otu_table.iterObservationData()])
        label_ix_ix = get_clusters(data, axis='column')

        sample_order += list(nonzero(label_ix)[0][array(label_ix_ix)])
    return array(sample_order)


def get_order_from_tree(ids, tree_text):
    """Returns the indices that would sort ids by tree tip order"""
    tree = parse_newick(tree_text, PhyloNode) 
    ordered_ids = []
    for tip in tree.iterTips():
        if tip.Name in ids:
            ordered_ids.append(tip.Name)
    return names_to_indices(ids, ordered_ids)

def make_otu_labels(otu_ids, lineages, n_levels=1):
    """Returns 'pretty' OTU labels: 'Lineage substring (OTU ID)'
        
       Lineage substring includes the last n_levels lineage levels 
    """

    if len(lineages[0]) > 0:
        otu_labels = []
        for i, lineage in enumerate(lineages):
            if n_levels > len(lineage):
                otu_label = '%s (%s)' %(';'.join(lineage),otu_ids[i])
            else:
                otu_label = '%s (%s)' \
                    %(';'.join(lineage[-n_levels:]), otu_ids[i])
            otu_labels.append(otu_label)
        otu_labels = [lab.replace('"','') for lab in otu_labels]
    else:
        otu_labels = otu_ids
    return otu_labels

def names_to_indices(names, ordered_names):
    """Returns the indices that would sort 'names' like 'ordered_names'
    """
    indices = []
    names_list = list(names)
    for ordered_name in ordered_names:
        if ordered_name in names_list:
            indices.append(names_list.index(ordered_name))
    return array(indices)
    

def get_log_transform(otu_table, eps=None):
    """Returns log10 of the data, setting zero values to eps.
    
       If eps is None, eps is set to 1/2 the smallest nonzero value.
    """
    ## NOTE: compared with qiime.make_otu_heatmap_html, this function does
    ##  *not* do a data = data - (data).min() before returning value.
    ##  This behavior is correct according to Dan Knights, but consider if
    ##  we ever merge these two scripts

    # explicit conversion to float: transform
    def f(s_v, s_id, s_md):
        return float64(s_v)
    float_otu_table = otu_table.transformSamples(f)

    if eps is None:
        # get the minimum among nonzero entries and divide by two
        eps = inf
        for (obs, sam) in float_otu_table.nonzero():
            eps = minimum(eps, float_otu_table.getValueByIds(obs,sam))
        if eps == inf:
            raise ValueError('All values in the OTU table are zero!')
        else:
            eps = eps / 2

    # set zero entries to eps/2 using a transform
    g = lambda x : x if (x != 0) else eps
    def g_m(s_v, s_id, s_md):
        return asarray(map(g,s_v))

    eps_otu_table = float_otu_table.transformSamples(g_m)

    # take log of all values
    def h(s_v, s_id, s_md):
        return log10(s_v)
    log_otu_table = eps_otu_table.transformSamples(h)

    return log_otu_table
    
def get_clusters(x_original, axis=['row','column'][0]):
    """Performs UPGMA clustering using euclidean distances"""
    x = x_original.copy()
    if axis=='column':
        x = x.T
    nr = x.shape[0]
    metric_f = get_nonphylogenetic_metric('euclidean')
    row_dissims = metric_f(x)
    # do upgma - rows
    BIG = 1e305
    row_nodes = map(PhyloNode, map(str,range(nr)))
    for i in range(len(row_dissims)):
        row_dissims[i,i] = BIG
    row_tree = UPGMA_cluster(row_dissims, row_nodes, BIG)
    row_order = [int(tip.Name) for tip in row_tree.iterTips()]
    return row_order


def get_fontsize(numrows):
    """Returns the fontsize needed to make text fit within each row.
    """
    thresholds = [25, 50, 75, 100, 125]
    sizes =      [ 5,  4,  3,   2, 1.5,   1]
    i = 0
    while numrows > thresholds[i]:
        i += 1
        if i == len(thresholds):
            break
    return sizes[i]


def plot_heatmap(otu_table, row_labels, col_labels, filename='heatmap.pdf',
        width=5, height=5, textborder=.25):
    """Create a heatmap plot, save as a pdf.
    
        'width', 'height' are in inches
        
        'textborder' is the fraction of the figure allocated for the
        tick labels on the x and y axes
    """
    nrow = len(otu_table.ObservationIds)
    ncol = len(otu_table.SampleIds)

    # determine appropriate font sizes for tick labels
    row_fontsize = get_fontsize(nrow)
    col_fontsize = get_fontsize(ncol)

    # create figure and plot heatmap
    fig = figure(figsize=(width, height))
    my_cmap=get_cmap('gist_gray')
    # numpy magic: [:,::-1] actually means fliplr()
    #imshow(x[:,::-1],interpolation='nearest', aspect='auto', cmap=my_cmap)

    data = [val for val in otu_table.iterObservationData()]
    imshow(fliplr(data),interpolation='nearest', aspect='auto', cmap=my_cmap)
    ax = fig.axes[0]

    # imshow is offset by .5 for some reason
    xlim(-.5, ncol-.5)
    ylim(-.5, nrow-.5)
    
    # add ticklabels to axes
    xticks(arange(ncol), col_labels[::-1], fontsize=col_fontsize)
    yticks(arange(nrow), row_labels, fontsize=row_fontsize)

    # turn off tick marks
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')

    # rotate x ticklabels
    for label in ax.xaxis.get_ticklabels():
        label.set_rotation(90)

    # add space for tick labels
    fig.subplots_adjust(left=textborder, bottom=textborder)
    cb = colorbar() # grab the Colorbar instance
    # set colorbar tick labels to a reasonable value (normal is large)
    for t in cb.ax.get_yticklabels():
        t.set_fontsize(5)
    fig.savefig(filename)
