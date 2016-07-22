from __future__ import division

__author__ = "Dan Knights"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Dan Knights", "Greg Caporaso", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Dan Knights"
__email__ = "daniel.knights@colorado.edu"

import numpy as np
import matplotlib
matplotlib.use('Agg', warn=False)
import matplotlib.pyplot as plt
from mpl_toolkits.axes_grid1 import make_axes_locatable
from scipy.cluster.hierarchy import linkage
from skbio.tree import TreeNode
from skbio.diversity.beta import pw_distances

from qiime.parse import parse_newick, PhyloNode
from qiime.filter import filter_samples_from_otu_table


def get_overlapping_samples(map_rows, otu_table):
    """Extracts only samples contained in otu table and mapping file.

       Returns: new_map_rows, new_otu_table
    """
    map_sample_ids = zip(*map_rows)[0]
    shared_ids = set(map_sample_ids) & set(otu_table.ids())

    otu_table = filter_samples_from_otu_table(otu_table, shared_ids, -np.inf,
                                              np.inf)

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

    for i, sample_id in enumerate(sample_ids):
        if sample_id in map_sample_ids:
            row_ix = map_sample_ids.index(sample_id)
            entry = metadata[0][row_ix][col_ix]
            category_labels.append(entry)
    return category_labels


def get_order_from_categories(otu_table, category_labels):
    """Groups samples by category values; clusters within each group"""
    category_labels = np.array(category_labels)
    sample_order = []

    for label in np.unique(category_labels):
        label_ix = category_labels == label
        selected = [s for (i, s) in zip(label_ix, otu_table.ids()) if i]
        sub_otu_table = filter_samples_from_otu_table(otu_table, selected,
                                                      -np.inf, np.inf)
        data = np.asarray(list(sub_otu_table.iter_data(axis='observation')))
        label_ix_ix = get_clusters(data, axis='column')

        sample_order += list(np.nonzero(label_ix)[0][np.array(label_ix_ix)])
    return np.array(sample_order)


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
                otu_label = '%s (%s)' % (';'.join(lineage), otu_ids[i])
            else:
                otu_label = '%s (%s)' \
                    % (';'.join(lineage[-n_levels:]), otu_ids[i])
            otu_labels.append(otu_label)
        otu_labels = [lab.replace('"', '') for lab in otu_labels]
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
    return np.array(indices)


def get_log_transform(otu_table):
    """Returns log10 of the data"""
    if otu_table.nnz == 0:
        raise ValueError('All values in the OTU table are zero!')

    # take log of all values
    def h(s_v, s_id, s_md):
        return np.log10(s_v)

    return otu_table.transform(h, axis='sample', inplace=False)


def get_clusters(x_original, axis='row'):
    """Performs UPGMA clustering using euclidean distances"""
    x = x_original.copy()
    if axis == 'column':
        x = x.T
    nr = x.shape[0]
    row_dissims = pw_distances(x, ids=map(str, range(nr)), metric='euclidean')
    # do upgma - rows
    # Average in SciPy's cluster.hierarchy.linkage is UPGMA
    linkage_matrix = linkage(row_dissims.condensed_form(), method='average')
    tree = TreeNode.from_linkage_matrix(linkage_matrix, row_dissims.ids)
    return [int(tip.name) for tip in tree.tips()]


def get_fontsize(numrows):
    """Returns the fontsize needed to make text fit within each row.
    """
    thresholds = [25, 50, 75, 100, 125]
    sizes = [5, 4, 3, 2, 1.5, 1]
    i = 0
    while numrows > thresholds[i]:
        i += 1
        if i == len(thresholds):
            break
    return sizes[i]


def plot_heatmap(otu_table, row_labels, col_labels, filename, imagetype='pdf',
                 width=5, height=5, dpi=None, textborder=.25,
                 color_scheme='YlGn'):
    """Create a heatmap plot, save as a pdf by default.

        'width', 'height' are in inches

        'textborder' is the fraction of the figure allocated for the
        tick labels on the x and y axes

        color_scheme: choices can be found at
         http://matplotlib.org/examples/color/colormaps_reference.html
    """
    nrow = otu_table.length(axis='observation')
    ncol = otu_table.length(axis='sample')

    # determine appropriate font sizes for tick labels
    row_fontsize = get_fontsize(nrow)
    col_fontsize = get_fontsize(ncol)

    # create figure and plot heatmap
    fig, ax = plt.subplots(figsize=(width, height))
    data = list(otu_table.iter_data(axis='observation'))
    im = plt.imshow(np.fliplr(data), interpolation='nearest', aspect='auto',
                    cmap=color_scheme)

    # imshow is offset by .5 for some reason
    plt.xlim(-.5, ncol - .5)
    plt.ylim(-.5, nrow - .5)

    # add ticklabels to axes
    plt.xticks(np.arange(ncol), col_labels[::-1], fontsize=col_fontsize,
               rotation=90)
    plt.yticks(np.arange(nrow), row_labels, fontsize=row_fontsize)

    # turn off tick marks
    ax.xaxis.set_ticks_position('none')
    ax.yaxis.set_ticks_position('none')

    # add space for tick labels
    fig.subplots_adjust(left=textborder, bottom=textborder)

    # create colorbar (legend) in its own axes so that tight_layout will
    # respect both the heatmap and colorbar when it makes room for everything.
    # code based on example in:
    #     http://matplotlib.org/users/tight_layout_guide.html
    divider = make_axes_locatable(ax)
    cax = divider.append_axes("right", "5%", pad="3%")
    cb = plt.colorbar(im, cax=cax)

    # set colorbar tick labels to a reasonable value (normal is large)
    for t in cb.ax.get_yticklabels():
        t.set_fontsize(5)

    plt.tight_layout()
    fig.savefig(filename, format=imagetype, dpi=dpi)
