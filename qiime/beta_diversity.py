#!/usr/bin/env python

__author__ = "Justin Kuzynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Rob Knight", "Greg Caporaso",
               "Jose Carlos Clemente Litran", "Jai Ram Rideout",
               "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

"""Contains code for performing beta diversity analyses of different samples.

command line usage help: python beta_diversity.py -h

This module has the responsibility for relating samples to one another. This
is performed using the following inputs:
    - a tab-delimited table of taxon x sample counts (taxa can be OTUs).
    - future: will be able to supply
      data in (taxon, sample, count) format, sparse representation is
      important for large datasets where full table is impractical to build.
      Note that this requires adapting the beta diversity metrics to work
      with sparse matrices, not yet done.

The output is a sample x sample matrix of distances, incl. row/col headers.
    Note that parser expects first field to be blank, i.e. first char of file
    is expected to be a tab.
"""

from StringIO import StringIO
from sys import exit, stderr
import os.path
import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')

from numpy import asarray
import cogent.maths.distance_transform as distance_transform
from biom import load_table

from qiime.util import (FunctionWithParams, TreeMissingError,
                        OtuMissingError)
from qiime.format import format_matrix, format_distance_matrix
from qiime.parse import parse_newick, PhyloNode
import qiime.beta_metrics


def get_nonphylogenetic_metric(name):
    """Gets metric by name from distance_transform.

    Metrics should be f(matrix) -> distances.
    """
    # looks for name, inserting possible dist_ to find functions
    # in distance_transform.py named e.g.:
    # binary_dist_chisq / dist_bray_curtis
    try:
        return getattr(distance_transform, 'dist_' + name.lower())
    except AttributeError:
        try:
            return getattr(distance_transform,
                           name.replace('binary', 'binary_dist').lower())
        except AttributeError:
            return getattr(distance_transform,
                           name.lower())


def get_phylogenetic_metric(name):
    """Gets metric by name from cogent.maths.unifrac.fast_unifrac

    Metrics should be f(matrix) -> distances.
    """
    # looks for name, inserting possible dist_ to find functions
    # in qiime.beta_metrics
    try:
        return getattr(qiime.beta_metrics, 'dist_' + name.lower())
    except AttributeError:
        try:
            return getattr(qiime.beta_metrics,
                           name.replace('binary', 'binary_dist').lower())
        except AttributeError:
            return getattr(qiime.beta_metrics,
                           name.lower())


def get_phylogenetic_row_metric(name):
    """Gets metric by name from qiime.beta_metrics

    Metrics should be f(matrix) -> distance_array.
    """
    # looks for name, inserting one_sample to find functions
    # in qiime.beta_metrics
    return getattr(qiime.beta_metrics, 'one_sample_' + name.lower())


def list_known_nonphylogenetic_metrics():
    """Lists known metrics by name from distance_transform.

    Assumes that functions starting with dist_ or binary_dist are metrics.
    """
    result = []
    for name in dir(distance_transform):
        if name.startswith('dist_'):
            result.append(name[5:])
        elif name.startswith('binary_dist_'):
            result.append('binary_' + name[12:])
    result.sort()
    return result


def list_known_phylogenetic_metrics():
    """Lists known phylogenetic metrics from cogent.maths.unifrac."""
    result = []
    for name in dir(qiime.beta_metrics):
        if name.startswith('dist_'):
            result.append(name[5:])
    result.sort()
    return result


def list_known_metrics():
    return list_known_nonphylogenetic_metrics() + \
        list_known_phylogenetic_metrics()


def single_file_beta(input_path, metrics, tree_path, output_dir,
                     rowids=None, full_tree=False):
    """ does beta diversity calc on a single otu table

    uses name in metrics to name output beta diversity files
    assumes input tree is already trimmed to contain only otus present in otu
    table, doesn't call getSubTree()
    inputs:
     input_path (str)
     metrics (str, comma delimited if more than 1 metric; or list)
     tree_path (str)
     output_dir (str)
     rowids (comma separated str)
    """
    metrics_list = metrics
    try:
        metrics_list = metrics_list.split(',')
    except AttributeError:
        pass

    otu_table = load_table(input_path)

    otumtx = asarray([v for v in otu_table.iter_data(axis='sample')])

    if tree_path:
        tree = parse_newick(open(tree_path, 'U'),
                            PhyloNode)
    else:
        tree = None

    input_dir, input_filename = os.path.split(input_path)
    input_basename, input_ext = os.path.splitext(input_filename)
    for metric in metrics_list:
        outfilepath = os.path.join(output_dir, metric + '_' +
                                   input_basename + '.txt')
        try:
            metric_f = get_nonphylogenetic_metric(metric)
            is_phylogenetic = False
        except AttributeError:
            try:
                metric_f = get_phylogenetic_metric(metric)
                is_phylogenetic = True
                if tree is None:
                    stderr.write("metric %s requires a tree, but none found\n"
                                 % (metric,))
                    exit(1)
            except AttributeError:
                stderr.write("Could not find metric %s.\n\nKnown metrics are: %s\n"
                             % (metric, ', '.join(list_known_metrics())))
                exit(1)
        if rowids is None:
            # standard, full way
            if is_phylogenetic:
                dissims = metric_f(otumtx, otu_table.ids(axis='observation'),
                                   tree, otu_table.ids(),
                                   make_subtree=(not full_tree))
            else:
                dissims = metric_f(otumtx)
            f = open(outfilepath, 'w')
            f.write(format_distance_matrix(otu_table.ids(), dissims))
            f.close()
        else:
            # only calc d(rowid1, *) for each rowid
            rowids_list = rowids.split(',')
            row_dissims = []  # same order as rowids_list
            for rowid in rowids_list:
                rowidx = otu_table.index(rowid, axis='sample')

                # first test if we can the dissim is a fn of only the pair
                # if not, just calc the whole matrix
                if metric_f.__name__ == 'dist_chisq' or \
                        metric_f.__name__ == 'dist_gower' or \
                        metric_f.__name__ == 'dist_hellinger' or\
                        metric_f.__name__ == 'binary_dist_chisq':
                    warnings.warn('dissimilarity ' + metric_f.__name__ +
                                  ' is not parallelized, calculating the whole matrix...')
                    row_dissims.append(metric_f(otumtx)[rowidx])
                else:
                    try:
                        row_metric = get_phylogenetic_row_metric(metric)
                    except AttributeError:
                        # do element by element
                        dissims = []
                        sample_ids = otu_table.ids()
                        observation_ids = otu_table.ids(axis='observation')
                        for i in range(len(sample_ids)):
                            if is_phylogenetic:
                                dissim = metric_f(
                                    otumtx[[rowidx, i], :], observation_ids,
                                    tree, [sample_ids[rowidx], sample_ids[i]],
                                    make_subtree=(not full_tree))[0, 1]
                            else:
                                dissim = metric_f(otumtx[[rowidx, i], :])[0, 1]
                            dissims.append(dissim)
                        row_dissims.append(dissims)
                    else:
                        # do whole row at once
                        dissims = row_metric(otumtx,
                                             otu_table.ids(axis='observation'),
                                             tree, otu_table.ids(), rowid,
                                             make_subtree=(not full_tree))
                        row_dissims.append(dissims)

            with open(outfilepath, 'w') as f:
                f.write(format_matrix(row_dissims, rowids_list,
                                      otu_table.ids(),
                                      convert_matching_names_to_zero=True))


def multiple_file_beta(input_path, output_dir, metrics, tree_path,
                       rowids=None, full_tree=False):
    """ runs beta diversity for each input file in the input directory

    performs minimal error checking on input args, then calls single_file_beta
    for each file in the input directory

    inputs:
     input_path (str)
     metrics (str, comma delimited if more than 1 metric; or list)
     tree_path (str)
     output_dir (str)

    """
    metrics_list = metrics
    try:
        metrics_list = metrics_list.split(',')
    except AttributeError:
        pass

    file_names = os.listdir(input_path)
    file_names = [fname for fname in file_names if not fname.startswith('.')]
    try:
        os.makedirs(output_dir)
    except OSError:
        pass  # hopefully dir exists
    for metric in metrics_list:
        try:
            metric_f = get_nonphylogenetic_metric(metric)
        except AttributeError:
            try:
                metric_f = get_phylogenetic_metric(metric)
                if not tree_path:
                    raise ValueError("a tree file is required for " + metric)
            except AttributeError:
                raise ValueError(
                    "Could not find metric %s.\n\nKnown metrics are: %s\n"
                    % (metric, ', '.join(list_known_metrics())))

    for fname in file_names:
        single_file_beta(os.path.join(input_path, fname),
                         metrics, tree_path, output_dir, rowids, full_tree)
