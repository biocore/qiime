#!/usr/bin/env python

__author__ = "Justin Kuzynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Rob Knight", "Greg Caporaso",
               "Jose Carlos Clemente Litran", "Jai Ram Rideout",
               "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
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
from biom.parse import parse_biom_table
from biom.table import Table
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


class BetaDiversityCalc(FunctionWithParams):

    """A BetaDiversityCalc takes taxon x sample counts, returns distance mat.

    Some BetaDiversityCalcs also require a tree.

    typical usage is:
    calc = BetaDiversityCalc(my_metric, 'metric name', \
        IsPhylogenetic=False)
    calc(data_path=stuff, result_path=stuff2)
    """
    _tracked_properties = ['IsPhylogenetic', 'Metric']

    def __init__(self, metric, name, is_phylogenetic=False, params=None):
        """Return new BetaDiversityCalc object with specified params.

        Note: metric must conform to the following interface:
            must return 2d sample distance matrix with sample order preserved.
            * phylogenetic methods: take samples by taxa array, taxon_names,
            cogent PhyloNode tree
            * nonphylo methods: take only samples by taxa array, numpy 2d

        Not sure what params does
        """
        self.Metric = metric  # should be f(table, tree) -> dist matrix
        self.Name = name
        self.IsPhylogenetic = is_phylogenetic
        self.Params = params or {}

    def getResult(self, data_path, tree_path=None):
        """Returns distance matrix from (indcidence matrix and optionally tree).

        Parameters:

        data_path: path to data file, matrix (samples = cols, taxa = rows)
        in tab-delimited text format

        tree_path: path or object.
        if method is phylogenetic, must supply tree_path.
        if path, path to
        Newick-format tree file where taxon ids match taxon ids in the
        input data file.

        returns 2d dist matrix, list of sample names ordered as in dist mtx
        """
        # if it's a phylogenetic metric, read the tree
        if self.IsPhylogenetic:
            tree = self.getTree(tree_path)
        else:
            tree = None

        otu_table = parse_biom_table(open(data_path, 'U'))
        otumtx = asarray([v for v in otu_table.iter_sample_data()])

        # get the 2d dist matrix from beta diversity analysis
        if self.IsPhylogenetic:
            return (self.Metric(otumtx, otu_table.observation_ids, tree,
                                otu_table.sample_ids),
                    list(otu_table.sample_ids))
        else:
            return self.Metric(otumtx), list(otu_table.sample_ids)

    def formatResult(self, result):
        """Generate formatted distance matrix. result is (data, sample_names)"""
        data, sample_names = result
        return format_distance_matrix(sample_names, data)


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

    otu_table = parse_biom_table(open(input_path, 'U'))

    otumtx = asarray([v for v in otu_table.iter_sample_data()])

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
                dissims = metric_f(otumtx, otu_table.observation_ids,
                                   tree, otu_table.sample_ids, make_subtree=(not full_tree))
            else:
                dissims = metric_f(otumtx)
            f = open(outfilepath, 'w')
            f.write(format_distance_matrix(otu_table.sample_ids, dissims))
            f.close()
        else:
            # only calc d(rowid1, *) for each rowid
            rowids_list = rowids.split(',')
            row_dissims = []  # same order as rowids_list
            for rowid in rowids_list:
                rowidx = otu_table.sample_ids.index(rowid)

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
                        for i in range(len(otu_table.sample_ids)):
                            if is_phylogenetic:
                                dissim = metric_f(otumtx[[rowidx, i], :],
                                                  otu_table.observation_ids, tree,
                                                  [otu_table.sample_ids[rowidx],
                                                   otu_table.sample_ids[i]],
                                                  make_subtree=(not full_tree))[0, 1]
                            else:
                                dissim = metric_f(otumtx[[rowidx, i], :])[0, 1]
                            dissims.append(dissim)
                        row_dissims.append(dissims)
                    else:
                        # do whole row at once
                        dissims = row_metric(otumtx,
                                             otu_table.observation_ids, tree,
                                             otu_table.sample_ids, rowid,
                                             make_subtree=(not full_tree))
                        row_dissims.append(dissims)

            # rows_outfilepath = os.path.join(output_dir, metric + '_' +\
            #     '_'.join(rowids_list) + '_' + os.path.split(input_path)[1])
            f = open(outfilepath, 'w')
            f.write(
                format_matrix(
                    row_dissims,
                    rowids_list,
                    otu_table.sample_ids))
            f.close()


def single_object_beta(otu_table, metrics, tr, rowids=None,
                       full_tree=False):
    """mod of single_file_beta to recieve and return otu obj, tree str

    uses name in metrics to name output beta diversity files
    assumes input tree is already trimmed to contain only otus present
    in otu_table, doesn't call getSubTree()
    inputs:
                otu_table -- a otu_table in the biom format
                metrics -- metrics (str, comma delimited if more than 1 metric)
                tr -- a phylonode cogent tree object if needed by the chosen beta
                                        diversity metric
                rowids -- comma seperated string
    """
    otumtx = asarray([v for v in otu_table.iter_sample_data()])

    if tr:
        tree = tr
    else:
        tree = None

    metrics_list = metrics.split(',')

    for metric in metrics_list:
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
                dissims = metric_f(otumtx, otu_table.observation_ids, tree,
                                   otu_table.sample_ids, make_subtree=(not full_tree))
            else:
                dissims = metric_f(otumtx)

            return (
                format_distance_matrix(
                    otu_table.sample_ids,
                    dissims).split(
                    '\n')
            )
        else:
            # only calc d(rowid1, *) for each rowid
            rowids_list = rowids.split(',')
            row_dissims = []  # same order as rowids_list
            for rowid in rowids_list:
                rowidx = otu_table.sample_ids.index(rowid)

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
                        for i in range(len(otu_table.sample_ids)):
                            if is_phylogenetic:
                                dissim = metric_f(otumtx[[rowidx, i], :],
                                                  otu_table.observation_ids, tree,
                                                  [otu_table.sample_ids[rowidx],
                                                   otu_table.sample_ids[i]],
                                                  make_subtree=(not full_tree))[0, 1]
                            else:
                                dissim = metric_f(otumtx[[rowidx, i], :])[0, 1]
                            dissims.append(dissim)
                        row_dissims.append(dissims)
                    else:
                        # do whole row at once
                        dissims = row_metric(otumtx,
                                             otu_table.observation_ids, tree,
                                             otu_table.sample_ids, rowid,
                                             make_subtree=(not full_tree))
                        row_dissims.append(dissims)

            return format_matrix(row_dissims, rowids_list, otu_table.sample_ids)


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
