#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Rob Knight", "Greg Caporaso",
               "William Van Treuren", "Jose Antonio Navas Molina",
               "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

"""Contains code for performing alpha diversity analyses within each sample.

This module has the responsibility for calculating diversity within each of a
set of samples.
designed to be called from the command line, but see also unit tests
for how to interface directly with the classes
see python alpha_diversity.py -h for more help

primary class is AlphaDiversityCalcs


log files aren't yet supported

fails uglily when called from an empty otu table file, make sure your files
aren't just headers
e.g.:
#../tests/scratch/rare_11_0.txt
#OTU ID [nothing else here]
tax1    [or here]

* param passing to metrics isn't supported
"""

import os.path
import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')

import skbio.diversity.alpha as alph
from numpy import array, zeros

from qiime.util import FunctionWithParams
from qiime.format import format_matrix
from sys import exit, stderr


class AlphaDiversityCalc(FunctionWithParams):

    """An AlphaDiversityCalc takes taxon x sample counts, returns diversities.

    Some AlphaDiversityCalcs also require a tree.

    Note: sample ids should be preserved during this process.
    """
    _tracked_properties = ['IsPhylogenetic', 'Metric']

    def __init__(self, metric, is_phylogenetic=False, params=None):
        """Return new AlphaDiversityCalc object with specified params.

        Note: expect params to contain both generic and per-method (e.g. for
        PD vs. Michaelis-Menten vs. whatever) params, so leaving it as a dict
        rather than setting attributes.

        Parameters:
        IMPORTANT:
        metric: typically function such that f(count_vector) -> float.
        can return multiple values, if so, f.return_names must be set to
        a tuple of strings, len = len of return tuple
        is_phylogenetic: set to True if metric is phylogenetic (and therefore
            needs tree file)

        Entries in params will be passed on directly to functions, may need
        to adapt to allow this to work using make_safe_f. **not yet


        ### this wasn't done, in favor of using metric.return_names
        Note: some methods produce confidence intervals, others don't.
        Handle by passing in tuple instead of string to names: expect same
        number of return rows as there are values in self.Name.
        Incidentally, this solution also allows you to make compound
        AlphaDiversityCalc objects that combine multiple calculations, might
        be worth naming these if there are particular combinations that are
        routinely useful.
        """
        self.Metric = metric
        self.Name = metric.__name__
        self.IsPhylogenetic = is_phylogenetic
        self.Params = params or {}

    def getResult(self, data_path, taxon_names=None, sample_names=None,
                  tree_path=None):
        """Returns per-sample diversity from incidence matrix and optional tree.

        Parameters:

        data_path: can be either a file path or an array,
        if array: either numpy array or list of numpy arrays where each row is a
        sample, contents are counts of each taxon, must be dense to allow
        phylogenetic calcs (where the taxon you have matters).
        must be 2d.  for one sample just do [sample_array]

        taxon_names: list of names of taxa, same order as in row (required for
        phylogenetic methods)

        tree: cogent.tree.PhyloNode object, or file path

        output:
        1d/2d array containing diversity of each sample, preserving order from
        input data  sample by (metric name or metric.return_name)
        1d: [(metric on sample1), (metric on sample2),...
        2d: [(return val 1 from sample1),(return val 2)...]
            [(return val 1 on sample2),...]
        """
        otu_table = self.getBiomData(data_path)
        data = otu_table.iter_data(axis='sample')
        if self.IsPhylogenetic:
            tree = self.getTree(tree_path)
            # build envs dict: envs = {otu_id:{sample_id:count}}
            envs = {}
            sample_ids = otu_table.ids()
            for obs_v, obs_id, obs_md in otu_table.iter(axis='observation'):
                obs = {}
                for sample_id, v in zip(sample_ids, obs_v):
                    obs[sample_id] = v
                envs[obs_id] = obs

            new_sample_names, result = self.Metric(tree, envs, **self.Params)
            ordered_res = zeros(len(sample_names), 'float')
            for i, sample in enumerate(sample_names):
                try:
                    # idx is sample's index in result from metric
                    idx = new_sample_names.index(sample)
                    ordered_res[i] = result[idx]
                except ValueError:
                    pass  # already is zero
            return array(ordered_res)

        else:
            def metric(row):
                return self.Metric(row.astype(int), **self.Params)
            result = map(metric, data)

            return array(result)

    def formatResult(self, result):
        """Generate formatted vector, here just tab-delimited text.
        WARNING: does this work on metrics with multiple return values?
        command line doesn't use this"""
        return '\t'.join(map(str, result))


class AlphaDiversityCalcs(FunctionWithParams):

    """Coordinates list of AlphaDiversityCalc objects to apply to same data.
    also used by command line if only one metric is being run on data
    """

    def __init__(self, alpha_calcs, params=None):
        """Initialize with list of AlphaDiversityCalc objects."""
        self.Params = params or {}
        self.Calcs = alpha_calcs

    def getResult(self, data_path, tree_path=None):
        """computes matrix of (samples vs alpha diversity methods) on data

        Specificially, will produce matrix where each col is one of the alpha
        diversity calculations loaded on __init__ and each row is a
        sample -- will also return/write row and col headers.

        inputs:
        * data_path: file path, tab delimited, otu table format --OR --
        tuple: (sample_names, taxon_names, data (2d numpy), lineages)
        * tree: newick tree path --OR-- cogent.core.tree.PhyloNode object

        output:
        result: a matrix of sample by alpha diversity method, sample_names,
        calc_names
        """
        otu_table = self.getBiomData(data_path)

        calc_names = []
        for calc in self.Calcs:
            # add either calc's multiple return value names, or fn name
            calc_names.extend(getattr(calc.Metric, 'return_names',
                                      (calc.Metric.__name__,)))
        needs_tree = max([c.IsPhylogenetic for c in self.Calcs])
        if needs_tree:
            tree = self.getTree(tree_path)
        else:
            tree = None
        # calculations
        res = []
        for c in self.Calcs:
            # add either calc's multiple return value names, or fn name
            metric_res = c(data_path=data_path,
                           taxon_names=otu_table.ids(axis='observation'),
                           tree_path=tree,
                           sample_names=otu_table.ids())
            if len(metric_res.shape) == 1:
                res.append(metric_res)
            elif len(metric_res.shape) == 2:
                for met in metric_res.T:
                    res.append(met)
            else:
                raise RuntimeError("alpha div shape not as expected")
        res_data = array(res).T

        return res_data, otu_table.ids(), calc_names

    def formatResult(self, result):
        """Generate formatted distance - result is (data, sample_names)"""
        data, sample_names, calc_names = result
        res = format_matrix(data, sample_names, calc_names)
        return res


def get_nonphylogenetic_metric(name):
    """Gets metric by name from list in this module
    """
    for metric in nonphylogenetic_metrics:
        if metric.__name__.lower() == name.lower():
            return metric

    raise AttributeError


def get_phylogenetic_metric(name):
    """Gets metric by name from list in this module
    """
    for metric in phylogenetic_metrics:
        if metric.__name__.lower() == name.lower():
            return metric

    raise AttributeError


def list_known_metrics():
    nonphylo = [metric.__name__ for metric in nonphylogenetic_metrics]
    phylo = [metric.__name__ for metric in phylogenetic_metrics]
    return nonphylo + phylo


import cogent.maths.unifrac.fast_unifrac as fast_unifrac
# change name of function if returns multiple values:
alph.chao1_ci.return_names = ('chao1_lower_bound', 'chao1_upper_bound')
alph.osd.return_names = ('observed', 'singles', 'doubles')
alph.esty_ci.return_names = ('esty_lower_bound', 'esty_upper_bound')
alph.lladser_ci.return_names = ('lladser_lower_bound', 'lladser_upper_bound')

# hand curated lists of metrics, these either return one value, or
# are modified above
phylogenetic_metrics = [fast_unifrac.PD_whole_tree]

# maintain additional aliases for backwards compatibility
def observed_species(counts):
    return alph.observed_otus(counts)


def simpson_reciprocal(counts):
    return alph.enspie(counts)

nonphylogenetic_metrics = [
    alph.ace,
    alph.berger_parker_d,
    alph.brillouin_d,
    alph.chao1,
    alph.chao1_ci,
    alph.dominance,
    alph.doubles,
    alph.enspie,
    alph.equitability,
    alph.esty_ci,
    alph.fisher_alpha,
    alph.gini_index,
    alph.goods_coverage,
    alph.heip_e,
    alph.kempton_taylor_q,
    alph.margalef,
    alph.mcintosh_d,
    alph.mcintosh_e,
    alph.menhinick,
    alph.michaelis_menten_fit,
    alph.observed_otus,
    observed_species,
    alph.osd,
    simpson_reciprocal,
    alph.robbins,
    alph.shannon,
    alph.simpson,
    alph.simpson_e,
    alph.singles,
    alph.strong]

cup_metrics = [alph.lladser_pe, alph.lladser_ci]


def single_file_alpha(infilepath, metrics, outfilepath, tree_path):
    metrics_list = metrics
    try:
        metrics_list = metrics_list.split(',')
    except AttributeError:
        pass

    calcs = []
    for metric in metrics_list:
        try:
            metric_f = get_nonphylogenetic_metric(metric)
            is_phylogenetic = False
        except AttributeError:
            try:
                metric_f = get_phylogenetic_metric(metric)
                is_phylogenetic = True
            except AttributeError:
                stderr.write(
                    "Could not find metric %s.\n Known metrics are: %s\n"
                    % (metric, ', '.join(list_known_metrics())))
                exit(1)
        c = AlphaDiversityCalc(metric_f, is_phylogenetic)
        calcs.append(c)

    all_calcs = AlphaDiversityCalcs(calcs)

    try:
        result = all_calcs(data_path=infilepath, tree_path=tree_path,
                           result_path=outfilepath, log_path=None)
        if result:  # can send to stdout instead of file
            print all_calcs.formatResult(result)
    except IOError as e:
        stderr.write("Failed because of missing files.\n")
        stderr.write(str(e) + '\n')
        exit(1)


def multiple_file_alpha(input_path, output_path, metrics, tree_path=None):
    """ performs minimal error checking on input args, then calls os.system
    to execute single_file_alpha for each file in the input directory

    this is to facilitate future task farming - replace os.system with
    write to file, each command is independant

    """
    file_names = os.listdir(input_path)
    file_names = [fname for fname in file_names if not fname.startswith('.')]
    if not os.path.exists(output_path):
        os.makedirs(output_path)

    metrics_list = metrics
    try:
        metrics_list = metrics_list.split(',')
    except AttributeError:
        pass

    for metric in metrics_list:
        try:
            metric_f = get_nonphylogenetic_metric(metric)
        except AttributeError:
            try:
                metric_f = get_phylogenetic_metric(metric)
                # bail if we got a phylo metric but no tree file
                if tree_path is None:
                    raise ValueError("phylogenetic metric supplied, but no " +
                                     "phylogenetic tree supplied")
                elif not os.path.exists(tree_path):
                    raise ValueError("phylogenetic metric supplied, but no " +
                                     "phylogenetic tree found in specified location")
            except AttributeError:
                raise ValueError(
                    "could not find metric.  %s.\n Known metrics are: %s\n"
                    % (metric, ', '.join(list_known_cup_metrics())))

    for fname in file_names:
        # future: try to make sure fname is a valid otu file

        output_fname = 'alpha_' + os.path.splitext(fname)[0] + '.txt'
        output_fp = os.path.join(output_path, output_fname)

        single_file_alpha(os.path.join(input_path, fname), metrics_list,
                          output_fp, tree_path)


def single_file_cup(otu_filepath, metrics, outfilepath, r):
    """Compute variations of the conditional uncovered probability.

    otu_filepath: path to otu_table file
    metrics: comma separated list of required metrics; or list
    outfilepath: path to output file
    r: Number of new colors that are required for the next prediction

    The opposite of uncovered probability is sometimes called coverage.
    """
    calcs = []

    params = {'r': r}

    metrics_list = metrics
    try:
        metrics_list = metrics_list.split(',')
    except AttributeError:
        pass

    for metric in metrics_list:
        try:
            metric_f = get_cup_metric(metric)
        except AttributeError:
            stderr.write(
                "Could not find metric %s.\n Known metrics are: %s\n"
                % (metric, ', '.join(list_known_cup_metrics())))
            exit(1)

        c = AlphaDiversityCalc(metric_f, params=params)
        calcs.append(c)

    all_calcs = AlphaDiversityCalcs(calcs)

    try:
        result = all_calcs(data_path=otu_filepath,
                           result_path=outfilepath, log_path=None)
        if result:  # can send to stdout instead of file
            print all_calcs.formatResult(result)
    except IOError as e:
        stderr.write("Failed because of missing files.\n")
        stderr.write(str(e) + '\n')
        exit(1)


def get_cup_metric(name):
    """Gets metric by name from list in this module
    """
    for metric in cup_metrics:
        if metric.__name__.lower() == name.lower():
            return metric
    raise AttributeError


def list_known_cup_metrics():
    """Show the names of available metrics."""
    return [metric.__name__ for metric in cup_metrics]
