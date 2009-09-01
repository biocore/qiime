#!/usr/bin/env python

__author__ = "Justin Kuzynski"
__copyright__ = "Copyright 2009, the PyCogent Project"
__credits__ = ["Justin Kuczynski", "Rob Knight"]
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Prototype"

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
from qiime.util import FunctionWithParams, TreeMissingError, OtuMissingError
from qiime.parse import parse_otus
from qiime.format import format_distance_matrix
from optparse import OptionParser
import cogent.maths.distance_transform #avoid hard-coding metrics
import qiime.beta_metrics
import qiime.beta_diversity
from sys import exit, stderr
import os.path

def get_nonphylogenetic_metric(name):
    """Gets metric by name from cogent.maths.distance_transform.
    
    Metrics should be f(matrix) -> distances.
    """
    try:
        return getattr(cogent.maths.distance_transform, name.lower())
    except AttributeError:
        return getattr(cogent.maths.distance_transform, 'dist_'+name.lower())

def get_phylogenetic_metric(name):
    """Gets metric by name from cogent.maths.unifrac.fast_unifrac
    
    Metrics should be f(matrix) -> distances.
    """
    return getattr(qiime.beta_metrics, name.lower())

def list_known_nonphylogenetic_metrics():
    """Lists known metrics by name from cogent.maths.distance_transform.

    Assumes that functions starting with dist_ or binary_dist are metrics.
    """
    result = []
    for name in dir(cogent.maths.distance_transform):
        if name.startswith('dist_') or name.startswith('binary_dist'):
            result.append(name)
    result.sort()
    return result

def list_known_phylogenetic_metrics():
    """Lists known phylogenetic metrics from cogent.maths.unifrac."""
    result = []
    for name in dir(qiime.beta_metrics):
        if name.startswith('dist_'):
            result.append(name)
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
        self.Metric = metric    #should be f(table, tree) -> dist matrix
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
        #if it's a phylogenetic metric, read the tree
        if self.IsPhylogenetic:
            tree = self.getTree(tree_path)
        else:
            tree = None
        
        sample_names, taxon_names, data, lineages = \
            self.getOtuTable(data_path)
        #prune tree so we only look at otus actually present
        if tree:
            tree = tree.getSubTree(taxon_names, ignore_missing=True)

        # get the 2d dist matrix from beta diversity analysis
        if self.IsPhylogenetic:
            return self.Metric(data.T, taxon_names, tree), sample_names
        else:
            return self.Metric(data.T), sample_names

    def formatResult(self, result):
        """Generate formatted distance matrix - result is (data, sample_names)"""
        data, sample_names = result
        return format_distance_matrix(sample_names, data)



def make_cmd_parser():
    """returns command-line options"""
    usage = """ %prog [options]
use %prog -h for help.

example: 
python beta_diversity.py -i TEST/otu_table.txt -m dist_unweighted_unifrac -o TEST/beta_unweighted_unifrac.txt -t TEST/repr_set.tre
or batch example: 
python beta_diversity.py -i TEST/beta_rare/ -m dist_unweighted_unifrac -t TEST/repr_set.tre -o TEST/rare_unifrac
processes every file in beta_rare, and creates a file "beta_" + fname
in results folder
-o is mandatory here


use -s to see metric options.
Output will be a sample by sample distance matrix. 
Default is to write to stdout for non-batch file processing."""
    
    parser = OptionParser(usage=usage)
    parser.add_option('-s', '--show_metrics', action='store_true', 
        dest="show_metrics",
        help='show available beta diversity metrics and quit')
    parser.add_option('-t', '--tree', dest='tree_path', default=None,
        help='path to newick tree file, required for phylogenetic metrics')
    parser.add_option('-o', '--out_path', dest='output_path', default=None,
        help='output path')
    parser.add_option('-i', '--in_path', dest='input_path', default=None,
        help='input path')
    parser.add_option('-m', '--metric', dest='metric', default=None,
        help='metric to use')    

    options, args = parser.parse_args()
    if (len(args) != 0 and options.show_metrics != True) :
        parser.error("incorrect number of arguments")
    return options, args


def single_file_beta(options, args):
    metric = options.metric
    try:
        metric_f = get_nonphylogenetic_metric(metric)
        is_phylogenetic = False
    except AttributeError:
        try:
            metric_f = get_phylogenetic_metric(metric)
            is_phylogenetic = True
        except AttributeError:
            stderr.write("Could not find metric %s.\n\nKnown metrics are: %s\n" \
                % (metric, ', '.join(list_known_metrics())))
            exit(1)

    c = BetaDiversityCalc(metric_f, metric, is_phylogenetic)

    try:
        result = c(data_path=options.input_path, tree_path=options.tree_path, 
            result_path=options.output_path, log_path=None)
        if result: #can send to stdout instead of file
            print c.formatResult(result)
    except IOError, e:
        stderr.write("Failed because of missing files.\n")
        stderr.write(str(e)+'\n')
        exit(1)

def multiple_file_beta(options, args):
    beta_script = qiime.beta_diversity.__file__
    file_names = os.listdir(options.input_path)
    for fname in file_names:
        beta_div_cmd = 'python ' + beta_script + ' -i '+\
            os.path.join(options.input_path, fname) + " -m " + options.metric\
            + ' -o ' + os.path.join(options.output_path,'beta_'+fname)
        if options.tree_path:
            beta_div_cmd += ' -t ' + options.tree_path
        os.system(beta_div_cmd)


if __name__ == '__main__':
    options, args = make_cmd_parser()
    if options.show_metrics:
        print("Known metrics are: %s\n" \
                % (', '.join(list_known_metrics()),))
        exit(0)

    if os.path.isdir(options.input_path):
        multiple_file_beta(options, args)
    elif os.path.isfile(options.input_path):
        single_file_beta(options, args)
    else:
        print("io error")
        exit(1)
