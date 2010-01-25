#!/usr/bin/env python

__author__ = "Justin Kuzynski"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Rob Knight"]
__license__ = "GPL"
__status__ = "1.0-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Pre-release"

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
import sys
import os.path

def get_nonphylogenetic_metric(name):
    """Gets metric by name from cogent.maths.distance_transform.
    
    Metrics should be f(matrix) -> distances.
    """
    # looks for name, inserting possible dist_ to find functions
    # in distance_transform.py named e.g.:
    # binary_dist_chisq / dist_bray_curtis
    try:
        return getattr(cogent.maths.distance_transform, 'dist_' + name.lower())
    except AttributeError:
        try:
            return getattr(cogent.maths.distance_transform, 
              name.replace('binary', 'binary_dist').lower())
        except AttributeError:
            return getattr(cogent.maths.distance_transform, 
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

def list_known_nonphylogenetic_metrics():
    """Lists known metrics by name from cogent.maths.distance_transform.

    Assumes that functions starting with dist_ or binary_dist are metrics.
    """
    result = []
    for name in dir(cogent.maths.distance_transform):
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
        """Generate formatted distance matrix. result is (data, sample_names)"""
        data, sample_names = result
        return format_distance_matrix(sample_names, data)

def single_file_beta(options, args):
    """ does beta diversity calc on a single otu table

    uses name in options.metrics to name output beta diversity files"""
    metrics_list = options.metrics.split(',')
    for metric in metrics_list:
        outfilepath = os.path.join(options.output_dir, metric + '_' +
            os.path.split(options.input_path)[1])        
        try:
            metric_f = get_nonphylogenetic_metric(metric)
            is_phylogenetic = False
        except AttributeError:
            try:
                metric_f = get_phylogenetic_metric(metric)
                is_phylogenetic = True
            except AttributeError:
                stderr.write("Could not find metric %s.\n\nKnown metrics are: %s\n"\
                    % (metric, ', '.join(list_known_metrics())))
                exit(1)
        calc = BetaDiversityCalc(metric_f, metric, is_phylogenetic)

        try:
            result = calc(data_path=options.input_path, 
                tree_path=options.tree_path, 
                result_path=outfilepath, log_path=None)
            if result: #can send to stdout instead of file
                print c.formatResult(result)
        except IOError, e:
            stderr.write("Failed because of missing files.\n")
            stderr.write(str(e)+'\n')
            exit(1)

def multiple_file_beta(options, args):
    """ runs beta diversity for each input file in the input directory 
    
    performs minimal error checking on input args, then calls os.system
    to execute single_file_beta for each file in the input directory

    this is to facilitate future task farming - replace os.system with 
    write to file, each command is independant

    """
    metrics_list = options.metrics.split(',')
    beta_script = qiime.beta_diversity.__file__
    file_names = os.listdir(options.input_path)
    file_names = [fname for fname in file_names if not fname.startswith('.')]
    try:
        os.makedirs(options.output_dir)
    except OSError:
        pass # hopefully dir exists
    for metric in metrics_list:
        try:
            metric_f = get_nonphylogenetic_metric(metric)
        except AttributeError:
            try:
                metric_f = get_phylogenetic_metric(metric)
                if not options.tree_path:
                    raise ValueError("a tree file is required for " + metric)
            except AttributeError:
                raise ValueError(
                    "Could not find metric %s.\n\nKnown metrics are: %s\n" \
                    % (metric, ', '.join(list_known_metrics())))
    
    for fname in file_names:
            outfilepath = options.output_dir
            beta_div_cmd = 'python ' + beta_script + ' -i '+\
                os.path.join(options.input_path, fname) + " -m " + options.metrics\
                + ' -o ' + outfilepath
            if options.tree_path:
                beta_div_cmd += ' -t ' + options.tree_path
            os.system(beta_div_cmd)

usage_str = """ %prog [options] {-i INPUT_PATH -o OUTPUT_DIR -m METRICS} or {-s}

[] indicates optional input (order unimportant)
{} indicates required input (order unimportant)

Example:

python %prog -i otu_table.txt -m bray_curtis,unweighted_unifrac -o outdir -t repr_set.tre

this creates two files: outdir/bray_curtis_otu_table.txt etc.

or batch example: 
python %prog -i beta_rare_dir -m bray_curtis,unweighted_unifrac -o outdir -t repr_set.tre
processes every file in beta_rare_dir, and creates a file "metric_" + infilename
in results folder

use -s to see metric options.  (prepending dist_ to most names works as well)
Output will be a sample by sample distance matrix. 
"""
def parse_command_line_parameters():
    """returns command-line options"""

    if len(sys.argv) == 1:
        sys.argv.append('-h')
    usage = usage_str
    version = '%prog ' + str(__version__)
    parser = OptionParser(usage=usage, version=version)

    parser.add_option('-i', '--input_path',
        help='input path [REQUIRED]')
        
    parser.add_option('-o', '--output_dir',
        help='output directory [REQUIRED]')

    parser.add_option('-m', '--metrics',
        help='metrics to use, comma delimited if >1 metric, '+\
            'no spaces [REQUIRED]')  
        
    parser.add_option('-s', '--show_metrics', action='store_true', 
        dest="show_metrics",
        help='show available beta diversity metrics and quit')

    parser.add_option('-t', '--tree_path', default=None,
        help='path to newick tree file, required for phylogenetic metrics'+\
        ' [default: %default]')  

    opts, args = parser.parse_args()
        
    if len(args) != 0:
        parser.error("positional argument detected.  make sure all"+\
         ' parameters are identified.' +\
         '\ne.g.: include the \"-m\" in \"-m MINIMUM_LENGTH\"')
         
    required_options = ['input_path','output_dir']
    if not opts.show_metrics:
        for option in required_options:
            if eval('opts.%s' % option) == None:
                parser.error('Required option --%s omitted.' % option)

    return opts, args
    
if __name__ == '__main__':
    options, args = parse_command_line_parameters()
    if options.show_metrics:
        print("Known metrics are: %s\n" \
                % (', '.join(list_known_metrics()),))
        exit(0)
    if options.output_dir.endswith('.txt'):
        stderr.write('output must be a directory, files will be named'+\
            ' automatically.  And we refuse to make .txt directories\n')
        exit(1)
    try: 
        os.makedirs(options.output_dir)
    except OSError:
        pass # hopefully dir already exists 
    if os.path.isdir(options.input_path):
        multiple_file_beta(options, args)
    elif os.path.isfile(options.input_path):
        single_file_beta(options, args)
    else:
        print("io error, input path not valid.  Does it exist?")
        exit(1)
