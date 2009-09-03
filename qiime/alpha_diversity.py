#!/usr/bin/env python

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2009, the PyCogent Project" #consider project name
__credits__ = ["Justin Kuczynski", "Rob Knight"] #remember to add yourself
__license__ = "GPL"
__version__ = "0.1"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"
__status__ = "Prototype"

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

#note: might want to use make_safe_f to strip out additional params passed on.
import numpy
import os.path
from optparse import OptionParser

import cogent.maths.stats.alpha_diversity as alph

from qiime.parse import parse_otus, make_envs_dict
from qiime.util import FunctionWithParams, make_safe_f
from qiime.format import format_matrix
from sys import exit, stderr
import os.path
import qiime.alpha_diversity

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
        me
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
        data = self.getData(data_path)
        if self.IsPhylogenetic:
            tree = self.getTree(tree_path)
            envs = make_envs_dict(data, sample_names, taxon_names)
            new_taxon_names, result = self.Metric(tree, envs, **self.Params)
            return numpy.array(result)
        else:
            def metric(row):
                return self.Metric(row, **self.Params)
            result = map(metric, data)
            
            return numpy.array(result)

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
        
        sample_names, taxon_names, data, lineages = \
            self.getOtuTable(data_path)

        calc_names = []
        for calc in self.Calcs:
            # add either calc's multiple return value names, or fn name
            calc_names.extend(getattr(calc.Metric, 'return_names', \
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
            metric_res = c(data_path=data.T, taxon_names=taxon_names, 
                tree_path=tree, sample_names=sample_names)
            if len(metric_res.shape) == 1:
                res.append(metric_res)
            elif len(metric_res.shape) == 2:
                for met in metric_res.T:
                    res.append(met) 
            else:
                raise RuntimeError, "alpha div shape not as expected"
        res_data = numpy.array(res).T
        
        return res_data, sample_names, calc_names
    
    def formatResult(self, result):
        """Generate formatted distance - result is (data, sample_names)"""
        data, sample_names, calc_names = result
        res = format_matrix(data, sample_names, calc_names)
        return res


def make_cmd_parser():
    """returns command-line options"""
    usage = """ %prog [options]
use %prog -h for help.

example: 
python alpha_diversity.py -i TEST/otu_table.txt -m osd,PD_whole_tree -o TEST/alpha_osd_PD.txt -t TEST/repr_set.tre

or batch example: 
python alpha_diversity.py -i TEST/alpha_rare -m chao1,PD_whole_tree -t TEST/repr_set.tre -o TEST/rare_chao1_PD
processes every file in alpha_rare, and creates a file "alpha_" + fname
in results folder
-o is mandatory here

Metrics is comma delimited, use -s to see options.
Output will be a sample by metric matrix. 
Default is to write to stdout for non-batch file processing"""
    
    parser = OptionParser(usage=usage)
    parser.add_option('-s', '--show_metrics', action='store_true', 
        dest="show_metrics",
        help='show available alpha diversity metrics and quit')
    parser.add_option('-t', '--tree', dest='tree_path', default=None,
        help='path to newick tree file, required for phylogenetic metrics')
    parser.add_option('-o', '--out_path', dest='output_path', default=None,
        help='output path')
    parser.add_option('-i', '--in_path', dest='input_path', default=None,
        help='output path')
    parser.add_option('-m', '--metrics', dest='metrics', default=None,
        help='output path')    

    options, args = parser.parse_args()
    if (len(args) != 0 and options.show_metrics != True) :
        parser.error("incorrect number of arguments")
    return options, args

def get_nonphylogenetic_metric(name):
    """Gets metric by name from list in this module
    """
    for metric in nonphylogenetic_metrics:
        if metric.__name__.lower() == name.lower(): return metric
    
    raise AttributeError

def get_phylogenetic_metric(name):
    """Gets metric by name from list in this module
    """
    for metric in phylogenetic_metrics:
        if metric.__name__.lower() == name.lower(): return metric
    
    raise AttributeError

def list_known_metrics():
    nonphylo = [metric.__name__ for metric in nonphylogenetic_metrics]
    phylo = [metric.__name__ for metric in phylogenetic_metrics]
    return nonphylo + phylo
    

import cogent.maths.unifrac.fast_unifrac as fast_unifrac
# change name of function if returns multiple values:
alph.chao1_confidence.return_names = ('chao1_lower_bound', 'chao1_upper_bound')
alph.osd.return_names = ('observed', 'singles', 'doubles')

#hand curated lists of metrics, these either return one value, or
#are modified above
phylogenetic_metrics = [fast_unifrac.PD_whole_tree]
nonphylogenetic_metrics = [
alph.berger_parker_d,
alph.brillouin_d,
alph.chao1,
alph.chao1_confidence,
alph.dominance,
alph.doubles,
alph.equitability,
alph.fisher_alpha,
alph.heip_e,
alph.kempton_taylor_q,
alph.margalef,
alph.mcintosh_d,
alph.mcintosh_e,
alph.menhinick,
alph.michaelis_menten_fit,
alph.observed_species,
alph.osd,
alph.reciprocal_simpson,
alph.robbins,
alph.shannon,
alph.simpson,
alph.simpson_e,
alph.singles,
alph.strong]

def single_file_alpha(infilepath, metrics, outfilepath, tree_path):
    metrics_list = metrics.split(',')
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
                stderr.write("Could not find metric %s.\n Known metrics are: %s\n" \
                    % (metric, ', '.join(list_known_metrics())))
                exit(1)
        c = AlphaDiversityCalc(metric_f, is_phylogenetic)
        calcs.append(c)
    
    all_calcs = AlphaDiversityCalcs(calcs)

    try:
        result = all_calcs(data_path=infilepath, tree_path=tree_path,
            result_path=outfilepath, log_path=None)
        if result:  #can send to stdout instead of file
            print all_calcs.formatResult(result)
    except IOError, e:
        stderr.write("Failed because of missing files.\n")
        stderr.write(str(e)+'\n')
        exit(1)


def multiple_file_alpha(options, args):
    """ performs minimal error checking on input args, then calls os.system
    to execute single_file_alpha for each file in the input directory

    this is to facilitate future task farming - replace os.system with 
    write to file, each command is independant

    """
    alpha_script = qiime.alpha_diversity.__file__
    file_names = os.listdir(options.input_path)
    if not os.path.exists(options.output_path):
        raise ValueError("output path does not appear to exist")

    for metric in options.metrics:
        try:
            metric_f = get_nonphylogenetic_metric(metric)
        except AttributeError:
            try:
                metric_f = get_phylogenetic_metric(metric)
            except AttributeError:
                raise ValueError("could not find metric.  %s.\n Known metrics are: %s\n" \
                    % (metric, ', '.join(list_known_metrics())))
    for fname in file_names:
        # future: try to make sure fname is a valid otu file
        #~ f = open(os.path.join(dir,fname))
        #~ try:
            #~ parse_otus(f)
        #~ except:
            #~ print fname + " not a valid otu table" # write this to log file?
            #~ continue
        #~ finally:
            #~ f.close()
            #~ print 'finally', fname
#~ 
        #~ print 'running w/ ' + fname
        alpha_div_cmd = 'python ' + alpha_script + ' -i '+\
            os.path.join(options.input_path, fname) + " -m " + options.metrics\
            + ' -o ' + os.path.join(options.output_path,'alpha_'+fname)
        if options.tree_path:
            alpha_div_cmd += ' -t ' + options.tree_path
        os.system(alpha_div_cmd)

if __name__ == '__main__':
    options, args = make_cmd_parser()
    if options.show_metrics:
        print("Known metrics are: %s\n" \
                % (', '.join(list_known_metrics()),))
        exit(0)

    if os.path.isdir(options.input_path):
        multiple_file_alpha(options, args)
    elif os.path.isfile(options.input_path):
        single_file_alpha(options.input_path, options.metrics, 
            options.output_path, options.tree_path)
    else:
        print("io error, input path not valid. does it exist?")
        exit(1)
