#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Dan Knight"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Dan Knights"]
__license__ = "GPL"
__version__ = "1.1.0-dev"
__maintainer__ = "Dan Knights"
__email__ = "daniel.knights@colorado.edu"
__status__ = "Development"
 

from optparse import make_option
from os import makedirs, rmdir, listdir
from os.path import join
from cogent.util.misc import remove_files
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.parse import parse_mapping_file
from qiime.supervised_learning import run_R_supervised_learner,\
    R_format_map_file, R_format_otu_table
options_lookup = get_options_lookup()
valid_methods = ['random_forest']

script_info={}
script_info['brief_description']="""Run supervised classification using \
OTUs as predictors and a mapping file category as class labels."""
script_info['script_description']="""This script trains a supervised classifier using OTUs \
(or other continuous input sample x observation data) as predictors, and a \
mapping file column containing discrete values as the class labels.

    Outputs:
        predictions.txt: the labels predicted by the classifier for the given
            samples. Each sample is predicted by a model that was trained 
            without it. 
        probabilities.txt: the label probabilities for each of the given 
            samples. (if available)
        summary.txt: a summary of the results, including the expected
            generalization error of the classifier
        features.txt: a list of discriminative OTUs with their associated
            importance scores (if available)
        params.txt: a list of any non-default parameters used in training
            the model.
    
By default, this script removes OTUs that are present in less than 10% of the \
samples. Run with "--show_params" to see how to set this manually. For an \
overview of the application of supervised classification to microbiota, \
see PubMed ID 21039646.

This script requires that R is installed and in the search path. To install R \
visit: http://www.r-project.org/. Once R is installed, run R and excecute the \
command "install.packages('randomForest')", then type q() to exit."""

script_info['script_usage']=[]
script_info['script_usage'].append(("""Simple example of random forests classifier""","""""","""supervised_learning.py -i otutable.txt -m map.txt -c 'Individual' -o ml"""))
script_info['script_usage'].append(("""Getting a sample params file for the random forests classifier""","""""","""supervised_learning.py -i otutable.txt -m map.txt -c 'Individual' -o ml --show_params"""))
script_info['script_usage'].append(("""Running with a user-specified params file for the random forests classifier""","""""","""supervised_learning.py -i otutable.txt -m map.txt -c 'Individual' -o ml -p params.txt"""))
script_info['output_description']="""Outputs a ranking of features (e.g. OTUs) by importance, an estimation of the generalization error of the classifier, and the predicted class labels and posterior class probabilities \
according to the classifier."""
script_info['required_options'] = [\
    make_option('-i', '--input_data', help='Input data file containing predictors (e.g. otu table)'),
    make_option('-m', '--mapping_file', help='File containing meta data (response variables)'),
    make_option('-c', '--category', help='Name of meta data category to predict'),
    make_option('-o','--output_dir',\
            help='the output directory [REQUIRED]'),\
]
script_info['optional_options']=[\
    make_option('-s', '--method', default='random_forest',
        help= 'Comma-separated list of supervised learning methods to apply. '\
            + 'Currently one option is available: "random_forest" '\
            + '[default: %default].'),

    make_option('-f','--force',action='store_true',\
        dest='force',help='Force overwrite of existing output directory'+\
        ' (note: existing files in output_dir will not be removed)'+\
        ' [default: %default]'),\
    make_option('-p','--param_file',type='string',\
        help='file containing parameters for the ' +\
             'supervised learning model inference [default: %default]',
             default=None),\
    make_option('--show_params',action="store_true",\
        help='show sample parameters file for a given method [default: %default]'),\
    make_option('-k','--keepfiles',action='store_true',\
        help='Keep R-formatted input files [default: %default]')\
]
script_info['version'] = __version__

sample_params = {
    'random_forest':\
"""
# example params file for random forests classifier
# commented lines (starting with "#") are ignored

# number of trees in forest; more improves generalization error estimate
params$ntree = 500

# seed integer for the random number generator (between 0 and 1e9)
# can be used to replicate results for stochastic processes
params$seed = 0

# specify the minimum number of samples in which an OTU must be present
# if both min.num.samples and min.percent.samples are specified, the larger
# of the two is used
params$min.num.samples = 50

# specify the minimum percent of samples in which an OTU must be present
params$min.percent.samples = .1

"""
}

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # create the output directories
    try:
        data_output_dir = join(opts.output_dir, "data")
        makedirs(data_output_dir)
    except OSError:
        if opts.force:
            pass
        else:
            # This check helps users avoid overwriting previous output.
            print "Output directory already exists. Please choose "+\
             "a different directory, or force overwrite with -f."
            exit(1)

    # make sub-dir for each learning method
    model_names = opts.method.split(',')
    for method in model_names:
        # ensure method is valid
        if not method in valid_methods:
            print "Method", method, "not found in valid methods (" +\
                ','.join(valid_methods) + ")."
            exit(1)
        # print model information if requested
        if opts.show_params:
            print sample_params[method]
            exit(0)
        subdir = join(opts.output_dir, method)
        try:
            makedirs(subdir)
        except OSError:
            pass

    # verify that category is in mapping file
    map_list = parse_mapping_file(open(opts.mapping_file,'U').readlines())
    if not opts.category in map_list[1][1:]:
        print "Category '%s' not found in mapping file columns:" %(opts.category)
        print map_list[1][1:]
        exit(1)

    # convert otu table and mapping file to R format
    otu_fp = R_format_otu_table(opts.input_data, data_output_dir)
    map_fp = R_format_map_file(opts.mapping_file, data_output_dir)

    # run the supervised learning algorithm
    results = run_R_supervised_learner(otu_fp, map_fp,
        opts.category, model_names, opts.output_dir,
        param_file=opts.param_file)
    
    # delete R-formatted otu table and map file
    if not opts.keepfiles:
        remove_files([otu_fp, map_fp])
        # attempt to remove data dir only if it is empty
        if len(listdir(data_output_dir)) == 0:
            rmdir(data_output_dir)
    
if __name__ == "__main__":
    main()
