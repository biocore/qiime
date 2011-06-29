#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Dan Knights"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Dan Knights"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Dan Knights"
__email__ = "daniel.knights@colorado.edu"
__status__ = "Release"
 

from qiime.util import make_option
from os import makedirs, rmdir, listdir
from os.path import join
from cogent.util.misc import remove_files
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.parse import parse_mapping_file
from qiime.supervised_learning import RSupervisedLearner,\
    RSupervisedLearnerFilter, R_format_table
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
    
It is strongly recommended that you remove low-depth samples and rare OTUs \
before running this script. This can drastically reduce the run-time, and in \
many circumstances will not hurt performance. It is also recommended to perform \
rarefaction to control for sampling effort before running this \
script. For example, to rarefy at depth 200, then remove OTUs present in \
< 10 samples run:

single_rarefaction.py -i otu_table_filtered.txt -d 200 -o otu_table_rarefied200.txt
filter_otu_table.py -i otu_table_rarefied200.txt -s 10

Run this script with "--show_params" to see how to set any model-specific parameters. \
For an overview of the application of supervised classification to microbiota, \
see PubMed ID 21039646.

This script requires that R is installed and in the search path. To install R \
visit: http://www.r-project.org/. Once R is installed, run R and excecute the \
command "install.packages("randomForest")", then type q() to exit."""

script_info['script_usage']=[]
script_info['script_usage'].append(("""Simple example of random forests classifier""","""""","""supervised_learning.py -i otutable.txt -m map.txt -c 'Individual' -o ml"""))
script_info['script_usage'].append(("""Simple example, filter OTU table first""","""""",\
"""single_rarefaction.py -i otu_table_filtered.txt -d 200 -o otu_table_rarefied200.txt
 filter_otu_table.py -i otu_table_rarefied200.txt -s 10
 supervised_learning.py -i otutable_filtered_rarefied200.txt -m map.txt -c 'Individual' -o ml
"""))

script_info['script_usage'].append(("""Getting a sample params file for the random forests classifier""","""""","""supervised_learning.py -i otutable.txt -m map.txt -c 'Individual' -o ml --show_params"""))
script_info['script_usage'].append(("""Running with a user-specified params file for the random forests classifier""","""""","""supervised_learning.py -i otutable.txt -m map.txt -c 'Individual' -o ml -p params.txt"""))
script_info['output_description']="""Outputs a ranking of features (e.g. OTUs) by importance, an estimation of the generalization error of the classifier, and the predicted class labels and posterior class probabilities \
according to the classifier."""
script_info['required_options'] = [\
    make_option('-i', '--input_data', help='Input data file containing predictors (e.g. otu table)'),
    make_option('-m', '--mapping_file', help='File containing meta data (response variables)'),
    make_option('-c', '--category', help='Name of meta data category to predict'),
]
script_info['optional_options']=[\
    make_option('-o','--output_dir',default='.',\
            help='the output directory [deafult: %default]'),
    make_option('-s', '--method', default='random_forest',
        help= 'Comma-separated list of supervised learning methods to apply. '\
            + 'Currently one option is available: "random_forest" '\
            + '[default: %default].'),
    make_option('-f','--force',action='store_true',\
        dest='force',help='Force overwrite of existing output directory'+\
        ' (note: existing files in output_dir will not be removed)'+\
        ' [default: %default]'),
    make_option('-p','--param_file',type='string',\
        help='file containing parameters for the ' +\
             'supervised learning model inference [default: %default]',
             default=None),
    make_option('--show_params',action="store_true",default=False,\
        help='show sample parameters file for a given method [default: %default]'),
    make_option('--filter_type',type="string",default=None,\
        help='type of filter to use. Currently one is available: BSSWSS. [default: %default]'),
    make_option('--filter_min',type='int',default=2,\
        help='minimum number of features to try with filter [default: %default]'),
    make_option('--filter_max',type='int',default=20,\
        help='maximum number of features to try with filter [default: %default]'),
    make_option('--filter_step',type='int',default=1,\
        help='step increment for number of features to try with filter [default: %default]'),
    make_option('--filter_reps',type='int',default=10,\
        help='Number of models to train for estimating filter error [default: %default]'),
    make_option('-k','--keepfiles',action='store_true',\
        help='Keep R-formatted input files [default: %default]'),
]
script_info['version'] = __version__

sample_params = {
    'random_forest':\
"""
# example params file for random forests classifier
# commented lines (starting with "#") are ignored

# number of trees in forest; more improves generalization error estimate
params$ntree = 1000

# seed integer for the random number generator (between 0 and 1e9)
# can be used to replicate results for stochastic processes
params$seed = 0
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
    if opts.verbose:
        print 'Formatting otu table...'
    otu_fp = join(data_output_dir, 'otus_R_format.txt')
    R_format_table(opts.input_data, output_filepath=otu_fp)
    if opts.verbose:
        print 'Formatting mapping file...'
    map_fp = join(data_output_dir, 'map_R_format.txt')
    R_format_table(opts.mapping_file, output_filepath=map_fp)

    # run the supervised learning algorithm
    if opts.verbose:
        print 'Running R...'

    if not opts.filter_type is None:
        learner = RSupervisedLearnerFilter()
    else:
        learner = RSupervisedLearner()
    results = learner(otu_fp, map_fp, opts.category, 
                        model_names, opts.output_dir,
                        param_file=opts.param_file, filter=opts.filter_type,
                        filter_min=opts.filter_min,filter_max=opts.filter_max,
                        filter_step=opts.filter_step, filter_reps=opts.filter_reps,
                        verbose=opts.verbose)
    
    # delete R-formatted otu table and map file
    if not opts.keepfiles:
        remove_files([otu_fp, map_fp])
        # attempt to remove data dir only if it is empty
        if len(listdir(data_output_dir)) == 0:
            rmdir(data_output_dir)
    
if __name__ == "__main__":
    main()
