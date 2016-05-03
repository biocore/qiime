#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Dan Knights"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Dan Knights", "Luke Ursell"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Dan Knights"
__email__ = "daniel.knights@colorado.edu"


from qiime.util import make_option
from os import makedirs, listdir
from os.path import join, isdir
from glob import glob
from numpy import mean
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.parse import parse_mapping_file
from qiime.supervised_learning import (
    run_supervised_learning, assemble_results,
    calc_baseline_error_to_observed_error, pooled_standard_deviation)
options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = """Run supervised classification using \
OTUs as predictors and a mapping file category as class labels."""
script_info['script_description'] = """This script trains a supervised classifier using OTUs \
(or other continuous input sample x observation data) as predictors, and a \
mapping file column containing discrete values as the class labels.

Outputs:
    * cv_probabilities.txt: the label probabilities for each of the \
        given samples. (if available)
    * mislabeling.txt: A convenient presentation of cv_probabilities \
        for mislabeling detection.
    * confusion_matrix.txt: confusion matrix for hold-out predictions.
    * summary.txt: a summary of the results, including the expected \
        generalization error of the classifier
    * feature_importance_scores.txt: a list of discriminative OTUs with their \
        associated importance scores (if available)

It is recommended that you remove low-depth samples and rare OTUs \
before running this script. This can drastically reduce the run-time, and in \
many circumstances will not hurt performance. It is also recommended to perform \
rarefaction to control for sampling effort before running this \
script. For example, to rarefy at depth 200, then remove OTUs present in \
< 10 samples run:

single_rarefaction.py -i otu_table.biom -d 200 -o otu_table_rarefied200.biom
filter_otus_from_otu_table.py -i otu_table_rarefied200.biom -s 10 -o otu_table_rarefied200.present10.biom

For an overview of the application of supervised classification to microbiota, \
see PubMed ID 21039646.

This script also has the ability to collate the supervised learning results \
produced on an input directory. For example, in order to reduce any variation \
introduced through producing a rarefied OTU table, the user can run \
multiple_rarefactions_even_depth.py on the OTU table, and then pass that directory \
into supervised_learning.py. The user can then pass a -w collate_results filepath \
to produce a single results file that contains the average estimated generalization \
error of the classified, and the pooled standard deviation (for cv5 and cv10 errortypes).

This script requires that R be installed and in the search path. To install R \
visit: http://www.r-project.org/. Once R is installed, run R and excecute the \
command "install.packages('randomForest')", then type q() to exit."""

script_info['script_usage'] = []

script_info['script_usage'].append(
    ("""Simple example of random forests classifier""",
     """""",
     """%prog -i otu_table.biom -m Fasting_Map.txt -c BarcodeSequence -o ml"""))

script_info['script_usage'].append(
    ("""Running with 10-fold cross-validation for improved estimates of generalization error and feature importances""",
     """""",
     """%prog -i otu_table.biom -m Fasting_Map.txt -c BarcodeSequence -o ml_cv10 -e cv10"""))

script_info['script_usage'].append(
    ("""Running with 1,000 trees for improved generalization error""",
     """""",
     """%prog -i otu_table.biom -m Fasting_Map.txt -c BarcodeSequence -o ml_ntree1000 --ntree 1000"""))

script_info['script_usage'].append(
    ("Run 10-fold cross validation on a directory of OTU tables rarefied at an even depth""",
     """""",
     """%prog -i rarefied_tables/ -m Fasting_Map.txt -c Treatment -o sl_rarefied_tables_cv10 -e cv10"""))

script_info['script_usage'].append(
    ("Run 10-fold cross validation on a directory of OTU tables rarefied at an even depth and collate the results into a single file""",
     """""",
     """%prog -i rarefied_tables/ -m Fasting_Map.txt -c Treatment -o sl_rarefied_tables_cv10_sweep -e cv10 -w sl_cv10_sweep.txt"""))

script_info[
    'script_usage_output_to_remove'] = [
    'ml',
    'ml_cv10',
    'ml_ntree1000',
    'sl_rarefied_tables_cv10',
    'sl_rarefied_tables_cv10_sweep']

# this example is better suited for the tutorial as it's going to be difficult to use in
# automated testing
# script_info['script_usage'].append(("""Simple example, filter OTU table first""","""""","""
#  single_rarefaction.py -i otu_table_filtered.txt -d 200 -o otu_table_rarefied200.txt
#  filter_otus_from_otu_table.py -i otu_table_rarefied200.txt -s 10
# supervised_learning.py -i otutable_filtered_rarefied200.txt -m map.txt
# -c 'Individual' -o ml"""))


script_info['output_description'] = """Outputs a ranking of features (e.g. OTUs) by importance, an estimation of the generalization error of the classifier, and the predicted class labels and posterior class probabilities \
according to the classifier."""
script_info['required_options'] = [
    make_option('-i', '--input_data', type='existing_path',
                help='Input data file containing predictors (e.g. otu table) '
                'or a directory of otu tables'),
    make_option('-m', '--mapping_file', type='existing_filepath',
                help='File containing meta data (response variables)'),
    make_option('-c', '--category', type='string', help='Name of meta data '
                'category to predict'),
    make_option('-o', '--output_dir', type='new_dirpath',
                help='the output directory'),
]

errortype_choices = ['oob', 'loo', 'cv5', 'cv10']

script_info['optional_options'] = [
    make_option('-f', '--force', action='store_true',
                dest='force', help='Force overwrite of existing output directory' +
                ' (note: existing files in output_dir will not be removed)' +
                ' [default: %default]'),
    make_option('--ntree', type='int', default=500,
                help='Number of trees in forest (more is better but slower) '
                '[default: %default]'),
    make_option('-e', '--errortype', type='choice', default='oob',
                choices=errortype_choices,
                help='type of error estimation. Valid choices are: ' +
                ', '.join(errortype_choices) + '. ' +
                'oob: out-of-bag, fastest, only builds one classifier, use for '
                'quick estimates; cv5: 5-fold cross validation, provides mean and '
                'standard deviation of error, use for good estimates on very '
                'large data sets; cv10: 10-fold cross validation, provides mean and '
                'standard deviation of error, use for best estimates; '
                'loo: leave-one-out cross validation, use for small data sets '
                '(less than ~30-50 samples) [default %default]'),
    make_option(
        '-w', '--collate_results_fp', default=None, type='new_filepath',
        help='When passing in a directory of OTU tables that are rarefied '
        'at an even depth, this option will collate the results into a single '
        'specified output file, averaging the estimated errors and standard deviations. '
        '[default: %default]')
]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    input_data = opts.input_data
    mapping_file = opts.mapping_file
    category = opts.category
    ntree = opts.ntree
    errortype = opts.errortype
    output_dir = opts.output_dir
    verbose = opts.verbose
    force = opts.force
    collate_results_fp = opts.collate_results_fp

    # create the output directories
    try:
        makedirs(opts.output_dir)
    except OSError:
        if force:
            pass
        else:
            # This check helps users avoid overwriting previous output.
            option_parser.error("Output directory already exists. Please choose"
                                " a different directory, or force overwrite with -f.")

    # verify that category is in mapping file
    map_list = parse_mapping_file(open(mapping_file, 'U').readlines())
    if not category in map_list[1][1:]:
        option_parser.error(
            "Category '%s' not found in mapping file columns:" %
            (category))
        print map_list[1][1:]
        exit(1)

    # if input is a single otu table
    if isdir(input_data) is False:

        # run the supervised learning algorithm
        result = run_supervised_learning(input_data, mapping_file, category,
                                         ntree, errortype, output_dir, verbose)

    # if input is a directory of otu tables
    if isdir(input_data) is True:
        input_tables = glob('%s/*biom' % input_data)

        coll_est_error = []
        coll_est_error_stdev = []
        baseline_error = []

        for table_fp in input_tables:
            # create output dir on per-table basis with convention:
            # "sl_TABLENAME_CATEGORY/"
            output_basename = table_fp.split('/')[-1]
            output_basename = output_basename.replace('.biom', '')
            output_name = "sl_%s_%s/" % (output_basename, category)
            output_fp = join(output_dir, output_name)
                # create the output directories
            try:
                makedirs(output_fp)
            except OSError:
                if force:
                    pass
                else:
                    # This check helps users avoid overwriting previous output.
                    option_parser.error("Output directory already exists. Please choose"
                                        " a different directory, or force overwrite with -f.")

            result = run_supervised_learning(table_fp, mapping_file, category,
                                             ntree, errortype, output_fp, verbose)

            # retrieve the estimated error and baseline error
            est_error_line, baseline_error_line = \
                result['summary'].readlines()[2:4]

            est_error_line = est_error_line.split('\t')[1]
            coll_est_error.append(float(est_error_line.split(' ')[0]))

            # only collect standard deviations for cv5 and cv10 errortypes
            if errortype in ['cv5', 'cv10']:
                est_error_stdev = est_error_line.split(' ')[2].strip()
                coll_est_error_stdev.append(float(est_error_stdev))

            # make sure baseline error is the same across all tables (it should
            # be)
            if baseline_error == []:
                baseline_error.append(
                    float(baseline_error_line.split('\t')[1].strip()))

        if collate_results_fp:
            output_file = open(collate_results_fp, 'w')

            # get assembled results
            results = assemble_results(coll_est_error, coll_est_error_stdev,
                                       baseline_error[0], errortype, ntree)
            output_file.write('\n'.join(results))
            output_file.close()

if __name__ == "__main__":
    main()
