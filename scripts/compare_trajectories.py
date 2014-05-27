#!/usr/bin/env python
from __future__ import division

__author__ = "Jose Antonio Navas Molina"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jose Antonio Navas Molina", "Antonio Gonzalez Pena",
               "Yoshiki Vazquez Baeza"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Jose Antonio Navas Molina"
__email__ = "josenavasmolina@gmail.com"

from os import mkdir
from os.path import join, exists

import pandas as pd

from skbio.maths.stats.ordination import OrdinationResults

from qiime.util import parse_command_line_parameters, make_option
from qiime.trajectory_analysis import (run_trajectory_analysis,
                                       TRAJECTORY_ALGORITHMS)


script_info = {}
script_info['brief_description'] = """"""
script_info['script_description'] = """"""
script_info['script_usage'] = [
    ("Example", "", "%prog")
]
script_info['required_options'] = [
    make_option('-i', '--input_ordination_fp', type='existing_filepath',
                help="Input ordination results filepath"),
    make_option('-m', '--map_fp', type='existing_filepath',
                help="Input metadata mapping filepath"),
    make_option('-c', '--category', type='string',
                help="Category name of the mapping file to use to create the "
                     "vectors"),
    make_option('-o', '--output_dir', type='new_dirpath',
                help="Name of the output directory to save the results")
]
script_info['optional_options'] = [
    make_option('-s', '--sort_by', type='string', default=None,
                help="Category name of the mapping file to use to sort"),
    make_option('--algorithm', type='choice', default='avg',
                choices=TRAJECTORY_ALGORITHMS,
                help="The algorithm used to create the vectors. The method "
                     "used can be RMS (either using 'avg' or 'trajectory'); or"
                     " the first difference (using 'diff'), or 'wdiff' for a "
                     "modified first difference algorithm (see --window_size) "
                     "the aforementioned use all the dimensions and weights "
                     "them using their percentage explained; returns the norm "
                     "of the created vectors; and their confidence using ANOVA"
                     ". The Vectors are created as follows: for 'avg' it "
                     "calculates the average at each timepoint (averaging "
                     "within a group), then calculates the norm of each point;"
                     " for 'trajectory' calculates the norm for the 1st-2nd, "
                     "2nd-3rd, etc.; for 'diff', it calculates the norm for "
                     "all the time-points and then calculates the first "
                     "difference for each resulting point; for 'wdiff' it uses"
                     " the same procedure as the previous method but the "
                     "subtraction will be between the mean of the next number "
                     "of elements specified in --window_size and the current "
                     "element, both methods ('wdiff' and 'diff') will also "
                     "include the mean and the standard deviation of the "
                     "calculations [defautl: %default]"),
    make_option('--vectors_axes', dest='vectors_axes', type='int', default=3,
                help="The number of axes to account while doing the vector "
                     "specific calculations. We suggest using 3 because those "
                     "are the ones being displayed in the plots but you could "
                     "use any number between 1 and number of samples- 1. To "
                     "use all of them pass 0. [default: %default]"),
    make_option('-w', '--weight_by_vector', default=False, action='store_true',
                help="Use -w when you want the output to be weighted by the "
                     "space between samples in the --sort_by column, i. e. "
                     "days between samples [default: %default]"),
    make_option('--window_size', type='int', default=None,
                help="Use --window_size, when selecting the modified first "
                     "difference (\'wdiff\') option for --algorithm. This "
                     "integer determines the number of elements to be averaged"
                     " per element subtraction, the resulting vector. "
                     "[default: %default]")
]
script_info['version'] = __version__

if __name__ == '__main__':
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    ord_fp = opts.input_ordination_fp
    mapping_fp = opts.map_fp
    vector_category = opts.category
    output_dir = opts.output_dir
    sort_by = opts.sort_by
    algorithm = opts.algorithm
    vectors_axes = opts.vectors_axes
    weighted = opts.weight_by_vector
    window_size = opts.window_size

    # Parse the ordination results
    with open(ord_fp, 'U') as f:
        ord_res = OrdinationResults.from_file(f)

    # Parse the mapping file
    # Using dropna to remove comments - with how='all' removes a row
    # if all its column values are NA, which should occur only with in comments
    metamap = pd.read_csv(mapping_fp, sep='\t', index_col=0).dropna(how='all')

    if vector_category not in metamap.keys():
        option_parser.error("Vector category %s does not exist in the mapping "
                            "file" % vector_category)

    if sort_by:
        if sort_by == 'SampleID':
            sort_category = None
        elif sort_by not in metamap.keys():
            option_parser.error("Sort category %s does not exist in the "
                                "mapping file" % sort_by)
        else:
            sort_category = sort_by

    if vectors_axes < 0 or vectors_axes > len(ord_res.eigvals):
        option_parser.error("--vectors_axes should be between 0 and the max "
                            "number of axes available (%d), found: %d "
                            % (len(ord_res.eigvals), vectors_axes))
    elif vectors_axes == 0:
        # If vector_axes = 0, we generate all of them
        vectors_axes = len(ord_res.eigvals)

    if weighted:
        if sort_by == 'SampleID':
            option_parser.error("The weighting_vector can't be the SampleID, "
                                "please specify a numeric column in the "
                                "--sort_by option.")
        elif not sort_category:
            option_parser.error("You should provide --sort_by if you want to "
                                "weight the output")

    if algorithm == 'wdiff':
        if not window_size:
            option_parser.error("You should provide --window_size when using "
                                "the wdiff algorithm")
        if window_size < 1:
            option_parser.error("--window_size should be a positive integer, "
                                "%d found" % window_size)

    res = run_trajectory_analysis(ord_res, metamap, vector_category,
                                  sort_category, algorithm, vectors_axes,
                                  weighted, window_size)

    if not exists(output_dir):
        mkdir(output_dir)

    vectors_fp = join(output_dir, 'vectors.txt')
    vectors_raw_fp = join(output_dir, 'vectors_raw_values.txt')
    with open(vectors_fp, 'w') as f:
        f.write('\n'.join(res))
    with open(vectors_raw_fp, 'w') as f:
        f.write('\n'.join(res))