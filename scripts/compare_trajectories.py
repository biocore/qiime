#!/usr/bin/env python
from __future__ import division

__author__ = "Jose Antonio Navas Molina"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jose Antonio Navas Molina", "Antonio Gonzalez Pena",
               "Yoshiki Vazquez Baeza"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jose Antonio Navas Molina"
__email__ = "josenavasmolina@gmail.com"

from os import mkdir
from os.path import join, exists

import pandas as pd
from skbio.stats.ordination import OrdinationResults

from qiime.parse import parse_mapping_file_to_dict
from qiime.util import parse_command_line_parameters, make_option
from qiime.compare_trajectories import (run_trajectory_analysis,
                                        TRAJECTORY_ALGORITHMS)


script_info = {}
script_info['brief_description'] = ("Run analysis of volatility using a "
                                    "variety of algorithms")
script_info['script_description'] = (
    "This script mainly allows performing analysis of volatility on time "
    "series data, but they can be applied to any data that contains a gradient"
    ". The methods available are RMS (either using 'avg' or 'trajectory'); or "
    "the first difference (using 'diff'), or 'wdiff' for a modified first "
    "difference algorithm. The trajectories are computed as follows. For 'avg'"
    " it calculates the average point within a group and then computes the "
    "norm of the distance of each sample from the averaged value. For "
    "'trajectory' each component of the result trajectory is computed as "
    "taking the sorted list of samples in the group and taking the norm of the"
    " coordinates of the 2nd samples minus the 1st sample, 3rd sample minus "
    "2nd sample and so on. For 'diff' it calculates the norm for all the "
    "time-points and then calculates the first difference for each resulting "
    "point. For 'wdiff', it calculates the norm for all the time-points and "
    "substracts the mean of the next number of elements, specified using the "
    "'--window_size' parameters, and the current element.")
script_info['script_usage'] = [
    ("Average method",
     "Execute the analysis of volatility using the average method, grouping "
     "the samples using the 'Treatment' category",
     "%prog -i pcoa_res.txt -m map.txt -c 'Treatment' -o avg_output"),
    ("Trajectory method",
     "Execute the analysis of volatility using the trajectory method, grouping"
     " the samples using the 'Treatment' category and sorting them using the "
     "'time' category",
     "%prog -i pcoa_res.txt -m map.txt -c 'Treatment' --algorithm trajectory "
     "-o trajectory_output -s time"),
    ("First difference method",
     "Execute the analysis of volatility using the first difference method, "
     "grouping the samples using the 'Treatment' category, sorting them using "
     "the 'time' category and calculating the trajectory using the first "
     "four axes",
     "%prog -i pcoa_res.txt -m map.txt -c 'Treatment' --algorithm diff "
     "-o diff_output -s time --axes 4"),
    ("Window difference method",
     "Execute the analysis of volatility using the window difference method, "
     "grouping the samples using the 'Treatment' category, sorting them using "
     "the 'time' category, weighting the output by the space between "
     "samples in the 'time' category and using a window size of three.",
     "%prog -i pcoa_res.txt -m map.txt -c 'Treatment' --algorithm wdiff "
     "-o wdiff_output -s time --window_size 3 -w")
]
script_info['required_options'] = [
    make_option('-i', '--input_fp', type='existing_filepath',
                help="Input ordination results filepath"),
    make_option('-m', '--map_fp', type='existing_filepath',
                help="Input metadata mapping filepath"),
    make_option('-c', '--categories', type='str',
                help="Comma-separated list of category names of the mapping "
                     "file to use to create the trajectories"),
    make_option('-o', '--output_dir', type='new_dirpath',
                help="Name of the output directory to save the results")
]
script_info['optional_options'] = [
    make_option('-s', '--sort_by', type='str', default=None,
                help="Category name of the mapping file to use to sort"),
    make_option('--algorithm', type='choice', default='avg',
                choices=TRAJECTORY_ALGORITHMS,
                help="The algorithm to use. Available methods: "
                     + str(TRAJECTORY_ALGORITHMS) + ". [Default: %default]"),
    make_option('--axes', type='int', default=3,
                help="The number of axes to account while doing the trajectory"
                     " specific calculations. We suggest using 3 because those"
                     " are the ones being displayed in the plots but you could"
                     " use any number between 1 and number of samples - 1. To "
                     "use all of them pass 0. [default: %default]"),
    make_option('-w', '--weight_by_vector', default=False, action='store_true',
                help="Use -w when you want the output to be weighted by the "
                     "space between samples in the --sort_by column, i. e. "
                     "days between samples [default: %default]"),
    make_option('--window_size', type='int', default=None,
                help="Use --window_size, when selecting the modified first "
                     "difference ('wdiff') option for --algorithm. This "
                     "integer determines the number of elements to be averaged"
                     " per element subtraction, the resulting trajectory. "
                     "[default: %default]")
]
script_info['output_description'] = (
    "This script generates two files in the output directory, "
    "'trajectories.txt' and 'trajectories_raw_values.txt'. The "
    "'trajectories.txt' file includes the resulting statistics and a list of "
    "categories that did not passed the tests to run the analysis. The "
    "'trajectories_raw_values.txt' file includes the resulting trajectory for "
    "each used category.")
script_info['version'] = __version__

if __name__ == '__main__':
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    ord_fp = opts.input_fp
    mapping_fp = opts.map_fp
    categories = opts.categories.split(',')
    output_dir = opts.output_dir
    sort_by = opts.sort_by
    algorithm = opts.algorithm
    axes = opts.axes
    weighted = opts.weight_by_vector
    window_size = opts.window_size

    # Parse the ordination results
    with open(ord_fp, 'U') as f:
        ord_res = OrdinationResults.read(f)

    # Parse the mapping file
    with open(mapping_fp, 'U') as f:
        map_dict = parse_mapping_file_to_dict(f)[0]
    metamap = pd.DataFrame.from_dict(map_dict, orient='index')

    for category in categories:
        if category not in metamap.keys():
            option_parser.error("Category %s does not exist in the mapping "
                                "file" % categories)

    sort_category = None
    if sort_by:
        if sort_by == 'SampleID':
            sort_category = None
        elif sort_by not in metamap.keys():
            option_parser.error("Sort category %s does not exist in the "
                                "mapping file. Available categories are: %s" %
                                (sort_by, metamap.keys()))
        else:
            sort_category = sort_by

    if axes < 0 or axes > len(ord_res.eigvals):
        option_parser.error("--axes should be between 0 and the max "
                            "number of axes available (%d), found: %d "
                            % (len(ord_res.eigvals), axes))

    if weighted:
        if sort_by == 'SampleID':
            option_parser.error("The weighting_vector can't be the SampleID, "
                                "please specify a numeric column in the "
                                "--sort_by option.")
        elif not sort_category:
            option_parser.error("To weight the output, provide a metadata "
                                "category using the --sort_by option")

    if algorithm == 'wdiff':
        if not window_size:
            option_parser.error("You should provide --window_size when using "
                                "the wdiff algorithm")
        if window_size < 1:
            option_parser.error("--window_size should be a positive integer, "
                                "%d found" % window_size)

    res = run_trajectory_analysis(ord_res, metamap, categories,
                                  sort_category, algorithm, axes,
                                  weighted, window_size)

    if not exists(output_dir):
        mkdir(output_dir)

    trajectories_fp = join(output_dir, 'trajectories.txt')
    trajectories_raw_fp = join(output_dir, 'trajectories_raw_values.txt')

    with open(trajectories_fp, 'w') as out_f:
        with open(trajectories_raw_fp, 'w') as raw_f:
            res.to_files(out_f, raw_f)
