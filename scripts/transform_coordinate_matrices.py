#!/usr/bin/env python
# File created on 21 Feb 2010
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Greg Caporaso", "Justin Kuczynski",
               "Jose Carlos Clemente Litran", "Jose Antonio Navas Molina"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"

from qiime.util import make_option
from os.path import split, splitext, exists, join
from os import makedirs
from qiime.util import parse_command_line_parameters
from qiime.parse import fields_to_dict
from qiime.format import format_p_value_for_num_iters
from qiime.transform_coordinate_matrices import (procrustes_monte_carlo,
                                                 get_procrustes_results)

script_info = {}
script_info['brief_description'] = "Transform two or more coordinate matrices"
script_info['script_description'] = ("This script transforms two or more "
                                     "coordinate matrices (e.g., the output of"
                                     " principal_coordinates.py) using "
                                     "procrustes analysis to minimize the "
                                     "distances between corresponding points. "
                                     "The first coordinate matrix provided is "
                                     "treated as the reference, and all other "
                                     "coordinate matrices are transformed to "
                                     "minimize distances to the reference "
                                     "points. Monte Carlo simulations can "
                                     "additionally be performed (-r random "
                                     "trials are run) to estimate the "
                                     "probability of seeing an M^2 value as "
                                     "extreme as the actual M^2.")
script_info['script_usage'] = [
    ("Write the transformed procrustes matrices to file",
     "",
     "%prog -i unweighted_unifrac_pc.txt,weighted_unifrac_pc.txt "
     "-o procrustes_output"),
    ("Generate transformed procrustes matrices and monte carlo p-values for "
     "two principal coordinate matrices",
     "",
     "%prog -i unweighted_unifrac_pc.txt,weighted_unifrac_pc.txt "
     "-o mc_procrustes_output_2 -r 1000"),
    ("Generate transformed procrustes matrices and monte carlo p-values for "
     "four principal coordinate matrices",
     "",
     "%prog -i unweighted_unifrac_pc.txt,weighted_unifrac_pc.txt,euclidean_pc."
     "txt,bray_curtis_pc.txt -o mc_procrustes_output_4 -r 1000"),
    ("Generate transformed procrustes matrices and monte carlo p-values for "
     "three principal coordinate matrices where the sample ids must be "
     "mapped between matrices",
     "",
     "%prog -i s1_pc.txt,s2_pc.txt,s3_pc.txt -s s1_s2_map.txt,s1_s3_map.txt "
     "-o mc_procrustes_output_3 -r 1000",
     )]
script_info['output_description'] = ("Two transformed coordinate matrices "
                                     "corresponding to the two input "
                                     "coordinate matrices, and (if -r was "
                                     "specified) a text file summarizing the "
                                     "results of the Monte Carlo simulations.")
script_info['required_options'] = [
    make_option('-i', '--input_fps', type='existing_filepaths',
                help='comma-separated list of input coordinate matrices'),
    make_option('-o', '--output_dir', type='new_dirpath',
                help='the output directory'),
]
script_info['optional_options'] = [
    make_option('-r', '--random_trials', type='int',
                help='Number of random permutations of matrix2 to perform. '
                     '[default: (no Monte Carlo analysis performed)]',
                default=None),
    make_option('-d', '--num_dimensions', type='int', default=3,
                help='Number of dimensions to include in output matrices '
                     '[default: %default]'),
    make_option('-s', '--sample_id_map_fps', type='existing_filepaths',
                help='If sample id maps are provided, there must be exactly '
                     'one fewer files here than there are coordinate matrices '
                     '(as each nth sample id map will provide the mapping from'
                     ' the first input coordinate matrix to the n+1th '
                     'coordinate matrix) [default: %default]',
                default=None),
    make_option('--store_trial_details',
                help='Store PC matrices for individual trials '
                     '[default: %default]',
                default=False, action='store_true'),
]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    input_fps = opts.input_fps
    sample_id_map_fps = opts.sample_id_map_fps
    num_dimensions = opts.num_dimensions
    max_dims_str = str(num_dimensions or 'alldim')
    output_dir = opts.output_dir
    random_trials = opts.random_trials

    if random_trials is not None and random_trials < 10:
        option_parser.error(
            'Must perform >= 10 trails for Monte Carlo analysis.')

    if sample_id_map_fps and \
       (len(sample_id_map_fps) + 1) != len(opts.input_fps):
        option_parser.error('If providing sample id maps, there must be '
                            'exactly one fewer sample id maps than input '
                            'coordinate matrices.')

    if not exists(output_dir):
        makedirs(output_dir)

    reference_input_fp = input_fps[0]
    reference_input_fp_dir, input_fn1 = split(reference_input_fp)
    reference_input_fp_basename, reference_input_fp_ext = splitext(input_fn1)
    output_summary_fp = join(output_dir, 'procrustes_results.txt')
    summary_file_lines = ['#FP1\tFP2\tNum included dimensions\tMonte Carlo '
                          'p-value\tCount better\tM^2',
                          '#Warning: p-values in this file are NOT currently '
                          'adjusted for multiple comparisons.']

    for i, query_input_fp in enumerate(input_fps[1:]):
        query_input_fp_dir, query_input_fn = split(query_input_fp)
        query_input_fp_basename, query_input_fp_ext = splitext(query_input_fn)
        output_matrix1_fp = join(output_dir, '%s_transformed_reference.txt'
                                 % reference_input_fp_basename)
        output_matrix2_fp = join(output_dir, '%s_transformed_q%d.txt'
                                 % (query_input_fp_basename, i + 1))

        if sample_id_map_fps:
            with open(sample_id_map_fps[i], "U") as f:
                sample_id_map = dict([(k, v[0]) for k, v in
                                     fields_to_dict(f).iteritems()])
        else:
            sample_id_map = None

        with open(reference_input_fp, 'U') as ref_in_f:
            with open(query_input_fp, 'U') as query_in_f:
                transf_coords1, transf_coords2, m_squared, rand_coords2 =\
                    get_procrustes_results(ref_in_f, query_in_f,
                                           sample_id_map=sample_id_map,
                                           randomize=False,
                                           max_dimensions=num_dimensions)

        transf_coords1.write(output_matrix1_fp)
        transf_coords2.write(output_matrix2_fp)

        if random_trials:
            if opts.store_trial_details:
                trial_output_dir = join(output_dir, 'trial_details_%d' % i + 2)
            else:
                trial_output_dir = None
            coords_f1 = open(reference_input_fp, 'U')
            coords_f2 = open(query_input_fp, 'U')
            actual_m_squared, trial_m_squareds, count_better, mc_p_value =\
                procrustes_monte_carlo(coords_f1,
                                       coords_f2,
                                       trials=random_trials,
                                       max_dimensions=num_dimensions,
                                       sample_id_map=sample_id_map,
                                       trial_output_dir=trial_output_dir)
            # truncate the p-value to the correct number of significant
            # digits
            mc_p_value_str = format_p_value_for_num_iters(mc_p_value,
                                                          random_trials)
            summary_file_lines.append('%s\t%s\t%s\t%s\t%d\t%1.3f' %
                                      (reference_input_fp, query_input_fp,
                                       max_dims_str, mc_p_value_str,
                                       count_better, actual_m_squared))
        else:
            summary_file_lines.append('%s\t%s\t%s\tNA\tNA\t%1.3f' %
                                      (reference_input_fp, query_input_fp,
                                       max_dims_str, m_squared))
    # Write output summary
    with open(output_summary_fp, 'w') as f:
        f.write('\n'.join(summary_file_lines))
        f.write('\n')


if __name__ == "__main__":
    main()
