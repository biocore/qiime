#!/usr/bin/env python
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Jai Ram Rideout", "Michael Dwan", "Logan Knecht",
               "Damien Coy", "Levi McCracken", "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Jai Ram Rideout"
__email__ = "jai.rideout@gmail.com"

from os import path
from skbio.util import create_dir
from qiime.parse import fields_to_dict, parse_distmat
from qiime.util import (parse_command_line_parameters,
                        get_options_lookup,
                        make_option)
from qiime.compare_distance_matrices import (run_mantel_correlogram,
                                             run_mantel_test)

options_lookup = get_options_lookup()

script_info = {}
script_info[
    'brief_description'] = """Computes Mantel correlation tests between sets of distance matrices"""
script_info['script_description'] = """
This script compares two or more distance/dissimilarity matrices for \
correlation by providing the Mantel, partial Mantel, and Mantel correlogram \
matrix correlation tests.

The Mantel test will test the correlation between two matrices. The data \
often represents the "distance" between objects or samples.

The partial Mantel test is a first-order correlation analysis that utilizes \
three distance (dissimilarity) matrices. This test builds on the traditional \
Mantel test which is a procedure that tests the hypothesis that distances \
between the objects within a given matrix are linearly independent of the \
distances withing those same objects in a separate matrix. It builds on the \
traditional Mantel test by adding a third "control" matrix.

Mantel correlogram produces a plot of distance classes versus Mantel \
statistics. Briefly, an ecological distance matrix (e.g. UniFrac distance \
matrix) and a second distance matrix (e.g. spatial distances, pH distances, \
etc.) are provided. The second distance matrix has its distances split into a \
number of distance classes (the number of classes is determined by Sturge's \
rule). A Mantel test is run over these distance classes versus the ecological \
distance matrix. The Mantel statistics obtained from each of these tests are \
then plotted in a correlogram. A filled-in point on the plot indicates that \
the Mantel statistic was statistically significant (you may provide what \
alpha to use).

For more information and examples pertaining to this script, please refer to \
the accompanying tutorial, which can be found at \
http://qiime.org/tutorials/distance_matrix_comparison.html.
"""
script_info['script_usage'] = []
script_info['script_usage'].append(("Partial Mantel",
                                    "Performs a partial Mantel test on two distance matrices, "
                                    "using a third matrix as a control. Runs 99 permutations to calculate the "
                                    "p-value.",
                                    "%prog --method partial_mantel -i "
                                    "weighted_unifrac_dm.txt,unweighted_unifrac_dm.txt -c PH_dm.txt "
                                    "-o partial_mantel_out -n 99"))
script_info['script_usage'].append(("Mantel",
                                    "Performs a Mantel test on all pairs of four distance matrices, "
                                    "including 999 permutations for each test.",
                                    "%prog --method mantel "
                                    "-i weighted_unifrac_dm.txt,unweighted_unifrac_dm.txt,"
                                    "weighted_unifrac_even100_dm.txt,unweighted_unifrac_even100_dm.txt "
                                    "-o mantel_out -n 999"))
script_info['script_usage'].append(("Mantel Correlogram",
                                    "This example computes a Mantel correlogram on two distance matrices "
                                    "using 999 permutations in each Mantel test. Output is written to the "
                                    "mantel_correlogram_out directory.",
                                    "%prog --method mantel_corr -i unweighted_unifrac_dm.txt,PH_dm.txt -o "
                                    "mantel_correlogram_out -n 999"))

script_info['output_description'] = """
Mantel: One file is created containing the Mantel 'r' statistic and p-value.

Partial Mantel: One file is created in the output directory, which contains \
the partial Mantel statistic and p-value.

Mantel Correlogram: Two files are created in the output directory: a text \
file containing information about the distance classes, their associated \
Mantel statistics and p-values, etc. and an image of the correlogram plot.
"""

script_info['required_options'] = [
    # All methods use these
    make_option('--method',
                help='matrix correlation method to use. Valid options: '
                '[mantel, partial_mantel, mantel_corr]',
                type='choice',
                choices=['mantel', 'partial_mantel', 'mantel_corr']),
    make_option('-i', '--input_dms', type='existing_filepaths',
                help='the input distance matrices, comma-separated. WARNING: Only '
                'symmetric, hollow distance matrices may be used as input. Asymmetric '
                'distance matrices, such as those obtained by the UniFrac Gain metric '
                '(i.e. beta_diversity.py -m unifrac_g), should not be used as input'),
    options_lookup['output_dir']
]
script_info['optional_options'] = [
    # All methods use these
    make_option('-n', '--num_permutations',
                help='the number of permutations to perform when calculating the '
                'p-value [default: %default]', default=100, type='int'),
    make_option('-s', '--sample_id_map_fp', type='existing_filepath',
                help='Map of original sample ids to new sample ids [default: '
                '%default]', default=None),
    # Standard Mantel specific, i.e., method == mantel
    make_option('-t', '--tail_type',
                help='the type of tail test to perform when calculating the p-value. '
                'Valid options: [two-sided, less, greater]. "two-sided" is a two-tailed '
                'test, while "less" tests for r statistics less than the observed r '
                'statistic, and "greater" tests for r statistics greater than the '
                'observed r statistic. Only applies when method is mantel [default: '
                '%default]', default='two-sided', type='choice',
                choices=['two-sided', 'greater', 'less']),
    # Mantel Correlogram specific, i.e., method == mantel_corr
    make_option('-a', '--alpha',
                help='the value of alpha to use when denoting significance in the '
                'correlogram plot. Only applies when method is mantel_corr',
                default=0.05, type='float'),
    make_option('-g', '--image_type',
                help='the type of image to produce. Valid options: [png, svg, pdf]. '
                'Only applies when method is mantel_corr [default: %default]',
                default='pdf', type='choice', choices=['pdf', 'png', 'svg']),
    make_option('--variable_size_distance_classes', action='store_true',
                help='if this option is supplied, each distance class will have an '
                'equal number of distances (i.e. pairwise comparisons), which may '
                'result in variable sizes of distance classes (i.e. each distance '
                'class may span a different range of distances). If this option is '
                'not supplied, each distance class will have the same width, but may '
                'contain varying numbers of pairwise distances in each class. This '
                'option can help maintain statistical power if there are large '
                'differences in the number of distances in each class. See '
                'Darcy et al. 2011 (PLoS ONE) for an example of this type of '
                'correlogram. Only applies when method is mantel_corr '
                '[default: %default]', default=False),
    # Partial Mantel specific, i.e., method == partial_mantel
    make_option('-c', '--control_dm',
                help='the control matrix. Only applies (and is *required*) when '
                'method is partial_mantel. [default: %default]', default=None,
                type='existing_filepath')
]
script_info['version'] = __version__

comment_mantel_pmantel = """\
# Number of entries refers to the number of rows (or cols) retained in each
# distance matrix after filtering the distance matrices to include only those
# samples that were in both distance matrices. p-value contains the correct
# number of significant digits.
"""

comment_corr = comment_mantel_pmantel[:-1] + \
    """
# Distance classes with values of None were in the second half of the distance
# classes and not all samples could be included in the distance class, so
# calculations were not performed.
"""


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if opts.num_permutations < 1:
        option_parser.error(
            "--num_permutations must be greater than or equal to 1.")

    # Create the output dir if it doesn't already exist.
    try:
        if not path.exists(opts.output_dir):
            create_dir(opts.output_dir)
    except:
        option_parser.error("Could not create or access output directory "
                            "specified with the -o option.")
    sample_id_map = None
    if opts.sample_id_map_fp:
        sample_id_map = dict([(k, v[0])
                              for k, v in fields_to_dict(open(opts.sample_id_map_fp, "U")).items()])
    input_dm_fps = opts.input_dms
    distmats = [parse_distmat(open(dm_fp, 'U')) for dm_fp in input_dm_fps]

    if opts.method == 'mantel':
        output_f = open(path.join(opts.output_dir, 'mantel_results.txt'), 'w')
        output_f.write(run_mantel_test('mantel', input_dm_fps, distmats,
                       opts.num_permutations, opts.tail_type,
                       comment_mantel_pmantel, sample_id_map=sample_id_map))
    elif opts.method == 'partial_mantel':
        output_f = open(path.join(opts.output_dir,
                        'partial_mantel_results.txt'), 'w')
        output_f.write(run_mantel_test('partial_mantel', input_dm_fps,
                       distmats, opts.num_permutations, opts.tail_type,
                       comment_mantel_pmantel, control_dm_fp=opts.control_dm,
                       control_dm=parse_distmat(open(opts.control_dm, 'U')),
                       sample_id_map=sample_id_map))
    elif opts.method == 'mantel_corr':
        output_f = open(path.join(opts.output_dir,
                        'mantel_correlogram_results.txt'), 'w')
        result_str, correlogram_fps, correlograms = run_mantel_correlogram(
            input_dm_fps, distmats, opts.num_permutations, comment_corr,
            opts.alpha, sample_id_map=sample_id_map,
            variable_size_distance_classes=opts.variable_size_distance_classes)

        output_f.write(result_str)
        for corr_fp, corr in zip(correlogram_fps, correlograms):
            corr.savefig(path.join(opts.output_dir, corr_fp + opts.image_type),
                         format=opts.image_type)
    output_f.close()

if __name__ == "__main__":
    main()
