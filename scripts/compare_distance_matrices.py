#!/usr/bin/env python
from __future__ import division

__author__ = "Michael Dwan"
__copyright__ = "Copyright 2012, The QIIME project"
__credits__ = ["Michael Dwan", "Jai Ram Rideout", "Logan Knecht",
               "Greg Caporaso"]
__license__ = "GPL"
__version__ = "1.4.0-dev"
__maintainer__ = "Michael Dwan"
__email__ = "mdwan.tgen@gmail.com"
__status__ = "Development"

from os import path

from cogent.util.misc import create_dir

from qiime.format import format_p_value_for_num_iters
from qiime.parse import fields_to_dict, parse_distmat
from qiime.util import (parse_command_line_parameters,
                        get_options_lookup,
                        make_compatible_distance_matrices,
                        make_option)
from qiime.stats import Mantel, MantelCorrelogram, PartialMantel
from qiime.util import DistanceMatrix

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = """
Computes Mantel correlation tests between sets of distance matrices
"""
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
"""
script_info['script_usage'] = []
script_info['script_usage'].append(("Partial Mantel",
"Performs a partial Mantel test on two distance matrices, "
"using a third matrix as a control. Runs 100 permutations to calculate the "
"p-value.",
"%prog --method partial_mantel -i "
"weighted_unifrac_dm.txt,unweighted_unifrac_dm.txt -c PH_dm.txt "
"-o mantel_out -n 100"))
script_info['script_usage'].append(("Mantel",
"Performs a Mantel test on all pairs of four distance matrices, "
"including 1000 permutations for each test.",
"%prog --method mantel "
"-i weighted_unifrac_dm.txt,unweighted_unifrac_dm.txt,"
"weighted_unifrac_even100_dm.txt,unweighted_unifrac_even100_dm.txt "
"-o mantel_out -n 1000"))
script_info['script_usage'].append(("Mantel Correlogram",
"This example computes a Mantel correlogram on two distance matrices "
"using 999 permutations in each Mantel test. Output is written to the "
"mantel_output directory.",
"%prog --method mantel_corr -i unweighted_unifrac_dm.txt,PH_dm.txt -o "
"mantel_output -n 999"))

script_info['output_description']= """
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
     make_option('-i','--input_dms',
        help='the input distance matrices, comma-separated'),
     options_lookup['output_dir']
]
script_info['optional_options'] = [
    # All methods use these
    make_option('-n','--num_permutations',
        help='the number of permutations to perform when calculating the '
        'p-value [default: %default]', default=100, type='int'),
     make_option('-s','--sample_id_map_fp',
        help='Map of original sample ids to new sample ids [default: '
        '%default]', default=None),
    # Standard Mantel specific, i.e., method == mantel
    make_option('-t','--tail_type',
        help='the type of tail test to perform when calculating the p-value. '
        'Valid options: [two sided, less, greater] Two-sided is a two-tailed '
        'test, while less tests for r statistics less than the observed r '
        'statistic, and greater tests for r statistics greater than the '
        'observed r statistic. Only applies when method is mantel [default: '
        '%default]', default='two sided'),
    # Mantel Correlogram specific, i.e., method == mantel_corr
    make_option('-a', '--alpha',
        help='the value of alpha to use when denoting significance in the '
        'correlogram plot. Only applies when method is mantel_corr',
        default=0.05, type='float'),
    make_option('-g', '--image_type',
        help='the type of image to produce. Valid options: [png, svg, pdf]. '
        'Only applies when method is mantel_corr [default: %default]',
        default='pdf', type="choice", choices=['pdf', 'png', 'svg']),
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

    # Create the output dir if it doesn't already exist.
    try:
        if not path.exists(opts.output_dir):
            create_dir(opts.output_dir)
    except:
        option_parser.error("Could not create or access output directory "
                            "specified with the -o option.")

    sample_id_map_fp = opts.sample_id_map_fp
    if sample_id_map_fp:
        sample_id_map = dict([(k, v[0]) \
        for k,v in fields_to_dict(open(sample_id_map_fp, "U")).items()])
    else:
        sample_id_map = None

    input_dm_fps = opts.input_dms.split(',')
    num_perms = opts.num_permutations
    method = opts.method

    if method == 'mantel':
        results_fp = 'mantel_results.txt'
        comment = comment_mantel_pmantel
        header = 'DM1\tDM2\tNumber of entries\tMantel r statistic\t' + \
                 'p-value\tNumber of permutations\tTail type\n'
    elif method == 'partial_mantel':
        results_fp = 'partial_mantel_results.txt'
        comment = comment_mantel_pmantel
        header = 'DM1\tDM2\tCDM\tNumber of entries\tMantel r statistic\t' + \
                 'p-value\tNumber of permutations\n'
    elif method == 'mantel_corr':
        results_fp = 'mantel_correlogram_results.txt'
        comment = comment_corr
        header = 'DM1\tDM2\tNumber of entries\tNumber of permutations\t' + \
                 'Class index\tNumber of distances\tMantel r statistic\t' + \
                 'p-value\tp-value (Bonferroni corrected)\n'

    # Write the comment string and header.
    output_f = open(path.join(opts.output_dir, results_fp), 'w')
    output_f.write(comment + header)

    # Loop over all pairs of dms.
    for i, fp1 in enumerate(input_dm_fps):
        for fp2 in input_dm_fps[i+1:]:
            # Open the current pair of distance matrices and make them
            # compatible by only keeping samples that match between them,
            # and ordering them by the same sample IDs.
            (dm1_labels, dm1_data), (dm2_labels, dm2_data) = \
                make_compatible_distance_matrices(
                    parse_distmat(open(fp1, 'U')),
                    parse_distmat(open(fp2, 'U')),
                    lookup=sample_id_map)
            if method == 'partial_mantel':
                if not opts.control_dm:
                    option_parser.error("You must provide a control matrix "
                            "when running the partial Mantel test.")

                # We need to intersect three sets (three matrices).
                (dm1_labels, dm1_data), (cdm_labels, cdm_data) = \
                    make_compatible_distance_matrices((dm1_labels, dm1_data),
                    parse_distmat(open(opts.control_dm, 'U')),
                        lookup=sample_id_map)
                (dm1_labels, dm1_data), (dm2_labels, dm2_data) = \
                    make_compatible_distance_matrices((dm1_labels, dm1_data),
                        (dm2_labels, dm2_data), lookup=sample_id_map)
            if len(dm1_labels) < 3:
                output_f.write('%s\t%s\t%d\tToo few samples\n' % (fp1,
                               fp2, len(dm1_labels)))
                continue
            # Create DistanceMatrix instances from our raw distance matrix
            # variables.
            dm1 = DistanceMatrix(dm1_data, dm1_labels, dm1_labels)
            dm2 = DistanceMatrix(dm2_data, dm2_labels, dm2_labels)

            # Create an instance of our correlation test and run it with the
            # specified number of permutations.
            if method == 'mantel':
                results = Mantel(dm1, dm2, opts.tail_type)(num_perms)
                p_str = format_p_value_for_num_iters(results['p_value'],
                                                     num_perms)
                output_f.write("%s\t%s\t%d\t%.5f\t%s\t%d\t%s\n" %
                        (fp1, fp2, len(dm1_labels), results['r_value'], p_str,
                         num_perms, opts.tail_type))
            elif method == 'partial_mantel':
                results = PartialMantel(dm1, dm2, DistanceMatrix(cdm_data,
                                        cdm_labels, cdm_labels))(num_perms)
                p_str = format_p_value_for_num_iters(results['mantel_p'],
                                                     num_perms)
                output_f.write("%s\t%s\t%s\t%d\t%.5f\t%s\t%d\n" %
                        (fp1, fp2, opts.control_dm, len(dm1_labels),
                         results['mantel_r'], p_str, num_perms))
            elif method == 'mantel_corr':
                results = MantelCorrelogram(dm1, dm2,
                                            alpha=opts.alpha)(num_perms)
                # Write the correlogram plot to a file.
                dm1_name = path.basename(fp1)
                dm2_name = path.basename(fp2)
                fig_file_name = '_'.join((dm1_name, 'AND', dm2_name,
                        'mantel_correlogram')) + '.' + opts.image_type
                results['correlogram_plot'].savefig(path.join(opts.output_dir,
                        fig_file_name), format=opts.image_type)

                # Iterate over the results and write them to the text file.
                first_time = True
                for class_idx, num_dist, r, p, p_corr in \
                        zip(results['class_index'], results['num_dist'],
                            results['mantel_r'], results['mantel_p'],
                            results['mantel_p_corr']):
                    p_str = None
                    if p is not None:
                        p_str = format_p_value_for_num_iters(p, num_perms)
                    p_corr_str = None
                    if p_corr is not None:
                        p_corr_str = format_p_value_for_num_iters(p_corr,
                                                                  num_perms)
                    if first_time:
                        output_f.write('%s\t%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\n'
                                % (fp1, fp2, len(dm1_labels), num_perms,
                                   class_idx, num_dist, r, p_str, p_corr_str))
                        first_time = False
                    else:
                        output_f.write('\t\t\t\t%s\t%d\t%s\t%s\t%s\n'
                                % (class_idx, num_dist, r, p_str, p_corr_str))
    output_f.close()


if __name__ == "__main__":
    main()
