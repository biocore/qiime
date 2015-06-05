#!/usr/bin/env python
# File created on 02 Nov 2012
from __future__ import division

__author__ = "Yoshiki Vazquez-Baeza"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Yoshiki Vazquez-Baeza", "Antonio Gonzalez-Pena"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Yoshiki Vazquez-Baeza"
__email__ = "yoshiki89@gmail.com"

from os.path import basename, splitext
from qiime.format import format_mapping_file
from qiime.parse import parse_matrix, parse_mapping_file
from qiime.util import parse_command_line_parameters, make_option
from qiime.add_alpha_to_mapping_file import (mean_alpha,
                                             add_alpha_diversity_values_to_mapping_file)

script_info = {}
script_info['brief_description'] = "Add alpha diversity data to a metadata " +\
    "mapping file"
script_info['script_description'] = "Add alpha diversity data to a mapping " +\
    "file for use with other QIIME scripts, i. e. make_emperor. The " +\
    "resulting mapping file will contain three new columns per metric in the " +\
    "alpha diversity data; the first column being the raw value, the second " +\
    "being a normalized raw value and the third one a label classifying " +\
    "the bin where this value fits based on the normalized value."
script_info['script_usage'] = [("Adding alpha diversity data:", "Add the alpha"
                                " diversity values to a mapping file and classify the normalized values "
                                "into 4 bins, where the limits will be  0 < x <= 0.25 for the first bin "
                                "0.25 < x <= 0.5 for the second bin, 0.5 < x <= 0.75 for the third bin "
                                "and 0.75 < x <= 1 for the fourth bin.", "%prog -i adiv_pd.txt -m mapping"
                                ".txt -b 4 -o alpha_mapping.txt")]
script_info['script_usage'].append(("Adding alpha diversity data with the "
                                    "quantile method:", "Add the alpha diversity values to a mapping file and"
                                    " classify the normalized values using the quartiles of the distribution "
                                    "of this values.", "%prog -i adiv_pd.txt -m mapping.txt -b 4 -o alpha_"
                                    "mapping_quantile.txt --binning_method=quantile"))
script_info['script_usage'].append(("Adding collated alpha diversity data",
                                    "Add the mean of the alpha diversity values at a specified rarefaction"
                                    " depth, this case is for use with the output of collate_alpha.py. It is "
                                    "recommended that the filenames are the name of the metric used in each "
                                    "file.", "%prog -i 'shannon.txt,chao1.txt' -m mapping.txt -b 4 -o collated_"
                                    "alpha_mapping.txt --depth=49 --collated_input"))
script_info['output_description'] = "The result of running this script is a " +\
    "metadata mapping file that will include 3 new columns per alpha " +\
    "diversity metric included in the alpha diversity file. For example, with" +\
    " an alpha diversity file with only PD_whole_tree, the new columns will " +\
    "PD_whole_tree_alpha, PD_whole_tree_normalized and PD_whole_tree_bin."
script_info['required_options'] = [
    make_option('-i', '--alpha_fps', type="existing_filepaths",
                help='alpha diversity data with one or multiple metrics i. e. the output'
                ' of alpha_diversity.py. This can also be a comma-separated list of colla'
                'ted alpha diversity file paths i. e. the output of collate_alpha.py, when'
                ' using collated alpha diversity data the --depth option is required'),
    make_option('-m', '--mapping_fp', type="existing_filepath",
                help='mapping file to modify by adding the alpha diversity data'),
]
script_info['optional_options'] = [
    make_option('-o', '--output_mapping_fp', type="new_filepath",
                help='filepath for the modified mapping file [default: %default]',
                default='mapping_file_with_alpha.txt'),
    make_option('-b', '--number_of_bins', type="int",
                help='number of bins [default: %default].', default=4),
    make_option('-x', '--missing_value_name', type="string",
                help='bin prefix name for the sample identifiers that exist in the '
                'mapping file (mapping_fp) but not in the alpha diversity file '
                '(alpha_fp) [default: %default].', default='N/A'),
    make_option(
        '--binning_method', type='choice', choices=['equal', 'quantile'],
        default='equal',
        help='Select the method name to create the bins, the options are '
        '\'equal\' and \'quantile\'. Both methods work over the normalized alpha '
        'diversity values. On the one hand \'equal\' will assign the bins on '
                'equally spaced limits, depending on the value of --number_of_bins i. e. '
                'if you select 4 the limits will be [0.25, 0.50, 0.75]. On the other hand '
                '\'quantile\' will select the limits based on the --number_of_bins i. e. '
                'the limits will be the quartiles if 4 is selected [default: %default].'),
    make_option('--depth', type='int', default=None, help='Select the rarefaction '
                'depth to use when the alpha_fps refers to collated alpha diversity file(s)'
                ' i. e. the output of collate_alpha.py. All the iterations contained at '
                'this depth will be averaged to form a single mean value [default: highest '
                'depth available].'),
    make_option('--collated_input', action='store_true', help='Use to specify that '
                'the -i option is composed of collated alpha diversity data.',
                default=False)
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    alpha_fps = opts.alpha_fps
    mapping_fp = opts.mapping_fp
    output_mapping_fp = opts.output_mapping_fp
    binning_method = opts.binning_method
    missing_value_name = opts.missing_value_name
    depth = opts.depth
    number_of_bins = opts.number_of_bins
    collated_input = opts.collated_input

    # if using collated data, make sure they specify a depth
    if collated_input:
        alpha_dict = {}

        # build up a dictionary with the filenames as keys and lines as values
        for single_alpha_fp in alpha_fps:
            alpha_dict[splitext(basename(single_alpha_fp))[0]] = open(
                single_alpha_fp, 'U').readlines()

        # format the collated data
        try:
            metrics, alpha_sample_ids, alpha_data = mean_alpha(alpha_dict,
                                                               depth)
        except ValueError as e:  # see mean_alpha for the possible exceptions
            option_parser.error(e.message)

    # when not using collated data, the user can only specify one input file
    else:
        if len(alpha_fps) > 1:
            option_parser.error('A comma-separated list of files should only be'
                                ' passed with the --alpha_fps option when using collated alpha '
                                'diversity data and also selecting a rarefaction depth with the'
                                ' --depth option.')
        else:
            metrics, alpha_sample_ids, alpha_data = parse_matrix(open(
                alpha_fps[0], 'U'))

    # parse the data from the files
    mapping_file_data, mapping_file_headers, comments = parse_mapping_file(
        open(mapping_fp, 'U'))

    # add the alpha diversity data to the mapping file
    out_mapping_file_data, out_mapping_file_headers = \
        add_alpha_diversity_values_to_mapping_file(metrics, alpha_sample_ids,
                                                   alpha_data, mapping_file_headers, mapping_file_data, number_of_bins,
                                                   binning_method, missing_value_name)

    # format the new data and write it down
    lines = format_mapping_file(
        out_mapping_file_headers,
        out_mapping_file_data)
    fd_out = open(output_mapping_fp, 'w')
    fd_out.writelines(lines)
    fd_out.close()

if __name__ == "__main__":
    main()
