#!/usr/bin/env python
# File created on 19 Jan 2011
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"


from qiime.util import make_option

import os

import warnings
warnings.filterwarnings('ignore', 'Not using MPI as mpi4py not found')

from qiime.parse import parse_distmat_to_dict, parse_mapping_file,\
    mapping_file_to_dict
from qiime.util import parse_command_line_parameters
import numpy
import matplotlib.pyplot as plt

from qiime.categorized_dist_scatterplot import get_avg_dists, get_sam_ids

script_info = {}
script_info[
    'brief_description'] = "Create a categorized distance scatterplot representing average distances between samples, broken down by categories"
script_info[
    'script_description'] = "Create a figure representing average distances between samples, broken down by categories. I call it a 'categorized distance scatterplot'. See script usage for more details. The mapping file specifies the relevant data - if you have e.g. 'N/A' values or samples you don't want included, first use filter_samples_from_otu_table.py to remove unwanted samples from the mapping file, and thus the analysis. Note that the resulting plot will include only samples in both the mapping file AND the distance matrix."
script_info['script_usage'] = [(
    "Canonical Example:",
    "Split samples by country. Within each country compare each child to all adults. Plot the average distance from that child to all adults, vs. the age of that child",
    "python categorized_dist_scatterplot.py -m map.txt -d unifrac_distance.txt -c Country -p AgeCategory:Child -s AgeCategory:Adult -a AgeYears -o fig1.png"),
    ("Example 2:",
     "Same as above, but compares Child with all other categories (e.g.: NA, Infant, etc.)",
     "python categorized_dist_scatterplot.py -m map.txt -d unifrac_distance.txt -c Country -p AgeCategory:Child -a AgeYears -o fig1.svg")]
script_info[
    'output_description'] = "a figure and the text dat for that figure "
script_info['required_options'] = [
    make_option('-m', '--map', type='existing_filepath',
                help='mapping file'),
    make_option('-d', '--distance_matrix', type='existing_filepath',
                help='distance matrix'),
    make_option('-p', '--primary_state', type='string',
                help="Samples matching this state will be plotted. E.g.: AgeCategory:Child . See qiime's filter_samples_from_otu_table.py for more syntax options"),
    make_option('-a', '--axis_category', type='string',
                help='this will form the horizontal axis of the figure, e.g.: AgeYears . Must be numbers'),
    make_option('-o', '--output_path', type='new_dirpath',
                help='output figure, filename extention determines format. E.g.: "fig1.png" or similar. A "fig1.txt" or similar will also be created with the data underlying the figure'),
]
script_info['optional_options'] = [
    make_option('-c', '--colorby', type='string',
                help='samples will first be separated by this column of the mapping file. They will be colored by this column of the mapping file, and all comparisons will be done only among samples with the same value in this column. e.g.: Country. You may omit -c, and the samples will not be separated'),
    make_option('-s', '--secondary_state', type='string',
                help='all samples matching the primary state will be compared to samples matcthing this secondary state. E.g.: AgeCategory:Adult'),
]
script_info['version'] = __version__


def main():
    option_parser, opts, args =\
        parse_command_line_parameters(**script_info)
    map_data, map_header, map_comments = parse_mapping_file(
        open(opts.map, 'U'))
    map_dict = mapping_file_to_dict(map_data, map_header)

    distdict = parse_distmat_to_dict(open(opts.distance_matrix, 'U'))

    if opts.colorby is None:
        colorby_cats = [None]
    else:
        colorby_idx = map_header.index(opts.colorby)
        colorby_cats = list(set([map_data[i][colorby_idx] for
                                 i in range(len(map_data))]))
    textfilename = os.path.splitext(opts.output_path)[0] + '.txt'
    text_fh = open(textfilename, 'w')
    text_fh.write(opts.axis_category + '\tdistance\tSampleID' + '\n')
    colorby_cats.sort()
    plt.figure()
    for cat_num, cat in enumerate(colorby_cats):
        # collect the primary and secondary samples within this category
        state1_samids, state2_samids = get_sam_ids(map_data, map_header,
                                                   opts.colorby, cat, opts.primary_state, opts.secondary_state)
        state1_samids =\
            list(set(state1_samids).intersection(set(distdict.keys())))
        state2_samids =\
            list(set(state2_samids).intersection(set(distdict.keys())))
        if state1_samids == [] or state2_samids == [] or \
                (len(state1_samids) == 1 and state1_samids == state2_samids):
            raise RuntimeError("one category of samples didn't have any valid" +
                               " distances. try eliminating samples from -p or -s, or changing" +
                               " your mapping file with filter_samples_from_otu_table.py")
        # go through dmtx
        state1_avg_dists = get_avg_dists(
            state1_samids,
            state2_samids,
            distdict)

        # plot
        xvals = [float(map_dict[sam][opts.axis_category]) for
                 sam in state1_samids]
        try:
            color = plt.cm.jet(cat_num / (len(colorby_cats) - 1))
        except ZeroDivisionError:  # only one cat
            color = 'b'
        plt.scatter(xvals, state1_avg_dists, edgecolors=color, alpha=.5,
                    facecolors='none')
        plt.xlabel(opts.axis_category)
        plt.ylabel('average distance')

        lines = [str(xvals[i]) + '\t' + str(state1_avg_dists[i]) +
                 '\t' + state1_samids[i] + '\n' for i in range(len(xvals))]
        text_fh.writelines(lines)

    if opts.colorby is not None:
        plt.legend(colorby_cats)
    plt.savefig(opts.output_path)

if __name__ == "__main__":
    main()
