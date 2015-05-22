#!/usr/bin/env python
# File created on 10 Mar 2011
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"


from qiime.util import parse_command_line_parameters
from qiime.util import make_option
from qiime.parse import parse_mapping_file_to_dict, parse_distmat
import os
import sys
import numpy

from qiime.cluster_quality import clust_qual_ratio

script_info = {}
script_info['brief_description'] = """compute the quality of a cluster"""
script_info[
    'script_description'] = """The input is a distance matrix (i.e. resulting file from beta_diversity.py)."""

script_info['script_usage'] = []

script_info['script_usage'].append(
    ("""cluster quality based on the treatment category:""",
     """to compute the quality of clusters, and print to stdout, use the following idiom:""",
     """%prog -i weighted_unifrac_otu_table.txt -m Fasting_Map.txt -c Treatment"""))

script_info['script_usage'].append(
    ("""cluster quality based on the DOB category:""",
     """to compute the quality of clusters, and print to stdout, use the following idiom:""",
     """%prog -i weighted_unifrac_otu_table.txt -m Fasting_Map.txt -c DOB"""))


script_info[
    'output_description'] = """The output is either a single number (with -s), or a more detailed output of the similarity between and within clusters."""
script_info['required_options'] = [
    make_option('-i', '--input_path', type='existing_filepath',
                help='input distance matrix file'),

    make_option('-m', '--map', type='existing_filepath',
                help='mapping file'),

    make_option('-c', '--category', type='string',
                help='column of mapping file delimiting clusters'),



]
script_info['optional_options'] = [

    make_option('-o', '--output_path', default=None, type='new_filepath',
                help='output path, prints to stdout if omitted'),
    make_option('-s', '--short', action="store_true",
                help='print only ' +
                'the ratio of mean dissimilarities between/within clusters' +
                ' instead of more detailed output'),

    make_option('--metric', default='ratio', type='choice',
                choices=['ratio'],
                help='choice of quality metric to apply. Currently only one option ' +
                'exists, the ratio of mean(distances between samples from different ' +
                'clusters) to mean(distances between samples from the same cluster) ' +
                'Default: %default'),
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if opts.output_path is not None:
        outf = open(opts.output_path, 'w')
    else:
        outf = sys.stdout

    dists = parse_distmat(open(opts.input_path, 'U'))
    map_data = parse_mapping_file_to_dict(open(opts.map, 'U'))
    diff_dists, same_dists = clust_qual_ratio(dists, map_data, opts.category)

    if opts.short:
        print >> outf, numpy.mean(diff_dists) / numpy.mean(same_dists)
    else:
        print >> outf, "dissimilarity ratio between/within (large for clustered data):"
        print >> outf, numpy.mean(diff_dists) / numpy.mean(same_dists)
        print >> outf, "dissimilarities between clusters: mean, std, num:"
        print >> outf, '\t'.join(map(str, [numpy.mean(diff_dists), numpy.std(diff_dists),
                                           len(diff_dists)]))
        print >> outf, "dissimilarities within clusters: mean, std, num:"
        print >> outf, '\t'.join(map(str, [numpy.mean(same_dists), numpy.std(same_dists),
                                           len(same_dists)]))


if __name__ == "__main__":
    main()
