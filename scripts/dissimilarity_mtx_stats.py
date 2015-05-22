#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski", "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"


from qiime.util import parse_command_line_parameters, matrix_stats
from qiime.util import make_option
import os
from qiime.parse import parse_distmat
from qiime.format import format_distance_matrix


script_info = {}
script_info[
    'brief_description'] = """Calculate mean, median and standard deviation from a set of distance matrices"""
script_info['script_description'] = """
This script reads in all (dis)similarity matrices from an input directory
(input_dir), then calculates and writes the mean, median, standdard deviation
(stdev) to an output folder.

The input_dir must contain only (dis)similarity matrices, and only those you
wish to perform statistical analyses on.
"""
script_info['script_usage'] = []
script_info['script_usage'].append(("Example",
                                    "This examples takes the \"dists/\" directory as input and returns the "
                                    "results in the \"dist_stats/\" directory.",
                                    "%prog -i dists/ -o dist_stats/"))
script_info['output_description'] = """
The outputs are in distance matrix format, where each value is the mean,
median, or stdev of that element in all the input distance matrices.
"""
script_info['required_options'] = [
    make_option('-i', '--input_dir', type='existing_dirpath',
                help='Path to input directory'),
    make_option('-o', '--output_dir', type='new_dirpath',
                help='Path to store result files')
]
script_info['optional_options'] = []
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    indir = opts.input_dir
    outdir = opts.output_dir
    if not os.path.exists(outdir):
        os.makedirs(outdir)

    # input
    file_names = os.listdir(indir)
    file_names = [fname for fname in file_names if not fname.startswith('.')]

    distmats = []
    headers_list = []
    for fname in file_names:
        f = open(os.path.join(indir, fname), 'U')
        headers, data = parse_distmat(f)
        f.close()
        distmats.append(data)
        headers_list.append(headers)

    # calcs
    headers, means, medians, stdevs = matrix_stats(headers_list, distmats)

    # output
    f = open(os.path.join(outdir, 'means.txt'), 'w')
    f.write(format_distance_matrix(headers, means))
    f.close()

    f = open(os.path.join(outdir, 'medians.txt'), 'w')
    f.write(format_distance_matrix(headers, medians))
    f.close()

    f = open(os.path.join(outdir, 'stdevs.txt'), 'w')
    f.write(format_distance_matrix(headers, stdevs))
    f.close()

if __name__ == "__main__":
    main()
