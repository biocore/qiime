#!/usr/bin/env python
# File created on 3 May 2011
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"

from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from qiime.nmds import nmds, multiple_file_nmds
from qiime.parse import parse_distmat
import os

options_lookup = get_options_lookup()

script_info = {}
script_info[
    'brief_description'] = """Nonmetric Multidimensional Scaling (NMDS)"""
script_info[
    'script_description'] = """Nonmetric Multidimensional Scaling (NMDS) is commonly used to compare groups of samples based on phylogenetic or count-based distance metrics (see section on beta_diversity.py)."""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""NMDS (Single File)""",
     """For this script, the user supplies a distance matrix (i.e. resulting file from beta_diversity.py), along with the output filename (e.g. beta_div_coords.txt), as follows:""",
     """%prog -i beta_div.txt -o beta_div_coords.txt"""))
script_info['script_usage'].append(
    ("""NMDS (Dimensions)""",
     """For this script, the user supplies a distance matrix (i.e. resulting file from beta_diversity.py), the number of dimensions of NMDS space and the output filename (e.g. beta_div_coords.txt), as follows:""",
     """%prog -i beta_div.txt -d 3 -o beta_div_3_coords.txt"""))
script_info['script_usage'].append(
    ("""NMDS (Multiple Files):""",
     """The script also functions in batch mode if a folder is supplied as input (e.g. from beta_diversity.py run in batch). No other files should be present in the input folder - only the distance matrix files to be analyzed. This script operates on every distance matrix file in the input directory and creates a corresponding nmds results file in the output directory, e.g.:""",
     """%prog -i beta_div_weighted_unifrac/ -o beta_div_weighted_nmds_results/"""))
script_info[
    'output_description'] = """The resulting output file consists of the NMDS axes (columns) for each sample (rows). Pairs of NMDS axes can then be graphed to view the relationships between samples. The bottom of the output file contains the stress of the ordination."""
script_info['required_options'] = [

    make_option('-i', '--input_path', type='existing_path',
                help='path to the input distance matrix file(s) (i.e., the output from beta_diversity.py). Is a directory for batch processing and a filename for a single file operation.'),

    make_option('-o', '--output_path', type='new_path',
                help='output path. directory for batch processing, ' +
                'filename for single file operation'),
]

script_info['optional_options'] = [
    make_option('-d', '--dimensions', default=3, type='int',
                help='number of dimensions of NMDS space' +
                'default: %default'),
]
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if os.path.isdir(opts.input_path):
        multiple_file_nmds(opts.input_path, opts.output_path, opts.dimensions)
    elif os.path.isfile(opts.input_path):
        f = open(opts.input_path, 'U')
        nmds_res_string = nmds(f, opts.dimensions)
        f.close()

        f = open(opts.output_path, 'w')
        f.write(nmds_res_string)
        f.close()
    else:
        print("io error, check input file path")
        exit(1)


if __name__ == "__main__":
    main()
