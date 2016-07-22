#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Justin Kuczynski"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Justin Kuczynski"
__email__ = "justinak@gmail.com"


from qiime.util import parse_command_line_parameters
from qiime.util import make_option
from qiime.hierarchical_cluster import single_file_upgma, multiple_file_upgma
import os

script_info = {}
script_info['brief_description'] = """Build a UPGMA tree comparing samples"""
script_info[
    'script_description'] = """In addition to using PCoA, it can be useful to cluster samples using UPGMA (Unweighted Pair Group Method with Arithmetic mean, also known as average linkage). As with PCoA, the input to this step is a distance matrix (i.e. resulting file from beta_diversity.py)."""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""UPGMA Cluster (Single File):""",
     """To perform UPGMA clustering on a single distance matrix (e.g.: beta_div.txt, a result file from beta_diversity.py) use the following idiom:""",
     """%prog -i $PWD/beta_div.txt -o $PWD/beta_div_cluster.tre"""))
script_info['script_usage'].append(
    ("""UPGMA Cluster (Multiple Files):""",
     """The script also functions in batch mode if a folder is supplied as input. This script operates on every file in the input directory and creates a corresponding upgma tree file in the output directory, e.g.:""",
     """%prog -i $PWD/beta_div_folder -o $PWD/beta_div_folder_results/"""))
script_info[
    'output_description'] = """The output is a newick formatted tree compatible with most standard tree viewing programs. Batch processing is also available, allowing the analysis of an entire directory of distance matrices."""
script_info['required_options'] = [
    make_option('-i', '--input_path',
                type='existing_path', help='input path.  directory for batch processing, ' +
                'filename for single file operation'),


    make_option('-o', '--output_path',
                type='new_path', help='output path. directory for batch processing, ' +
                'filename for single file operation'),
]
script_info['optional_options'] = []
script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    if os.path.isdir(opts.input_path):
        multiple_file_upgma(opts.input_path, opts.output_path)
    elif os.path.isfile(opts.input_path):
        single_file_upgma(opts.input_path, opts.output_path)
    else:
        print("io error, check input file path")
        exit(1)

if __name__ == "__main__":
    main()
