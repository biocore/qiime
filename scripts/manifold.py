#!/usr/bin/env python

__author__ = "Joshua Haas"
__copyright__ = "Copyright 2014, The QIIME Project"
__credits__ = ["Joshua Haas","Antonio Gonzalez Pena","Gregory Ditzler"]
__license__ = "GPL"
__version__ = "1.8.0-dev"
__maintainer__ = "Joshua Haas"
__email__ = "laptopdude2@gmail.com"


"""
Command line interface for the manifold script.
"""

from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from qiime.manifold import compute_manifold

import os

options_lookup = get_options_lookup()

script_info = {}
script_info['brief_description'] = """Manifold Learning Dimensionality Reduction"""
script_info[
    'script_description'] = """Manifold learning techniques are used to visualize patterns and structures inherent in high-dimensional data. By coloring the data using metadata, clustering can reveal pertinent relationships and dependencies."""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Manifold (Single File)""",
     """For this script, the user supplies a .biom file, along with the output filename (e.g. isomap.txt), as follows:""",
     """%prog -i otu_table.biom -o isomap.txt"""))
script_info['script_usage'].append(
    ("""Manifold (Multiple Files):""",
     """The script also functions in batch mode if a folder is supplied as input. No other files should be present in the input folder - only the .biom files to be analyzed. This script operates on every .biom file in the input directory and creates a corresponding results file in the output directory, e.g.:""",
     """%prog -i otu_tables/ -o manifold_results/"""))
script_info[
    'output_description'] = """The resulting output file consists of the computed axes (columns) for each sample (rows). Pairs of axes can then be graphed to view the relationships between samples. The bottom of the output file contains the eigenvalues and % variation explained for each PC."""
script_info['required_options'] = [

    make_option('-i', '--input_path', type='existing_path',
                help='path to the input .biom file(s). Is a directory for batch processing and a filename for a single file operation.'),

    make_option('-o', '--output_path', type='new_path',
                help='output path. Is a directory for batch processing and a ' +
                'filename for a single file operation'),

    make_option('-a', '--algorithm', type='string',
                help='algorithm. Which manifold technique to use. Implemented algorithms: ' +

                'isomap, ' +
                'lle, ' +
                'spectral, ' +
                'ltsa, ' +
                'mds. ' +
                
                'Note that lle allows "standard", "modified", and "hessian" via the "method" paramter.')
]

script_info['optional_options'] = [

    make_option('-p','--params', type='string',
                help='parameters for the chosen algorithm. For a full list of parameters see the Scikit Learn documentation for a full list.' +
                'This option should be in the form "param1=value1,param2=value2". Note that some options do not apply to certain algorithms.')
]

script_info['version'] = __version__


def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    #Individual parameters are separated by commas
    params_list = opts.params.split(",")
    params = dict()
    for pair in params_list:
        #The value is separated from the parameter by an equals sign
        key_value = pair.split("=")
        params[key_value[0]] = key_value[1]

    if os.path.isdir(opts.input_path):
        multiple_file_manifold(opts.input_path, opts.output_path, opts.algorithm, params)
    elif os.path.isfile(opts.input_path):
        manifold_res_string = compute_manifold(opts.input_path, opts.algorithm, params)

        f = open(opts.output_path, 'w')
        f.write(manifold_res_string)
        f.close()
    else:
        print("io error, check input file path")
        exit(1)


if __name__ == "__main__":
    main()
