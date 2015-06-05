#!/usr/bin/env python
from __future__ import division

__author__ = "Dan Knights"
__copyright__ = "Copyright 2012, The QIIME Project"
__credits__ = ["Dan Knights"]
__license__ = "GPL"
__version__ = "1.9.1-dev"
__maintainer__ = "Dan Knights"
__email__ = "danknights@gmail.com"

from qiime.util import make_option
from qiime.util import parse_command_line_parameters

from qiime.detrend import detrend_pcoa, validate_metadata

import warnings

script_info = {}
script_info['brief_description'] = """Detrend Principal Coordinates"""
script_info[
    'script_description'] = """Ordination plots (e.g. principal coordinates analysis) of samples that lay along a naturally occurring gradient (e.g. depth, time, pH) often exhibit a curved shape known as the "arch" or "horseshoe" effect. This can cause samples near the endpoints of the gradient to appear closer to one another than would be expected. This script will attempt to remove any (compounded) quadratic curvature in a set of 2D coordinates. If requested, it will also report an evaluation of the association of the transformed coordinates with a known gradient."""
script_info['script_usage'] = []
script_info['script_usage'].append(
    ("""Examples:""",
     """The simplest usage takes as input only a table of principal coordinates:""",
     """%prog -i $PWD/pcoa.txt -o detrending"""))
script_info['script_usage'].append(
    ("""""",
     """One may also include a metadata file with a known real-valued gradient as one of the columns. In this case, the output folder will include a text file providing a summary of how well the analysis fit with the hypothesis that the primary variation is due to the gradient (in this case, "DEPTH"):""",
     """%prog -i $PWD/pcoa.txt -m map.txt -c DEPTH -o detrending"""))
script_info['script_usage'].append(
    ("""""",
     """Note that if you provide a real-valued known gradient the script will prerotate the first two axes of the PCoA coords in order to achieve optimal alignment with that gradient. This can be disabled with "-r":""",
     """%prog -i $PWD/pcoa.txt -m map.txt -c DEPTH -o detrending -r"""))

script_info['required_options'] = [
    make_option('-i', '--input_fp', action='store',
                type='existing_filepath', dest='input_fp', help='Path to read ' +
                'PCoA/PCA/ordination table')
]

script_info['optional_options'] = [

    make_option('-o', '--output_dir', action='store', type='new_filepath',
                help='Path to output directory' +
                ' [default: .]'),

    make_option('-m', '--map_fp', action='store', type='existing_filepath',
                default=None, help='Path to metadata file' +
                ' [default: None]'),

    make_option('-c', '--gradient_variable', action='store',
                default=None, help='Column header for gradient variable in metadata table' +
                ' [default: None]'),

    make_option('-r', '--suppress_prerotate', action='store_true',
                default=False,
                help="Suppress pre-rotation of the coordinates for optimal detrending;" +
                " not pre-rotating assumes that the curvature is symmetrical across" +
                " the vertical axis [default: %default]"),
]
script_info['output_description'] = """Output is detrended PCoA matrices."""
script_info['version'] = __version__


def main():
    # parse command line
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    # can't pass map without column header
    if opts.map_fp is not None and opts.gradient_variable is None:
        option_parser.error(
            "Cannot pass -m/--map_fp without -c/--gradient_variable.")

    # ensure metadata is present and real-valued
    if opts.map_fp is not None:
        validate_metadata(opts.map_fp, opts.gradient_variable)

    # run detrending
    detrend_pcoa(opts.input_fp,
                 map_fp=opts.map_fp,
                 gradient_variable=opts.gradient_variable,
                 suppress_prerotate=opts.suppress_prerotate,
                 output_dir=opts.output_dir)

if __name__ == "__main__":
    main()
