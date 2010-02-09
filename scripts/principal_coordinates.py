#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Justin Kuczynski", "Rob Knight", "Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Pre-release"
 
from qiime.util import parse_command_line_parameters
from optparse import make_option
from qiime.principal_coordinates import pcoa
from qiime.parse import parse_distmat

script_description = """ Contains code for the pricipal coordinates analysis (PCoA).

This module has the responsibility of calculating the PCoA coordinates of a distance matrix."""

script_usage = """
principal_coordinates.py -i unifrac_dist_mtx.txt -o unifrac_pcoa.txt
"""

required_options = [\
 make_option('-i', '--input_path',
        help='input filepath. [REQUIRED]'),\
 make_option('-o', '--output_path',
        help='output filepath. [REQUIRED]'),\
]

optional_options = [\
]

def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)
      
    infilepath = opts.input_path
    f = open(infilepath,'U')
    pcoa_res_string = pcoa(f)
    f.close()

    f = open(opts.output_path, 'w')
    f.write(pcoa_res_string)
    f.close()


if __name__ == "__main__":
    main()