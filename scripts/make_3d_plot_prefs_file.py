#!/usr/bin/env python
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso", "Jesse Stombaugh","Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Greg Caporaso"
__email__ = "gregcaporaso@gmail.com"
__status__ = "Pre-release"

from optparse import make_option
from qiime.util import parse_command_line_parameters
from qiime.format import build_prefs_string

script_description = \
"""
This is a quick-and-dirty script to write prefs files to be passed via -p to 
make_3d_plots.py. The prefs file allow for gradient coloring of continuous 
values in the 3D plots. Currently there is only one color gradient: red to
blue, because, as mentioned, this is a quick-and-dirty script. If we decide to
stick with the pref file method for defining color gradients, we'll update
this script at that time.
"""
script_usage = \
"""
To make a prefs file to be used by make_3d_plots.py The -b value passed in is
the same as that passed in via -b to make_3D_plots.py: the command delimited
list of fields that data should be included for.  For example the -b string
could be "#SampleID,Individual" and output to the file "prefs_out.txt" 
using the -p parameter.

$ python qiime_dir/scripts/make_3d_plot_prefs_file.py -b "#SampleID,Individual" -p prefs_out.txt
"""

required_options = [\
    make_option('-b','--color_by',action='store',\
          type='string',dest='color_by',help='mapping fields to color by '+\
          '[default: %default]'),\
    make_option('-p','--output_prefs_fp',action='store',\
          type='string',dest='output_prefs_fp',\
          help='path to store output file [default: %default]'),\
    ]

optional_options = []

if __name__ == "__main__":
    option_parser, opts, args = parse_command_line_parameters(
        script_description=script_description,
        script_usage=script_usage,
        version=__version__,
        required_options=required_options,
        optional_options=optional_options)
    
    out = build_prefs_string(opts.color_by)
    f = open(opts.output_prefs_fp,'w')
    f.write(out)
    f.close()
