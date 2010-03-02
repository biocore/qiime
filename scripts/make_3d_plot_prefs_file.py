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

script_info={}
script_info['brief_description']="""Generate preferences file for 3D plots using Metadata Fields"""
script_info['script_description']="""This script generates a preferences (prefs) file, which can be passed via -p to make_3d_plots.py. The prefs file allows for gradient coloring of continuous values in the 3D plots. Currently there is only one color gradient: red to blue."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Example:""","""To make a prefs file to be used by make_3d_plots.py, the -b value should be passed in as the same as that passed in via -b to make_3d_plots.py. For multiple fields, the command should be a delimited list of fields. For example the -b string could be "#SampleID,Individual" and output to the file "prefs_out.txt" using the -p parameter.""","""make_3d_plot_prefs_file.py -b "#SampleID,Individual" -p prefs_out.txt"""))
script_info['output_description']="""The result of this script is a text file, containing coloring preferences to be used by make_3d_plots.py."""
script_info['optional_options']=[]


script_info['required_options']=[\
    make_option('-b','--color_by',action='store',\
          type='string',dest='color_by',help='mapping fields to color by '+\
          '[default: %default]'),\
    make_option('-p','--output_prefs_fp',action='store',\
          type='string',dest='output_prefs_fp',\
          help='path to store output file [default: %default]'),\
]

script_info['version']=__version__

if __name__ == "__main__":
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    out = build_prefs_string(opts.color_by)
    f = open(opts.output_prefs_fp,'w')
    f.write(out)
    f.close()
