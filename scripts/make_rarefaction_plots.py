#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Meg Pirrung"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Meg Pirrung"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Meg Pirrung"
__email__ = "meg.pirrung@colorado.edu"
__status__ = "Pre-release"
 
from optparse import make_option
from qiime.util import parse_command_line_parameters
from qiime.pycogent_backports.misc import get_random_directory_name
import sys
from sys import argv, exit, exc_info
from random import choice, randrange
from time import strftime
from qiime import parse, util
from qiime.make_rarefaction_plots import make_plots, make_output_files, \
parse_rarefaction_data
import os.path
from os.path import exists, splitext, split
import shutil


#make_rarefaction_plots.py
script_info={}
script_info['brief_description']="""Generate Rarefaction Plots"""
script_info['script_description']="""Once the batch alpha diversity files have been collated, you may want to compare the diversity using plots. Using the results from make_rarefaction_averages.py, you can plot the samples and or by category in the mapping file using this script.

This script creates an html file of rarefaction plots based on the supplied rarefaction files in the folder given (-i) from make_rarefaction_averages.py. The user may also supply optional arguments like an image type (-i), and a resolution (-d)."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Default Example:""","""For generated rarefaction plots using the default parameters, including the mapping file and one rarefaction file, you can use the following command:""","""make_rarefaction_plots.py -r chao1/"""))
script_info['script_usage'].append(("""Specify Image Type and Resolution:""","""Optionally, you can change the resolution ("-d") and the type of image created ("-i"), by using the following command:""","""make_rarefaction_plots.py -i chao1/ -d 180 -p pdf"""))
script_info['output_description']="""The result of this script produces a folder and within that folder there are sub-folders for each data file (metric) supplied as input. Within the sub-folders, there will be images for each of the categories specified by the user."""
script_info['required_options']=[\
make_option('-i', '--input_dir', help='name of folder containing rarefaction files, takes output from make_rarefaction_data.py [REQUIRED]')
]
script_info['optional_options']=[\
make_option('-t', '--rarefactionAve', help='name of overall average rarefaction file, takes output from make_rarefaction_data.py'),
make_option('-p', '--imagetype', type='string', help='extension for image type choose from (jpg, gif, png, svg, pdf). [default: %default]', default='png'),
make_option('-d', '--resolution', help='output image resolution, [default: %default]', type='int', default='75'),
make_option('-o', '--dir_path',help='directory prefix for all analyses [default: %default]', default='.'),
make_option('-y', '--ymax', help='maximum value for y axis, [default: %default] the default value will tell the script to calculate a y axis maximum depending on the data', type='int', default='0')
]
script_info['version'] = __version__

def main():
    option_parser, options, args = parse_command_line_parameters(**script_info)
      
    prefs = {}

    input_dir = options.input_dir
    rarenames = os.listdir(input_dir)
    rarenames = [r for r in rarenames if not r.startswith('.')]
    rares = dict()
    for r in rarenames:
        try:
             rarefl = open(input_dir + '/' + r, 'U').readlines()
             rares[r] = rarefl
        except(IOError):
            option_parser.error('Problem with rarefaction file. %s'%\
            sys.exc_info()[1])
            exit(0)
    prefs['rarefactions'] = rares
    
    if options.imagetype not in ['jpg','gif','png','svg','pdf']:
        option_parser.error('Supplied extension not supported.')
        exit(0)
    else:
        prefs['imagetype'] = options.imagetype
        
    try:
        prefs['resolution'] = int(options.resolution)
    except(ValueError):
        option_parser.error('Inavlid resolution.')
        exit(0)

    try:
        prefs['ymax'] = int(options.ymax)
    except(ValueError):
        option_parser.error('Inavlid maximum y axis value.')
        exit(0)

    #output directory check
    if options.dir_path != '.':
        if os.path.exists(options.dir_path):
            prefs['output_path'] = options.dir_path
        else:
            try:
                os.mkdir(options.dir_path)
                prefs['output_path'] = options.dir_path
            except(ValueError):
                option_parser.error('Could not create output directory.')
                exit(0)
    else:
        prefs['output_path'] = get_random_directory_name()
    
    graphNames = make_plots(prefs)
    make_output_files(prefs, util.get_qiime_project_dir(), graphNames)

if __name__ == "__main__":
    main()