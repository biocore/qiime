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
from qiime.util import parse_command_line_parameters, get_qiime_project_dir
from cogent.util.misc import get_random_directory_name
from sys import argv, exit, exc_info
from random import choice, randrange
from time import strftime
from qiime.colors import sample_color_prefs_and_map_data_from_options
from qiime.parse import parse_rarefaction_data
from qiime.make_rarefaction_plots import make_plots, make_output_files
from os.path import exists, splitext, split
from os import listdir, mkdir

script_info={}
script_info['brief_description']="""Generate Rarefaction Plots"""
script_info['script_description']="""Once the batch alpha diversity files have been collated, you may want to compare the diversity using plots. Using the results from make_rarefaction_averages.py, you can plot the samples and or by category in the mapping file using this script.

This script creates an html file of rarefaction plots based on the supplied rarefaction files in the folder given (-i) from make_rarefaction_averages.py. The user may also supply optional arguments like an image type (-i), and a resolution (-d)."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Default Example:""","""For generated rarefaction plots using the default parameters, including the mapping file and one rarefaction file, you can use the following command:""","""make_rarefaction_plots.py -r chao1/"""))
script_info['script_usage'].append(("""Specify Image Type and Resolution:""","""Optionally, you can change the resolution ("-d") and the type of image created ("-i"), by using the following command:""","""make_rarefaction_plots.py -i chao1/ -d 180 -p pdf"""))
script_info['output_description']="""The result of this script produces a folder and within that folder there are sub-folders for each data file (metric) supplied as input. Within the sub-folders, there will be images for each of the categories specified by the user."""
script_info['required_options']=[\
make_option('-i', '--input_dir', help='name of folder containing rarefaction files, takes output from collate_alpha.py [REQUIRED]'),
make_option('-m', '--map_fname', help='name of mapping file [REQUIRED]')
]
script_info['optional_options']=[\
make_option('-t', '--rarefactionAve', help='name of overall average rarefaction file, takes output from make_rarefaction_averages.py'),
make_option('-b', '--colorby', type='string', help='name of columns to make rarefaction graphs of, comma delimited no spaces.'),
make_option('-p', '--prefs_path', type='string', help='preferences file for coloring of columns.'),
make_option('-k', '--background_color', type='string', help='Background color for graphs.'),
make_option('-g', '--imagetype', type='string', help='extension for image type choose from (jpg, gif, png, svg, pdf). [default: %default]', default='png'),
make_option('-d', '--resolution', help='output image resolution, [default: %default]', type='int', default='75'),
make_option('-o', '--dir_path',help='directory prefix for all analyses [default: %default]', default='.'),
make_option('-y', '--ymax', help='maximum value for y axis, [default: %default] the default value will tell the script to calculate a y axis maximum depending on the data', type='int', default='0')
]
script_info['version'] = __version__

def main():
    option_parser, options, args = parse_command_line_parameters(**script_info)
      
    ops = {}

    input_dir = options.input_dir
    rarenames = listdir(input_dir)
    rarenames = [r for r in rarenames if not r.startswith('.')]
    rares = dict()
    for r in rarenames:
        try:
             rarefl = open(input_dir + '/' + r, 'U').readlines()
             rares[r] = parse_rarefaction_data(rarefl)
        except(IOError):
            option_parser.error('Problem with rarefaction file. %s'%\
            exc_info()[1])
            exit(0)
    ops['rarefactions'] = rares
    
    if options.imagetype not in ['jpg','gif','png','svg','pdf']:
        option_parser.error('Supplied extension not supported.')
        exit(0)
    else:
        ops['imagetype'] = options.imagetype
        
    try:
        ops['resolution'] = int(options.resolution)
    except(ValueError):
        option_parser.error('Inavlid resolution.')
        exit(0)

    try:
        ops['ymax'] = int(options.ymax)
    except(ValueError):
        option_parser.error('Inavlid maximum y axis value.')
        exit(0)

    #prefs check
    if(options.prefs_path):
        try:
            open(options.prefs_path, 'U').readlines()
        except(IOError):
            option_parser.error('Problem with prefs file. %s'%\
            sys.exc_info()[1])
            exit(0)
    ops['prefs_path'] = options.prefs_path
    ops['prefs'] = sample_color_prefs_and_map_data_from_options(options)

    ops['colorby'] = options.colorby

    #output directory check
    if options.dir_path != '.':
        if exists(options.dir_path):
            ops['output_path'] = options.dir_path
        else:
            try:
                mkdir(options.dir_path)
                ops['output_path'] = options.dir_path
            except(ValueError):
                option_parser.error('Could not create output directory.')
                exit(0)
    else:
        ops['output_path'] = get_random_directory_name()
    
    graphNames = make_plots(ops)
    make_output_files(ops, get_qiime_project_dir(), graphNames)

if __name__ == "__main__":
    main()