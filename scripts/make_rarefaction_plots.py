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
script_info['script_description']="""Once the batch alpha diversity files have been collated, you may want to compare the diversity using plots. Using the results from collate_alpha.py, you can plot the samples and or by category in the mapping file using this script.

This script creates an html file of rarefaction plots based on the supplied mapping file (-m) and the supplied rarefaction files (-r) from collate_alpha.py. The user may also supply optional arguments that will only create plots for supplied metadata columns from the mapping file (-p), an image type (-i), and a resolution (-d). If the user would like to suppress html output they can pass the -n flag, and output raw data with the -w flag. The -y option allows the user to supply a maximum value for the yaxis of the plots."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Default Example:""","""For generated rarefaction plots using the default parameters, including the mapping file and one rarefaction file, you can use the following command:""","""make_rarefaction_plots.py -m Mapping_file.txt -r chao1.txt"""))
script_info['script_usage'].append(("""Multiple File Example:""","""If you would like to generate plots for multiple files, you can use the following command:""","""make_rarefaction_plots.py -m Mapping_file.txt -r chao1.txt,PD_whole_tree.txt"""))
script_info['script_usage'].append(("""Category Specific Example:""","""In the case that you want to make plots for a specific category (i.e., pH), you can use the following command:""","""make_rarefaction_plots.py -m Mapping_file.txt -r chao1.txt -p"""))
script_info['script_usage'].append(("""Specify Image Type and Resolution:""","""Optionally, you can change the resolution ("-d") and the type of image created ("-i"), by using the following command:""","""make_rarefaction_plots.py -m Mapping_file.txt -r chao1.txt -p pH -d 180 -i pdf"""))
script_info['output_description']="""The result of this script produces a folder and within that folder there are sub-folders for each data file (metric) supplied as input. Within the sub-folders, there will be images for each of the categories specified by the user."""
script_info['required_options']=[]
script_info['optional_options']=[]

script_description = """Create an html file of rarefaction plots based on the supplied\
the supplied rarefaction files (-r) from make_rarefaction_data.py. \
The user may also supply optional arguments that will only create plots for \
supplied metadata columns from the mapping file (-p), an image type (-i), and a\
resolution (-d). The -y option allows the\
user to supply a maximum value for the yaxis of the plots."""

script_usage = """Usage: %prog [options] {-r RAREFACTION}

Example Usage:
Create an html file with all rarefaction plots for all relevant categories:
python %prog -r rare1.txt,rare2.txt

Create an html file with rarefaction plots for provided categories, images of type PNG, resolution of 150:
python %prog -r rare1.txt,rare2.txt -i png -d 150
"""

required_options = [\
make_option('-i', '--input_dir', help='name of folder containing rarefaction files, takes output from make_rarefaction_data.py [REQUIRED]')
]

optional_options = [\
make_option('-t', '--rarefactionAve', help='name of overall average rarefaction file, takes output from make_rarefaction_data.py'),
make_option('-p', '--imagetype', type='string', help='extension for image type choose from (jpg, gif, png, svg, pdf). [default: %default]', default='png'),
make_option('-d', '--resolution', help='output image resolution, [default: %default]', type='int', default='75'),
make_option('-o', '--dir_path',help='directory prefix for all analyses [default: %default]', default='.'),
make_option('-y', '--ymax', help='maximum value for y axis, [default: %default] the default value will tell the script to calculate a y axis maximum depending on the data', type='int', default='0')
]

def main():
    option_parser, options, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)
      
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
            option_parser.error('Problem with rarefaction file. %s'%sys.exc_info()[1])
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

    if '/' in argv[0]:
        prefs['output_path'] = argv[0].rsplit('/',1)[0]+'/'
    else:
        prefs['output_path'] = './'
    
    dir_path = options.dir_path
    dir_path = os.path.join(dir_path,'rarefaction_plots')

    alphabet = "ABCDEFGHIJKLMNOPQRSTUZWXYZ"
    alphabet += alphabet.lower()
    alphabet += "01234567890"
    data_file_path=''.join([choice(alphabet) for i in range(10)])
    prefs['output_path'] = os.path.join(dir_path,data_file_path)
    
    make_plots(prefs)
    make_output_files(prefs, util.get_qiime_project_dir())

if __name__ == "__main__":
    main()