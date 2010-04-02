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
from qiime.parse import parse_rarefaction
from cogent.util.misc import get_random_directory_name
import sys
from sys import argv, exit, exc_info
from random import choice, randrange
from time import strftime
from qiime.colors import sample_color_prefs_and_map_data_from_options
from qiime.parse import parse_mapping_file, parse_rarefaction
from qiime.make_rarefaction_averages import make_averages, \
is_max_category_ops, parse_rarefaction
import os.path
from os.path import exists, splitext, split
import shutil


#make_rarefaction_averages.py
script_info={}
script_info['brief_description']="""Generate Rarefaction Averages"""
script_info['script_description']="""Once the batch alpha diversity files have been collated, you may want to compare the diversity. Using the results from collate_alpha.py, you can average rarefaction values across sample metadata in the mapping file using this script.

This script creates a directory of average rarefaction series based on the supplied mapping file (-m) and the supplied rarefaction files (-r) from collate_alpha.py."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Default Example:""","""For generated rarefaction plots using the default parameters, including the mapping file and one rarefaction file, you can use the following command:""","""make_rarefaction_plots.py -m Mapping_file.txt -r chao1.txt"""))
script_info['script_usage'].append(("""Multiple File Example:""","""If you would like to generate plots for multiple files, you can use the following command:""","""make_rarefaction_plots.py -m Mapping_file.txt -r chao1.txt,PD_whole_tree.txt"""))
script_info['script_usage'].append(("""Category Specific Example:""","""In the case that you want to make plots for a specific category (i.e., pH), you can use the following command:""","""make_rarefaction_plots.py -m Mapping_file.txt -r chao1.txt -p pH"""))
script_info['output_description']="""The result of this script produces a folder and within that folder there are sub-folders for each data file (metric) supplied as input. Within the sub-folders, there will be text files of averages for each of the categories specified by the user."""

script_info['required_options']=[\
make_option('-m', '--map_fname', help='name of mapping file [REQUIRED]'),
make_option('-r', '--rarefaction', help='name of rarefaction file, takes output from collate_alpha OR tab delimited data from a previous run of this script. If using raw data from a previous run, set -x flag. [REQUIRED]')
]
script_info['optional_options']=[\
make_option('-k', '--background_color', type='string', help='Background color for graphs.'),
make_option('-b', '--colorby', type='string', help='name of columns to make rarefaction graphs of, comma delimited no spaces.'),
make_option('-p', '--prefs_path', type='string', help='preferences file for coloring of columns.'),
make_option('-o', '--dir_path',help='directory prefix for all analyses [default: %default]', default='.')
]
script_info['version'] = __version__



def main():
    option_parser, options, args = parse_command_line_parameters(**script_info)
      
    ops = {}

    #rarefaction data check
    rarenames = options.rarefaction.split(',')
    rares = dict()
    for r in rarenames:
        try:
             rarefl = open(r, 'U').readlines()
             rares[r] = parse_rarefaction(rarefl)
        except(IOError):
            option_parser.error('Problem with rarefaction file. %s'%\
            sys.exc_info()[1])
            exit(0)
    ops['rarefactions'] = rares

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
    # #colorby check
    # if options.colorby:
    #     if options.colorby.split(',')[0] == 'ALL':
    #         ops['categories'] = []
    #         temp = ops['map'][1]
    #         for p in temp:
    #             is_max, l = is_max_category_ops(ops['map'], p)
    #             if l != 1 and not is_max:
    #                 ops['categories'].append(p)
    #     else:
    #         suppliedcats =  set(options.ops.split(','))
    #         availablecats = set(ops['map'][1])
    #     
    #         if suppliedcats.issubset(availablecats):
    #             ops['categories'] = options.ops.split(',')
    #         else:
    #             option_parser.error('Categories %s not found in mapping file, \
    # please check spelling and syntax.'%list(suppliedcats.difference(availablecats)))
    #             exit(0)
    
    #output directory check
    if options.dir_path != '.':
        if os.path.exists(options.dir_path):
            ops['output_path'] = options.dir_path
        else:
            try:
                os.mkdir(options.dir_path)
                ops['output_path'] = options.dir_path
            except(ValueError):
                option_parser.error('Could not create output directory.')
                exit(0)
    else:
        ops['output_path'] = get_random_directory_name()
    
    make_averages(ops)

if __name__ == "__main__":
    main()