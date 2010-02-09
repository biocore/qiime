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
from qiime.make_rarefaction_plots import make_plots, make_output_files, is_max_category_ops
import os.path
from os.path import exists, splitext, split
import shutil

script_description = """Create an html file of rarefaction plots based on the supplied\
mapping file (-m) and the supplied rarefaction files (-r) from beta_diversity.py. \
The user may also supply optional arguments that will only create plots for \
supplied metadata columns from the mapping file (-p), an image type (-i), and a\
resolution (-d)"""

script_usage = """Usage: %prog [options] {-o OUTPUT_FP}

Example Usage:
Create an html file with all rarefaction plots for all relevant categories:
python %prog -m mappingfile.txt -r rare1.txt,rare2.txt

Create an html file with rarefaction plots for only the provided categories:
python %prog -m mappingfile.txt -r rare1.txt,rare2.txt -p DAY,DONOR

Create an html file with rarefaction plots for provided categories, images of type PNG,\
resolution of 150:
python %prog -m mappingfile.txt -r rare1.txt,rare2.txt -p DAY,DONOR -i png -d 150
"""

required_options = [\
make_option('-m', '--map', help='name of mapping file [REQUIRED]'),
make_option('-r', '--rarefaction', help='name of rarefaction file [REQUIRED]')
]

optional_options = [\
make_option('-p', '--prefs', type='string', help='name of columns to make rarefaction graphs of, \
comma delimited no spaces. Use \'ALL\' command to make graphs of all metadata columns. \
[default: %default]', default='ALL'),
make_option('-i', '--imagetype', type='string', help='extension for image type choose from \
(jpg, gif, png, svg, pdf). [default: %default]', default='png'),
make_option('-d', '--resolution', help='output image resolution, \
[default: %default]', type='int', default='75'),
make_option('-o', '--dir_path',help='directory prefix for all analyses \
[default: %default]', default='.')
]

def main():
    option_parser, options, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)
      
    prefs = {}
    
    try:
        prefs['mapfl'] = open(options.map, 'U').readlines()
    except(IOError):
        option_parser.error('Problem with mapping file. %s'%sys.exc_info()[1])
        exit(0)
    prefs['map'] = parse.parse_map(prefs['mapfl'], return_header=True, strip_quotes=True)
    prefs['map'][0][0] = [h.strip('#').strip(' ') for h in prefs['map'][0][0]]

    rarenames = options.rarefaction.split(',')
    rares = dict()
    for r in rarenames:
        try:
            rares[r] = open(r, 'U').readlines()
        except(IOError):
            option_parser.error('Problem with rarefaction file. %s'%sys.exc_info()[1])
            exit(0)
    prefs['rarefactions'] = rares

    if options.prefs.split(',')[0] == 'ALL':
        prefs['categories'] = []
        temp = prefs['map'][0][0]
        for p in temp:
            is_max, l = is_max_category_ops(prefs['map'], p)
            if l != 1 and not is_max:
                prefs['categories'].append(p)
    else:
        suppliedcats =  set(options.prefs.split(','))
        availablecats = set(prefs['map'][0][0])
        
        if suppliedcats.issubset(availablecats):
            prefs['categories'] = options.prefs.split(',')
        else:
            option_parser.error('Categories %s not found in mapping file, \
please check spelling and syntax.'%list(suppliedcats.difference(availablecats)))
            exit(0)

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

    if '/' in argv[0]:
        prefs['output_path'] = argv[0].rsplit('/',1)[0]+'/'
    else:
        prefs['output_path'] = './'
    
    dir_path = options.dir_path
    dir_path = os.path.join(dir_path,'rarefaction_graphs')

    alphabet = "ABCDEFGHIJKLMNOPQRSTUZWXYZ"
    alphabet += alphabet.lower()
    alphabet += "01234567890"
    data_file_path=''.join([choice(alphabet) for i in range(10)])
    prefs['output_path'] = os.path.join(dir_path,data_file_path)
    
    outputlines = make_plots(prefs)
    make_output_files(prefs, outputlines, util.get_qiime_project_dir())

if __name__ == "__main__":
    main()