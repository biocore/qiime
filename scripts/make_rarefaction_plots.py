#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Meg Pirrung"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Meg Pirrung", "Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.0.0-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Development"
 
from optparse import make_option
from qiime.util import parse_command_line_parameters, get_qiime_project_dir, \
                       create_dir,get_options_lookup
from cogent.util.misc import get_random_directory_name
from sys import argv, exit, exc_info
from qiime.colors import sample_color_prefs_and_map_data_from_options
from qiime.parse import parse_rarefaction_data,parse_rarefaction
from qiime.make_rarefaction_plots import make_averages
from os.path import exists, splitext, split,isdir
from os import listdir,path

options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""Generate Rarefaction Plots"""
script_info['script_description']="""Once the batch alpha diversity files have been collated, you may want to compare the diversity using plots. Using the results from collate_alpha.py, you can plot the samples and or by category in the mapping file using this script.

This script creates an html file of rarefaction plots based on the supplied collated alpha-diversity files in the folder given (-i) from collate_alpha.py. The user may also supply optional arguments like an image type (-i), and a resolution (-d)."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Default Example:""","""For generated rarefaction plots using the default parameters, including the mapping file and one rarefaction file, you can use the following command:""","""make_rarefaction_plots.py -r chao1/ -m mapping_file.txt"""))
script_info['script_usage'].append(("""Specify Image Type and Resolution:""","""Optionally, you can change the resolution ("-d") and the type of image created ("-i"), by using the following command:""","""make_rarefaction_plots.py -i chao1/ -m mapping_file.txt -d 180 -g pdf"""))
script_info['script_usage'].append(("""Use Prefs File:""","""You can also supply a preferences file "-p", as follows""","""make_rarefaction_plots.py -i chao1/ -m mapping_file.txt -d 180 -p prefs.txt"""))
script_info['script_usage'].append(("""Set Background Color:""","""Alternatively, you can set the plot background "-k", as follows: a preferences file "-p", as follows""","""make_rarefaction_plots.py -i chao1/ -m mapping_file.txt -k black"""))
script_info['output_description']="""The result of this script produces a folder and within that folder there is a sub-folder containing image files. Within the main folder, there is an html file."""
script_info['required_options']=[\
make_option('-i', '--input_dir', help='name of folder containing rarefaction files, takes output from collate_alpha.py.  The user can also supply a list of files, which are comma-separated. [REQUIRED]'),
make_option('-m', '--map_fname', help='name of mapping file [REQUIRED]')
]
script_info['optional_options']=[\
make_option('-b', '--colorby', type='string', help='name of columns to make rarefaction graphs of, comma delimited no spaces.'),
make_option('-p', '--prefs_path', type='string', help='preferences file for coloring of columns.'),
make_option('-k', '--background_color', type='string', help='Background color for graphs.',default='white'),
make_option('-g', '--imagetype', type='string', help='extension for image type choose from (png, svg, pdf).  WARNING: Some formats may not properly open in your browser! [default: %default]', default='png'),
make_option('-d', '--resolution', help='output image resolution, [default: %default]', type='int', default='75'),
make_option('-y', '--ymax', type='int', help='this is the ymax value to be used for the plots, so you can compare rarefaction plots between two different analyses [default: %default]'),
make_option('-w', '--webpage', action='store_false', help='this is allows to user to not create the webpage, which may be slow with large datasets [default: %default]', default=True),
options_lookup['output_dir']
]
script_info['version'] = __version__

#removed --rarefaction_average option and --ymax option, since it is calculated
# in the script



def main():
    option_parser, options, args = parse_command_line_parameters(**script_info)
      
    ops = {}
    input_dir = options.input_dir

    rares = {}
    if isdir(input_dir):
        rarenames = listdir(input_dir)
        rarenames = [r for r in rarenames if not r.startswith('.')]
        for r in rarenames:
            try:
                 rarefl = open(path.join(input_dir,r), 'U').readlines()
                 rares[r] = parse_rarefaction(rarefl)
            except(IOError):
                option_parser.error('Problem with rarefaction file. %s'%\
                exc_info()[1])
                exit(0)
    else:
        try:
             input_file=input_dir.split(',')
             for i in range(len(input_file)):
                 input_path=split(input_file[i])[-1]
                 rarefl = open(input_file[i], 'U').readlines()
                 rares[input_path] = parse_rarefaction(rarefl)
        except(IOError):
            option_parser.error('Problem with rarefaction file. %s'%\
            exc_info()[1])
            exit(0)
    if options.imagetype not in ['png','svg','pdf']:
        option_parser.error('Supplied extension not supported.')
        exit(0)
    else:
        imagetype = options.imagetype
        
    try:
        resolution = int(options.resolution)
    except(ValueError):
        option_parser.error('Inavlid resolution.')
        exit(0)
    
    #Get the command-line options.
    prefs, data, background_color, label_color = \
                    sample_color_prefs_and_map_data_from_options(options)
    
    #output directory check
    if options.output_dir != '.':
        if exists(options.output_dir):
            output_dir = options.output_dir
        else:
            try:
                create_dir(options.output_dir,False)
                output_dir = options.output_dir
            except(ValueError):
                option_parser.error('Could not create output directory.')
                exit(0)
    else:
        output_dir = get_random_directory_name()
    
    #Generate the plots and html text
    
    ymax=options.ymax
    make_webpage=options.webpage
    html_output = make_averages(prefs, data, background_color, label_color, \
                                rares, output_dir,resolution,imagetype,ymax,
                                make_webpage)

    if html_output:
        #Write the html file.
        outfile = open(path.join(output_dir,'rarefaction_plots.html'),'w')
        outfile.write(html_output)
        outfile.close()


if __name__ == "__main__":
    main()