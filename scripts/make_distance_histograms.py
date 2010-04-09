#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Jeremy Widmann","Rob Knight"]
__license__ = "GPL"
__version__ = "1.0.0"
__maintainer__ = "Jeremy Widmann"
__email__ = "Jeremy.Widmann@colorado.edu"
__status__ = "Release"
 

from qiime.util import parse_command_line_parameters, get_qiime_project_dir,\
    get_options_lookup
from optparse import make_option
from qiime.make_distance_histograms import group_distances, _make_path, \
    draw_all_histograms, _make_relative_paths, make_main_html, \
    monte_carlo_group_distances
from cogent.util.misc import get_random_directory_name
from qiime.colors import sample_color_prefs_and_map_data_from_options,\
    iter_color_groups
from os import mkdir
from string import strip

options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""Make distance histograms"""
script_info['script_description']="""To visualize the distance between samples and/or categories in the metadata mapping file, the user can generate histograms to represent the distances between samples. This script generates an HTML file, where the user can compare the distances between samples based on the different categories associated to each sample in the metadata mapping file. """
script_info['script_usage']=[]
script_info['script_usage'].append(("""Examples:""","""Distance Histograms are a way to compare different categories and see which tend to have larger/smaller distances than others. For example, in the hand study, you may want to compare the distances between hands to the distances between individuals (with the file "hand_distances.txt" using the parameter -d hand_distances.txt). The categories are defined in the metadata mapping file (specified using the parameter -m hand_map.txt). If you want to look at the distances between hands and individuals, choose the "Hand" field and "Individual" field (using the parameter --fields Hand,Individual (notice the fields are comma delimited)). For each of these groups of distances a histogram is made. The output is a HTML file ("QIIME_Distance_Histograms.html" when the parameter --html_output is specified) which is created in the "Distance_Histograms" directory (using the parameter -o Distance_Histograms to specify output directory) where you can look at all the distance histograms individually, and compare them between each other.

In the following command, the user only supplies a distance matrix (i.e. resulting file from beta_diversity.py), the user-generated metadata mapping file and one category (e.g. pH):""","""make_distance_histograms.py -d beta_div.txt -m Mapping_file.txt --fields pH"""))
script_info['script_usage'].append(("""""","""For comparison of multiple categories (e.g. pH, salinity), you can use the following command:""","""make_distance_histograms.py -d beta_div.txt -m Mapping_file.txt --fields pH,salinity"""))
script_info['script_usage'].append(("""""","""If the user would like to write the result to a dynamic HTML, you can use the following command:""","""make_distance_histograms.py -d beta_div.txt -m Mapping_file.txt --fields pH --html_output"""))
script_info['script_usage'].append(("""""","""In the case that the user generates their own preferences file (prefs.txt), they can use the following command:""","""make_distance_histograms.py -d beta_div.txt -m Mapping_file.txt -p prefs.txt"""))
script_info['script_usage'].append(("""""","""Note: In the case that a preferences file is passed, the user does not need to supply fields in the command-line.""",""""""))
script_info['output_description']="""The result of this script will be a folder containing images and/or an html file (with appropriate javascript files), depending on the user-defined parameters."""
script_info['required_options']=[\
    make_option('-d','--distance_matrix_file',dest='distance_matrix_file',\
        type='string',help='''Path to distance matrix file.'''),\
    make_option('-m', '--map_fname', dest='map_fname', \
         help='This is the metadata mapping file  [default=%default]'), \
]

script_info['optional_options']=[\
    make_option('-p', '--prefs_path',help='This is the user-generated preferences \
file. NOTE: This is a file with a dictionary containing preferences for the \
analysis.  This dict must have a "Fields" key mapping to a list of desired fields. [default: %default]'),
    make_option('-o', '--dir_path', dest='dir_path', default='.',\
        help='Directory to output data for all analyses. [default: %default]'
),\
    make_option('-k', '--background_color', default='white', help='This is the \
    background color to use in the plots (Options are \'black\' or \'white\'. \
    [default: %default]'),\
    make_option('--monte_carlo',dest='monte_carlo',default=False,\
        action='store_true',help='''Perform Monte Carlo analysis on distances.  [Default: %default]'''),\
    make_option('--html_output',dest='html_output',default=False,\
        action='store_true',help='''Write output in HTML format. [Default: %default]'''),\
    make_option('-f','--fields', dest='fields',\
        help='Comma delimited list of fields to compare.  This overwrites fields in prefs file.  If this is not provided, the first field in metadata mapping file will be used.  Usage: --fields Field1,Field2,Field3'),\
    make_option('--monte_carlo_iters', dest='monte_carlo_iters',type="int",\
        default=10,help='Number of iterations to perform for Monte Carlo analysis. [default: %default]'),\

]

script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    color_prefs, color_data, background_color, label_color= \
                            sample_color_prefs_and_map_data_from_options(opts)
    
    #list of labelname, groups, colors, data_colors, data_color_order    
    groups_and_colors=list(iter_color_groups(mapping=color_data['map'],\
        prefs=color_prefs))
    
    #dict mapping labelname to list of: [groups, colors, data_colors,
    # data_color_order]
    field_to_colors = {}
    for color_info in groups_and_colors:
        field_to_colors[color_info[0]]=color_info[1:]
    
    qiime_dir = get_qiime_project_dir()+'/qiime/support_files/'
    
    if opts.prefs_path:
        prefs = eval(open(opts.prefs_path, 'U').read())
    else:
        prefs=None
    
    fields = opts.fields
    if fields is not None:
        fields = map(strip,fields.split(','))
    elif opts.prefs_path is not None:
        prefs = eval(open(opts.prefs_path, 'U').read())
        fields = prefs.get('FIELDS',None)
    
    within_distances, between_distances, dmat = \
        group_distances(mapping_file=opts.map_fname,\
        dmatrix_file=opts.distance_matrix_file,\
        fields=fields,\
        dir_prefix=get_random_directory_name(output_dir=opts.dir_path,\
            prefix='distances'))
    
    if opts.html_output:
        #histograms output path
        histograms_path = \
            _make_path([opts.dir_path,'histograms'])
        try:
            mkdir(histograms_path)
        except OSError:     #raised if dir exists
            pass
        
        #draw all histograms
        distances_dict, label_to_histogram_filename = \
            draw_all_histograms(single_field=within_distances, \
                paired_field=between_distances, \
                dmat=dmat,\
                histogram_dir=histograms_path,\
                field_to_color_prefs=field_to_colors,\
                background_color=background_color)
        
        #Get relative path to histogram files.
        label_to_histogram_filename_relative = \
            _make_relative_paths(label_to_histogram_filename, opts.dir_path)
        
        outfile_name = 'QIIME_Distance_Histograms.html'
        make_main_html(distances_dict=distances_dict,\
            label_to_histogram_filename=label_to_histogram_filename_relative,\
            root_outdir=opts.dir_path, \
            outfile_name = outfile_name, \
            title='Distance Histograms')
        
        #Handle saving web resources locally.
        #javascript file
        javascript_path = \
            _make_path([opts.dir_path,'js'])
        try:
            mkdir(javascript_path)
        except OSError:     #raised if dir exists
            pass
        js_out = open(javascript_path+'/histograms.js','w')
        js_out.write(open(qiime_dir+'js/histograms.js').read())
        js_out.close()
        
        #Qiime logo
        logo_path = \
            _make_path([opts.dir_path,'web_resources'])
        try:
            mkdir(logo_path)
        except OSError:     #raised if dir exists
            pass
        logo_out = open(logo_path+'/qiime_header.png','w')
        logo_out.write(open(qiime_dir+'images/qiime_header.png').read())
        logo_out.close()
    
    if opts.monte_carlo:
        monte_carlo_group_distances(mapping_file=opts.map_fname,\
            dmatrix_file=opts.distance_matrix_file,\
            prefs=prefs, \
            dir_prefix = opts.dir_path,\
            fields=fields,\
            default_iters=opts.monte_carlo_iters)


if __name__ == "__main__":
    main()