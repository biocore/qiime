#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jeremy Widmann","Rob Knight"]
__license__ = "GPL"
__version__ = "1.2.1-dev"
__maintainer__ = "Jeremy Widmann"
__email__ = "Jeremy.Widmann@colorado.edu"
__status__ = "Development"
 

from qiime.util import parse_command_line_parameters, get_qiime_project_dir,\
    get_options_lookup
from qiime.util import make_option
from qiime.make_distance_histograms import group_distances, _make_path, \
    draw_all_histograms, _make_relative_paths, make_main_html, \
    monte_carlo_group_distances, monte_carlo_group_distances_within_between
from cogent.util.misc import get_random_directory_name
from qiime.colors import sample_color_prefs_and_map_data_from_options,\
    iter_color_groups
from qiime.parse import (parse_mapping_file, parse_distmat, 
    parse_prefs_file, QiimeParseError)
from os import mkdir
from string import strip

options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""Make distance histograms"""
script_info['script_description']="""To visualize the distance between samples and/or categories in the metadata mapping file, the user can generate histograms to represent the distances between samples. This script generates an HTML file, where the user can compare the distances between samples based on the different categories associated to each sample in the metadata mapping file. """
script_info['script_usage']=[]
script_info['script_usage'].append(("""Examples:""","""Distance Histograms are a way to compare different categories and see which tend to have larger/smaller distances than others. For example, in the hand study, you may want to compare the distances between hands to the distances between individuals (with the file "hand_distances.txt" using the parameter -d hand_distances.txt). The categories are defined in the metadata mapping file (specified using the parameter -m hand_map.txt). If you want to look at the distances between hands and individuals, choose the "Hand" field and "Individual" field (using the parameter --fields Hand,Individual (notice the fields are comma delimited)). For each of these groups of distances a histogram is made. The output is a HTML file ("QIIME_Distance_Histograms.html") which is created in the "Distance_Histograms" directory (using the parameter -o Distance_Histograms to specify output directory) where you can look at all the distance histograms individually, and compare them between each other.

In the following command, the user only supplies a distance matrix (i.e. resulting file from beta_diversity.py), the user-generated metadata mapping file and one category (e.g. pH):""","""make_distance_histograms.py -d beta_div.txt -m Mapping_file.txt --fields pH"""))
script_info['script_usage'].append(("""""","""For comparison of multiple categories (e.g. pH, salinity), you can use the following command:""","""make_distance_histograms.py -d beta_div.txt -m Mapping_file.txt --fields pH,salinity"""))
script_info['script_usage'].append(("""""","""HTML output is automatically generated. If the user would like to suppress the HTML output, you can use the following command:""","""make_distance_histograms.py -d beta_div.txt -m Mapping_file.txt --fields pH --suppress_html_output"""))
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
    make_option('-o', '--dir_path', dest='dir_path', type='string', \
        default='.',\
        help='Directory to output data for all analyses. [default: %default]'
),\
    make_option('-k', '--background_color', dest='background_color',\
        default='white', help='This is the \
    background color to use in the plots (Options are \'black\' or \'white\'. \
    [default: %default]'),\
    make_option('--monte_carlo',dest='monte_carlo',default=False,\
        action='store_true',help='''Perform Monte Carlo analysis on distances.  [Default: %default]'''),\
    make_option('--suppress_html_output',dest='suppress_html_output',\
        default=False,action='store_true',help='''Suppress HTML format output. [Default: %default]'''),\
    make_option('-f','--fields', dest='fields',\
        help='Comma delimited list of fields to compare.  Put list of fields in quotes.  This overwrites fields in prefs file.  If this is not provided, the first field in metadata mapping file will be used.  Usage: --fields "Field1,Field2,Field3"'),\
    make_option('--monte_carlo_iters', dest='monte_carlo_iters',type="int",\
        default=100,help='Number of iterations to perform for Monte Carlo analysis. [default: %default]'),\

]

script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    #Some code for error checking of input args:
    
    #Check if distance_matrix_file is valid:
    try:
        d_header, d_mat = parse_distmat(open(opts.distance_matrix_file,'U'))
    except:
        option_parser.error("This does not look like a valid distance matrix file.  Please supply a valid distance matrix file using the -d option.")
    
    #Check if map_fname is valid:
    try:
        mapping, m_header, m_comments = \
            parse_mapping_file(open(opts.map_fname,'U'))
    except QiimeParseError:
        option_parser.error("This does not look like a valid metadata mapping file.  Please supply a valid mapping file using the -m option.")
    
    #make sure background_color is valid
    if opts.background_color not in ['black','white']:
        option_parser.error("'%s' is not a valid background color.  Please pass in either 'black' or 'white' using the -k option."%(opts.background_color))
    
    #make sure prefs file is valid if it exists
    if opts.prefs_path is not None:
        try:
            prefs_file = open(opts.prefs_path, 'U').read()
        except IOError:
            option_parser.error("Provided prefs file, '%s', does not exist.  Please pass in a valid prefs file with the -p option."%(opts.prefs_path))
            
    if opts.prefs_path is not None:
        prefs = parse_prefs_file(prefs_file)
    else:
        prefs=None

    
    color_prefs, color_data, background_color, label_color, ball_scale, arrow_colors= \
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
    
        
    fields = opts.fields
    if fields is not None:
        fields = map(strip,fields.split(','))
        fields = [i.strip('"').strip("'") for i in fields]
    elif prefs is not None:
        fields = prefs.get('FIELDS',None)
    
    #Check that all provided fields are valid:
    if fields is not None:
        for f in fields:
            if f not in m_header:
                option_parser.error("The field, %s, is not in the provided mapping file.  Please supply correct fields (using the -f option or providing a 'FIELDS' list in the prefs file) corresponding to fields in mapping file."%(f))
    
    within_distances, between_distances, dmat = \
        group_distances(mapping_file=opts.map_fname,\
        dmatrix_file=opts.distance_matrix_file,\
        fields=fields,\
        dir_prefix=get_random_directory_name(output_dir=opts.dir_path,\
            prefix='distances'))
    
    if not opts.suppress_html_output:
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
        #Do Monte Carlo for all fields
        monte_carlo_group_distances(mapping_file=opts.map_fname,\
            dmatrix_file=opts.distance_matrix_file,\
            prefs=prefs, \
            dir_prefix = opts.dir_path,\
            fields=fields,\
            default_iters=opts.monte_carlo_iters)
            
        #Do Monte Carlo for within and between fields
        monte_carlo_group_distances_within_between(\
            single_field=within_distances,\
            paired_field=between_distances, dmat=dmat, \
            dir_prefix = opts.dir_path,\
            num_iters=opts.monte_carlo_iters)


if __name__ == "__main__":
    main()
