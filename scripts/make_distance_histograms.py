#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Jeremy Widmann","Rob Knight","Jesse Stombaugh",
               "Jai Ram Rideout"]
__license__ = "GPL"
__version__ = "1.8.0"
__maintainer__ = "Jeremy Widmann"
__email__ = "Jeremy.Widmann@colorado.edu"
__status__ = "Development"
 

from qiime.util import parse_command_line_parameters, get_qiime_project_dir,\
    get_options_lookup
from qiime.util import make_option,get_interesting_mapping_fields
from qiime.make_distance_histograms import group_distances, \
    draw_all_histograms, _make_relative_paths, make_main_html, \
    monte_carlo_group_distances, monte_carlo_group_distances_within_between
from cogent.util.misc import get_random_directory_name
from qiime.colors import sample_color_prefs_and_map_data_from_options,\
    iter_color_groups
from qiime.parse import (parse_mapping_file, parse_distmat, 
    parse_prefs_file, QiimeParseError)
from qiime.pycogent_backports.test import is_symmetric_and_hollow
from os import mkdir,path
from string import strip

options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""Make distance histograms"""
script_info['script_description']="""
To visualize the distance between samples and/or categories in the metadata
mapping file, the user can generate histograms to represent the distances
between samples. This script generates an HTML file, where the user can compare
the distances between samples based on the different categories associated to
each sample in the metadata mapping file.

Distance histograms provide a way to compare different categories and see which
tend to have larger/smaller distances than others. For example, in a hand
study, you may want to compare the distances between hands to the distances
between individuals (with the file "hand_distances.txt" using the parameter -d
hand_distances.txt). The categories are defined in the metadata mapping file
(specified using the parameter -m hand_map.txt). If you want to look at the
distances between hands and individuals, choose the "Hand" field and
"Individual" field (using the parameter --fields Hand,Individual (notice the
fields are comma-delimited)). For each of these groups of distances a
histogram is made. The output is an HTML file which is created in the
"Distance_Histograms" directory (using the parameter -o Distance_Histograms to
specify output directory) where you can look at all the distance histograms
individually, and compare them between each other.
"""

script_info['script_usage']=[]
script_info['script_usage'].append(("Distance histograms example",
"In the following command, the user supplies a distance matrix (i.e. the "
"resulting file from beta_diversity.py), the user-generated metadata mapping "
"file and one category \"Treatment\" to generate distance histograms.",
"%prog -d unweighted_unifrac_dm.txt -m Fasting_Map.txt --fields Treatment "
"-o example1"))
script_info['script_usage'].append(("Multiple categories",
"For comparison of multiple categories (e.g. Treatment, DOB), you can use the "
"following command (separating each category with a comma).",
"%prog -d unweighted_unifrac_dm.txt -m Fasting_Map.txt --fields "
"Treatment,DOB -o example2"))
script_info['script_usage'].append(("Suppress HTML output",
"By default, HTML output is automatically generated. If the user would like "
"to suppress the HTML output, you can use the following command.",
"%prog -d unweighted_unifrac_dm.txt -m Fasting_Map.txt --fields Treatment "
"--suppress_html_output -o example3"))
script_info['script_usage'].append(("Preferences file",
"You can provide your own preferences file (prefs.txt) with the following "
"command. If a preferences file is supplied, you do not need to supply fields "
"on the command-line.",
"%prog -d unweighted_unifrac_dm.txt -m Fasting_Map.txt -p prefs.txt -o "
"example4"))

script_info['output_description']="""
The result of this script will be a folder containing images and/or an HTML
file (with appropriate javascript files), depending on the user-defined
parameters.
"""

script_info['required_options']=[
    make_option('-d','--distance_matrix_file',
        help='Input distance matrix filepath (i.e. the result of '
        'beta_diversity.py). WARNING: Only symmetric, hollow distance '
        'matrices may be used as input. Asymmetric distance matrices, such as '
        'those obtained by the UniFrac Gain metric (i.e. beta_diversity.py '
        '-m unifrac_g), should not be used as input',
        type='existing_filepath'),
    make_option('-m', '--map_fname', dest='map_fname',
         help='Input metadata mapping filepath.',
         type='existing_filepath')
]

script_info['optional_options']=[
    make_option('-p', '--prefs_path',
        help='Input user-generated preferences filepath. NOTE: This is a '
        'file with a dictionary containing preferences for the analysis. '
        'This dictionary must have a "Fields" key mapping to a list of '
        'desired fields. [default: %default]',
        type='existing_filepath'),
    make_option('-o', '--dir_path',
        default='./',help='Output directory. [default: %default]',
        type='new_dirpath'),
    make_option('-k', '--background_color', dest='background_color',
        default='white', type='choice',choices=['black','white'],
        help='Background color for use in the plots '
        '(black or white) [default: %default]'),
    make_option('--monte_carlo',dest='monte_carlo',default=None,
        action='store_true',
        help='Deprecated: pass --monte_carlo_iters > 0 to enable'),
    make_option('--suppress_html_output',dest='suppress_html_output',
        default=False,action='store_true',
        help='Suppress HTML output. [default: %default]'),
    make_option('-f','--fields', default=None, type='string',
        help='Comma-separated list of fields to compare, where the list of '
        'fields should be in quotes (e.g. "Field1,Field2,Field3"). '
        'Note: if this option is passed on the command-line, it will '
        'overwrite the fields in prefs file. [default: first field in mapping '
        'file is used]'),
    make_option('--monte_carlo_iters', dest='monte_carlo_iters',type="int",
        default=0,
        help='Number of iterations to perform for Monte Carlo analysis. '
        '[default: %default; No monte carlo simulation performed]')
]
script_info['option_label']={'distance_matrix_file':'Distance matrix filepath',
                             'map_fname':'QIIME-formatted mapping filepath',
                             'prefs_path': 'Preferences filepath',
                             'dir_path': 'Output directory',
                             'background_color': 'Background color',
                             'monte_carlo': 'Perform Monte Carlo',
                             'monte_carlo_iters':'# of Monte Carlo iterations',
                             'suppress_html_output': 'Suppress HTML',
                             'fields':'Categories to compare'}

script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    #Some code for error checking of input args:
    
    #Check if distance_matrix_file is valid:
    try:
        d_header, d_mat = parse_distmat(open(opts.distance_matrix_file,'U'))
    except:
        option_parser.error("This does not look like a valid distance matrix file.  Please supply a valid distance matrix file using the -d option.")

    if not is_symmetric_and_hollow(d_mat):
        option_parser.error("The distance matrix must be symmetric and "
                            "hollow.")

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

    
    color_prefs, color_data, background_color, label_color, ball_scale,\
     arrow_colors=sample_color_prefs_and_map_data_from_options(opts)
    
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
    else:
        fields = get_interesting_mapping_fields(mapping, m_header)
    
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
        histograms_path = path.join(opts.dir_path,'histograms')
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
        
        dm_fname=path.split(opts.distance_matrix_file)[-1]
        basename=path.splitext(dm_fname)[0]
        outfile_name = basename+'_distance_histograms.html'
        make_main_html(distances_dict=distances_dict,\
            label_to_histogram_filename=label_to_histogram_filename_relative,\
            root_outdir=opts.dir_path, \
            outfile_name = outfile_name, \
            title='Distance Histograms')
        
        #Handle saving web resources locally.
        #javascript file
        javascript_path = path.join(opts.dir_path,'js')
        try:
            mkdir(javascript_path)
        except OSError:     #raised if dir exists
            pass
        js_out = open(javascript_path+'/histograms.js','w')
        js_out.write(open(qiime_dir+'js/histograms.js').read())
        js_out.close()
        
    monte_carlo_iters = opts.monte_carlo_iters
    if monte_carlo_iters > 0:
        #Do Monte Carlo for all fields
        monte_carlo_group_distances(mapping_file=opts.map_fname,\
            dmatrix_file=opts.distance_matrix_file,\
            prefs=prefs, \
            dir_prefix = opts.dir_path,\
            fields=fields,\
            default_iters=monte_carlo_iters)
            
        #Do Monte Carlo for within and between fields
        monte_carlo_group_distances_within_between(\
            single_field=within_distances,\
            paired_field=between_distances, dmat=dmat, \
            dir_prefix = opts.dir_path,\
            num_iters=monte_carlo_iters)


if __name__ == "__main__":
    main()
