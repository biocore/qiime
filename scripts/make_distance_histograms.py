#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Jeremy Widmann"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Jeremy Widmann","Rob Knight"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Jeremy Widmann"
__email__ = "Jeremy.Widmann@colorado.edu"
__status__ = "Pre-release"
 

from qiime.util import parse_command_line_parameters, get_qiime_project_dir
from optparse import make_option
from qiime.make_distance_histograms import group_distances, _make_path, \
    draw_all_histograms, _make_relative_paths, make_main_html, \
    monte_carlo_group_distances
from qiime.make_3d_plots import create_dir
from os import mkdir
from string import strip

script_description = \
"""To visualize the distance between samples and/or categories in the mapping
file, the user can generate histograms to represent the distances between
samples.  This script generates an HTML file, where the user can compare the
distances between samples based on the different categories associated to each
sample in the mapping file."""

script_usage = \
"""
Distance Histograms are a way to compare different categories and see which
tend to have larger/smaller distances than others.  For example, in the hand
study, you may want to compare the distances between hands to the distances
between individuals (with the file "hand_distances.txt" using the parameter -d
hand_distances.txt). The categories are defined in the mapping file  (specified
using the parameter -m hand_map.txt).  If you want to look at the distances
between hands and individuals, choose the "Hand" field and "Individual" field
(using the parameter --fields Hand,Individual (notice the fields are comma
delimited)).  For each of these groups of distances a histogram is made.  The
output is a HTML file ("QIIME_Distance_Histograms.html" when the parameter
--html_output is specified) which is created in the "Distance_Histograms"
directory (using the parameter -o Distance_Histograms to specify output
directory) where you can look at all the distance histograms individually, and
compare them between each other.

$ python qiime_dir/make_distance_histograms.py -d hand_distances.txt -m 
hand_map.txt --fields Treatment,Individual -o Distance_Histograms --html_output
"""

required_options = [\
    make_option('-d','--distance_matrix_file',dest='distance_matrix_file',\
        type='string',help='''Path to distance matrix file.'''),\
    make_option('-m','--mapping_file',dest='mapping_file',type='string',\
        help='''Path to environment mapping file.'''),\
]

optional_options = [\
    make_option('-p','--prefs_file',dest='prefs_file',type='string',\
        help='''File containing prefs for analysis.  NOTE: This is a file with a dict containing preferences for the analysis.  This dict must have a "Fields" key mapping to a list of desired fields.[default: %default]'''),\
    make_option('-o', '--dir_path', dest='dir_path',\
        help='Directory to output data for all analyses. [default: %default]',\
        default='.'),\
    make_option('--monte_carlo',dest='monte_carlo',default=False,\
        action='store_true',help='''Perform Monte Carlo on distances.  [Default: %default]'''),\
    make_option('--html_output',dest='html_output',default=False,\
        action='store_true',help='''Write output in HTML format. [Default: %default]'''),\
    make_option('--fields', dest='fields',\
        help='Comma delimited list of fields to compare.  This overwrites fields in prefs file.  If this is not provided, the first field in mapping file will be used.  Usage: --fields Field1,Field2,Field3'),\
]


def main():
    option_parser, opts, args = parse_command_line_parameters(
        script_description=script_description,
        script_usage=script_usage,
        version=__version__,
        required_options=required_options,
        optional_options=optional_options)
    
    qiime_dir = get_qiime_project_dir()+'/support_files/'
    
    if opts.prefs_file:
        prefs = eval(open(opts.prefs_file, 'U').read())
    else:
        prefs=None
    
    fields = opts.fields
    if fields is not None:
        fields = map(strip,fields.split(','))
    elif opts.prefs_file is not None:
        prefs = eval(open(opts.prefs_file, 'U').read())
        fields = prefs['FIELDS']
    else:
        fields = None
    
    within_distances, between_distances, dmat = \
        group_distances(mapping_file=opts.mapping_file,\
        dmatrix_file=opts.distance_matrix_file,\
        fields=fields,\
        dir_prefix=create_dir(opts.dir_path,'distances'))
    
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
                histogram_dir=histograms_path)
        
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
        logo_out.write(open(qiime_dir+'qiime_header.png').read())
        logo_out.close()
    
    if opts.monte_carlo:
        monte_carlo_group_distances(mapping_file=opts.mapping_file,\
            dmatrix_file=opts.distance_matrix_file,\
            prefs=prefs, \
            dir_prefix = opts.dir_path,\
            fields=fields)


if __name__ == "__main__":
    main()