#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Antonio Gonzalez Pena, Kyle Patnode"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Development"

from qiime.colors import data_colors, data_color_order
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.util import make_option
from qiime.plot_semivariogram import fit_semivariogram, FitModel
from qiime.parse import parse_distmat, parse_mapping_file
from qiime.filter import filter_samples_from_distance_matrix, sample_ids_from_metadata_description
from pylab import plot, xlabel, ylabel, title, savefig, ylim, xlim, legend, show, figure
from numpy import asarray
import os
from StringIO import StringIO

options_lookup = get_options_lookup()

script_info={}

script_info['brief_description']="""Fits a model between two distance matrices and plots the result"""
script_info['script_description']="""Fits a spatial autocorrelation model between two matrices and plots the result. This script will work with two distance matrices but will ignore the 0s at the diagonal and the values that go to N/A"""
script_info['script_usage']=[]
script_info['script_usage'].append(
      ("Fitting",
       "For this script, the user supplies two distance matrices (i.e. resulting file from beta_diversity.py), along with the output filename (e.g. semivariogram), and the model to fit, as follows:",
       "%prog -x distance.txt -y unifrac.txt -o semivariogram_exponential.png"))
script_info['script_usage'].append(
      ("",
       "Modify the the default method to gaussian",
       "%prog -x distance.txt -y unifrac.txt -m gaussian -o semivariogram_gaussian.png"))
script_info['output_description']="""The resulting output file consists of a pdf image containing the plot between the two distances matrices and the fitted model"""

script_info['required_options']=[\
 make_option('-x', '--input_path_x',type='existing_filepath',\
     help='path to distance matrix to be displayed in the x axis'),\
 make_option('-y', '--input_path_y',type='existing_filepath',\
     help='path to distance matrix to be displayed in the y axis'),\
 make_option('-o', '--output_path',type='new_path',
     help='output path. directory for batch processing, '+\
       'filename for single file operation'),

]
script_info['optional_options']=[\
 make_option('-b', '--binning', type='string',\
     default=None, help='binning ranges. Format: [increment,top_limit], when ' +\
     'top_limit is -1=infinitum; you can specify several ranges using the same ' +\
     'format, i.e. [2.5,10][50,-1] will set two bins, one from 0-10 using 2.5 ' +\
     'size steps and from 10-inf using 50 size steps. Note that the binning is ' +\
     'used to clean the plots (reduce number of points) but ignored to fit the ' +\
     'model. [default: %default]'),
 make_option('--ignore_missing_samples', help='This will overpass the error raised ' +\
     'when the matrices have different sizes/samples', action='store_true', default=False),         
 make_option('--x_max', type='float', help='x axis max limit [default: auto]', default=None),         
 make_option('--x_min', type='float', help='x axis min limit [default: auto]', default=None),         
 make_option('--y_max', type='float', help='y axis max limit [default: auto]', default=None),         
 make_option('--y_min', type='float', help='y axis min limit [default: auto]', default=None),
 make_option('-X', '--x_label', default='Distance Dissimilarity (m)',type='string',\
     help='Label for the x axis [default: %default]'),
 make_option('-Y', '--y_label', default='Community Dissimilarity',type='string',\
     help='Label for the y axis [default: %default]'),
 make_option('-t', '--fig_title', default='Semivariogram',type='string',\
     help='Title of the plot [default: %default]'),       
 make_option('--dot_color', type='string', help='dot color for plot, more info:' +\
    ' http://matplotlib.sourceforge.net/api/pyplot_api.html' +\
    ' [default: %default]', default="white"), 
 make_option('--dot_marker', type='string', help='dot color for plot, more info:' +\
    ' http://matplotlib.sourceforge.net/api/pyplot_api.html' +\
    ' [default: %default]', default="o"), 
 make_option('--line_color', type='string', help='line color for plot, more info:' +\
    ' http://matplotlib.sourceforge.net/api/pyplot_api.html' +\
    ' [default: %default]', default="blue"), 
 make_option('--dot_alpha', type='float', help='alpha for dots, more info:' +\
    ' http://matplotlib.sourceforge.net/api/pyplot_api.html' +\
    ' [default: %default]', default=1),
 make_option('--line_alpha', type='float', help='alpha for dots, more info:' +\
    ' http://matplotlib.sourceforge.net/api/pyplot_api.html' +\
    ' [default: %default]', default=1),
 make_option('--model', type='choice',\
     choices=FitModel.options, default='exponential', 
     help='model to be fitted to the data. Valid ' +\
     'choices are:' + ', '.join(FitModel.options) + '. [default: %default]'),\
 make_option('-p', '--print_model', action='store_true',
     help='Print in the title of the plot the function of the fit. ' +\
     '[default: %default]',default=False),
 make_option('-c', '--category', type='string', help='category to color the '
    'plots by [default: %default]', default=None),
 make_option('-m', '--mapping_fp', type='existing_filepath', help='metadata '
    'mapping file, only used when coloring by a category [default: %default]',
    default=None)
]

script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    category = opts.category
    mapping_fp = opts.mapping_fp

    if (category and mapping_fp == None) or (category == None and mapping_fp):
        option_parser.error('If coloring by a metadata category, both the '
            'category and the mapping file must be supplied.')
    elif mapping_fp and category:
        mapping_data, mapping_headers, _ = parse_mapping_file(open(mapping_fp,
            'U'))
        if category not in mapping_headers:
            option_parser.error("The category supplied must exist in the "
                "metadata mapping file, '%s' does not exist." % category)
        index = mapping_headers.index(category)
        categories = list(set([line[index] for line in mapping_data]))
    list_of_plots = []

    if opts.binning is None:
        ranges = []
    else:
        # simple ranges format validation
        if opts.binning.count('[')!=opts.binning.count(']') or\
          opts.binning.count('[')!=opts.binning.count(','):
            raise ValueError, "The binning input has an error: '%s'; " % +\
             "\nthe format should be [increment1,top_limit1][increment2,top_limit2]" 
        # spliting in ranges
        rgn_txt = opts.binning.split('][')
        # removing left [ and right ]
        rgn_txt[0] = rgn_txt[0][1:]
        rgn_txt[-1] = rgn_txt[-1][:-1]
        # converting into int
        ranges = []
        max = 0
        
        for i,r in enumerate(rgn_txt):
            try:
                values = map(float,r.split(','))
            except ValueError:
                raise ValueError, "Not a valid format for binning %s" % opts.binning 
            if len(values)!=2:
                raise ValueError, "All ranges must have only 2 values: [%s]" % r
            elif i+1!=len(rgn_txt): 
                if values[0]>values[1]:
                    raise ValueError, "The bin value can't be greater than the max value: [%s]" % r
                elif values<0:
                    raise ValueError, "This value can not be negative: [%s]" % r
                elif max>values[1]:
                    raise ValueError, "This value can not smaller than the previous one: [%s]" % r
                else:
                    max=values[1]
            
            ranges.append(values)
    
    x_samples, x_distmtx = parse_distmat(open(opts.input_path_x,'U'))
    y_samples, y_distmtx = parse_distmat(open(opts.input_path_y,'U'))
    
    if opts.ignore_missing_samples:
        ignoring_from_x = list(set(x_samples)-set(y_samples))
        ignoring_from_y = list(set(y_samples)-set(x_samples))
        
        if opts.verbose:
            print '\nFrom %s we are ignoring: %s\n' % (opts.input_path_x, ignoring_from_x)
            print '\nFrom %s we are ignoring: %s\n' % (opts.input_path_y, ignoring_from_y)
            print '\nOnly using: %s\n' % (list(set(x_samples) & set(y_samples)))
        
        x_file = StringIO(\
            filter_samples_from_distance_matrix((x_samples, x_distmtx), ignoring_from_x))
        x_samples, x_distmtx = parse_distmat(x_file)
        
        y_file = StringIO(\
            filter_samples_from_distance_matrix((y_samples, y_distmtx), ignoring_from_y))
        y_samples, y_distmtx = parse_distmat(y_file)
    else:
        if x_distmtx.shape!=y_distmtx.shape:
            raise ValueError, 'The distance matrices have different sizes. ' +\
                'You can cancel this error by passing --ignore_missing_samples'

    figure(figsize=(3,3))
    if category == None:
        x_val, y_val, x_fit, y_fit, func_text = fit_semivariogram(
            (x_samples,x_distmtx), (y_samples,y_distmtx), opts.model, ranges)
        
        plot(x_val, y_val, color=opts.dot_color, marker=opts.dot_marker, linestyle="None", alpha=opts.dot_alpha)
        plot(x_fit, y_fit, linewidth=2.0, color=opts.line_color, alpha=opts.line_alpha)
    else:
        for single_category, color_key in zip(categories, data_color_order):
            good_sample_ids = sample_ids_from_metadata_description(
                open(mapping_fp), '%s:%s' % (category, single_category))

            _y_samples, _y_distmtx = parse_distmat(StringIO(
                filter_samples_from_distance_matrix((y_samples, y_distmtx),
                good_sample_ids, negate=True)))
            _x_samples, _x_distmtx = parse_distmat(StringIO(
                filter_samples_from_distance_matrix((x_samples, x_distmtx),
                good_sample_ids, negate=True)))

            x_val, y_val, x_fit, y_fit, func_text = fit_semivariogram(
                (_x_samples, _x_distmtx), (_y_samples, _y_distmtx),
                opts.model,ranges)
            color_only = str(data_colors[color_key]).split(':')[1]
            plot(x_val, y_val, color=color_only,
                marker=opts.dot_marker, linestyle="None", alpha=opts.dot_alpha)
            plot(x_fit, y_fit, linewidth=2.0, color=color_only,
            alpha=opts.line_alpha, label=single_category)

    legend(loc=0, bbox_to_anchor=(0.5, -0.05),
          fancybox=True, shadow=True)


    if opts.x_min!=None and opts.x_max!=None:
        xlim([opts.x_min,opts.x_max])
    if opts.y_min!=None and opts.y_max!=None:
        ylim([opts.y_min,opts.y_max])
        
    x_label = opts.x_label
    y_label = opts.y_label
    fig_title = '%s (%s)' % (opts.fig_title, opts.model)
    
    xlabel(x_label)
    ylabel(y_label)
    if opts.print_model:
        title(fig_title + ' ' + func_text)
    else:
        title(fig_title)
    
    savefig(opts.output_path)

if __name__ == "__main__":
    main()
