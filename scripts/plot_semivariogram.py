#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Antonio Gonzalez Pena"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Antonio Gonzalez Pena"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Antonio Gonzalez Pena"
__email__ = "antgonza@gmail.com"
__status__ = "Release"
 
from qiime.util import parse_command_line_parameters, get_options_lookup
from optparse import make_option
from qiime.plot_semivariogram import fit_semivariogram, FitModel
from qiime.parse import parse_distmat
from pylab import plot, xlabel, ylabel, title, savefig
from numpy import asarray
import os

options_lookup = get_options_lookup()

script_info={}

script_info['brief_description']="""Fits a model between two distance matrices and plots the result"""
script_info['script_description']="""Fits a spatial autocorrelation model between two matrices and plots the result. This script will work with two distance matrices but will ignore the 0s at the diagonal and the values that go to N/A"""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Fitting""","""For this script, the user supplies two distance matrices (i.e. resulting file from beta_diversity.py), along with the output filename (e.g. semivariogram), and the model to fit, as follows:""","""plot_semivariogram.py -x distance.txt -y unifrac.txt -m exponential -o semivariogram.png"""))
script_info['output_description']="""The resulting output file consists of a pdf image containing the plot between the two distances matrices and the fitted model"""

script_info['required_options']=[\
 make_option('-x', '--input_path_x',\
     help='path to distance matrix to be displayed in the x axis'),\
 make_option('-y', '--input_path_y',\
     help='path to distance matrix to be displayed in the y axis'),\
 make_option('-m', '--model', type='choice',\
     choices=FitModel.options, default='exponential', 
     help='model to be fitted to the data. Valid ' +\
     'choices are:' + ', '.join(FitModel.options) + '. [default: %default]'),\
 make_option('-o', '--output_path',
     help='output path. directory for batch processing, '+\
       'filename for single file operation'),\
]
script_info['optional_options']=[\
 make_option('-b', '--binning', type='string',\
     default=None, help='binning ranges. Format: [increment,top_limit], when ' +\
     'top_limit is -1=infinitum; you can specify several ranges using the same ' +\
     'format, i.e. [2.5,10][50,-1] will set two bins, one from 0-10 using 2.5 ' +\
     'size steps and from 10-inf using 50 size steps. Note that the binning is ' +\
     'used to clean the plots (reduce number of points) but ignored to fit the ' +\
     'model. [default: %default]'),\
]

script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    if opts.binning is None:
        ranges = []
    else:
        # simple ranges format validation
        if opts.binning.count('[')!=opts.binning.count(']') or\
          opts.binning.count('[')!=opts.binning.count(','):
            raise ValueError, "The binning input has an error: '%s'; " % opts.binning +\
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
            values = map(float,r.split(','))
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
    (x_val,y_val,x_fit,y_fit) = fit_semivariogram(x_distmtx, y_distmtx, opts.model, ranges)
    
    plot(x_val, y_val, 'o', color="white")   
    plot(x_fit, y_fit, linewidth=2.0, color="blue")
    
    x_label = 'Distance (m)'
    y_label = 'Community Dissimilarity'
    fig_title = 'Semivariogram (%s)' % opts.model
    
    xlabel(x_label)
    ylabel(y_label)
    title(fig_title)
    
    savefig(opts.output_path)

if __name__ == "__main__":
    main()
