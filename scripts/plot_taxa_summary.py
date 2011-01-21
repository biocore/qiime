#!/usr/bin/env python
# File created on 19 Jan 2011
from __future__ import division

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Jesse Stombaugh","Julia Goodrich", "Justin Kuczynski"]
__license__ = "GPL"
__version__ = "1.2.0-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Development"
"""
This script generates taxonomy charts
"""

from qiime.util import parse_command_line_parameters, get_qiime_project_dir
from optparse import make_option
from qiime.plot_taxa_summary import make_all_charts
from cogent.util.misc import get_random_directory_name
from qiime.colors import taxonomy_color_prefs_and_map_data_from_options
import re
import matplotlib
import os
import shutil

script_info={}
script_info['brief_description']="""Make taxaonomy summary charts based on taxonomy assignment"""
script_info['script_description']="""This script automates the construction of pie and area charts showing the breakdown of taxonomy by given levels. The script creates an html file for easy visualization of all of the charts on the same page. It uses the taxonomy or category counts from summarize_taxa.py for combined samples by level (-i) and user specified labels for each file passed in (-l). Output will be in a randomly generated folder name within the user specified folder (-o) the default is the current working directory. There is also additional functionality that breaks each taxonomic level up by sample (-s). This will create a chart for each sample at each specified level. The user can also specify the number of categories displayed in a single pie charts the rest are grouped together as 'other category' (-n) default is 20.
"""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Examples:""","""If you wish to run the code using default parameters, you must supply a counts file (Class.txt) along with the taxon level label (Class), by using the following command:""","""plot_taxa_summary.py -i Class.txt -l Class"""))
script_info['script_usage'].append(("""""","""If you want to make pie charts for multiple levels at a time (phylum.txt,class.txt,genus.txt) use the following command:""","""plot_taxa_summary.py -i phylum.txt,class.txt,genus.txt -l phylum,class,genus"""))
script_info['script_usage'].append(("""""","""If you want specify an output directory (e.g. "pie_charts/", regardless of whether the directory exists, use the following command:""","""plot_taxa_summary.py -i Class.txt -l Class -o pie_charts/"""))
script_info['script_usage'].append(("""""","""Additionally, if you would like to display on a set number of taxa ("-n 10") and generate pie charts for all samples ("-s"), you can use the following command:""","""plot_taxa_summary.py -i Class.txt -l Class -o pie_charts/ -n 10 -s"""))
script_info['script_usage'].append(("""""","""If you would like to display generate pie charts for samples samples: 'sample1' and 'sample2' that are in the counts file header, you can use the following command:""","""plot_taxa_summary.py -i Class.txt -l Class -o pie_charts/ -b sample1,sample2"""))
script_info['output_description']="""The script generates an output folder, which contains several files. For each pie chart there is a png and a pdf file. The best way to view all of the pie charts is by opening up the file taxonomy_summary_pie_chart.html."""

script_info['required_options']=[\
make_option('-i', '--input_files', dest='counts_fname',\
    action='store',type='string',\
    help='list of files with sample counts by taxonomy [REQUIRED]'),
make_option('-l', '--labels', dest='labels',action='store',type='string',
    help='list of labels for pie chart (i.e. Phylum, Class) [REQUIRED]'),
make_option('-c', '--chart_type', dest='chart_type',\
    action='store',type='string',\
    help='type of chart to plot (i.e. pie or area). The user has the ability\
 to plot multiple types, by using a comma-separated list (e.g. area,pie)\
 [REQUIRED]')
]

script_info['optional_options']=[\
#make_option('-s', '--sample_flag', dest='do_sample',
#    help='if -s is passed, pie charts will be created for each sample',
#        default=False, action = 'store_true'),
make_option('-n', '--num', dest='num_categories', \
    help='Maximum number of individual categories in each pie chart. \
All additional categories are grouped into an "other" category. \
NOTE: this is only used for the pie charts. \
[default: %default]', default='20'),
make_option('-o', '--dir-prefix', dest='dir_path',\
    help='output folder'),
make_option('-b', '--colorby', dest='colorby',\
    help='This is the samples to make pie charts for in the counts files from\
summarize_taxa.py. The sample name must match the name of a sample id \
in the header of the counts file exactly and multiple categories can be \
list by comma separating them without spaces. If you want to see the pie charts\
 broken up by all samples -s is still funtional. If -s is set and -b is used \
 it will just be broken up by all samples. If neither -s or -b are set the \
 pie charts will be based on all samples put together, one for each level. \
 [default: %default]',default=None),
 make_option('-p', '--prefs_path',help='This is the user-generated preferences \
file. NOTE: This is a file with a dictionary containing preferences for the \
analysis. The label taxonomy_coloring is used for the coloring, see example \
prefs file preferences_file. [default: %default]'),
 make_option('-k', '--background_color',help='This is the background color to \
use in the plots. [default: %default]')
]

script_info['version']=__version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    #Check the version of Matplotlib
    matplotlib_version = re.split("[^\d]", matplotlib.__version__)
    matplotlib_version_info = tuple([int(i) for i in matplotlib_version if \
                            i.isdigit()])
    if matplotlib_version_info != (0,98,5,3) and \
        matplotlib_version_info != (0,98,5,2):
        print "This code was only tested with Matplotlib-0.98.5.2 and \
              Matplotlib-0.98.5.3"

    #get QIIME directory
    qiime_dir=get_qiime_project_dir()

    if not opts.counts_fname:
        option_parser.error("A list of input files must be specified")
    if not opts.labels:
        option_parser.error(\
            "A list of label names cooresponding to files must be specified")

    #get color preferences
    color_prefs, color_data, background_color, label_color= \
                   taxonomy_color_prefs_and_map_data_from_options(opts)
    
    # set the colorby options
    #colorby = opts.colorby
    #do_sample = opts.do_sample
    '''
    if colorby is None and not do_sample:
        colorby = None
    elif colorby is not None and not do_sample:
        colorby = colorby.strip().strip("'").split(',')
    else:
        colorby = []
        for c in color_data['counts'].values():
            colorby.extend(c[0])
        colorby = set(colorby)
        
    '''
    colorby = opts.colorby
    if colorby==None:
        colorby=[]
        for c in color_data['counts'].values():
            colorby.extend(c[0])
    else:
        colorby=colorby.strip().strip("'").split(',')

    #colorby = list(set(colorby)
    
    counts_fname = opts.counts_fname
    
    #Define labels to use
    labels = opts.labels
    data = [(label,f.strip()) \
            for f,label in zip(counts_fname.split(","),labels.split(","))]
    filepath=data[0][1]
    
    filename=filepath.strip().rpartition('/')[0]
    num_categories = int(opts.num_categories)

    #create directory path
    if opts.dir_path:
        if os.path.exists(opts.dir_path):
            dir_path=opts.dir_path
        else:
            try:
                os.mkdir(opts.dir_path)
                dir_path=opts.dir_path
            except OSError:
                pass
    else:
        dir_path='./'
        
    if dir_path == './':
        dir_path = os.getcwd()
    
    #data_dir_path = get_random_directory_name(output_dir=dir_path)
    
    #make javascript output directory
    javascript_path = os.path.join(dir_path,'js')
    try:
        os.mkdir(javascript_path)
    except OSError: #raised if dir exists
        pass
        
    # move javascript file to javascript output directory
    shutil.copyfile(os.path.join(qiime_dir,'qiime','support_files',\
                    'js/overlib.js'),\
                    os.path.join(javascript_path,'overlib.js'))

    #make css output directory
    css_path = os.path.join(dir_path,'css')
    try:
        os.mkdir(css_path)
    except OSError: #raised if dir exists
        pass
        
    # move css file to css output directory
    shutil.copyfile(os.path.join(qiime_dir,'qiime','support_files',\
                    'css/qiime_style.css'),\
                    os.path.join(css_path,'qiime_style.css'))

    plots_to_make=opts.chart_type.split(',')
    for i in plots_to_make:
        chart_type=i.lower().strip()
        #make pie chart output path
        charts_path = os.path.join(dir_path,'charts')
        try:
            os.mkdir(charts_path)
        except OSError: #raised if dir exists
            pass
        
        make_all_charts(data,dir_path,filename,num_categories, \
        colorby,args,color_data, color_prefs,background_color,label_color,\
        chart_type)
        
    
if __name__ == "__main__":
    main()
