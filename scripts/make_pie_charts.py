#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Julia Goodrich"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Julia Goodrich"]
__license__ = "GPL"
__version__ = "0.92-dev"
__maintainer__ = "Julia Goodrich"
__email__ = "julia.goodrich@colorado.edu"
__status__ = "Pre-release"
"""
This script generates taxonomy pie charts
"""

from qiime.util import parse_command_line_parameters, get_qiime_project_dir
from optparse import make_option
from qiime.make_pie_charts import make_all_pie_charts
from qiime.pycogent_backports.misc import get_random_directory_name
from qiime.colors import sample_color_prefs_and_map_data_from_options
import re
import matplotlib
import os
import shutil

script_info={}
script_info['brief_description']="""Make pie charts based on taxonomy assignment"""
script_info['script_description']="""This script automates the construction of pie charts showing the breakdown of taxonomy by given levels. The script creates an html file for easy visualization of all of the pie charts on the same page. It uses the taxonomy or category counts from summarize_taxa.py for combined samples by level (-i) and user specified labels for each file passed in (-l). Output will be in a randomly generated folder name within the user specified folder (-o) the default is the current working directory. There is also additional functionality that breaks each taxonomic level up by sample (-s). This will create a pie chart for each sample at each specified level. The user can also specify the number of categories displayed in a single pie charts the rest are grouped together as 'other category' (-n) default is 20.
"""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Examples:""","""If you wish to run the code using default parameters, you must supply a counts file (Class.txt) along with the taxon level label (Class), by using the following command:""","""make_pie_charts.py -i Class.txt -l Class"""))
script_info['script_usage'].append(("""""","""If you want to make pie charts for multiple levels at a time (phylum.txt,class.txt,genus.txt) use the following command:""","""make_pie_charts.py -i phylum.txt,class.txt,genus.txt -l phylum,class,genus"""))
script_info['script_usage'].append(("""""","""If you want specify an output directory (e.g. "pie_charts/", regardless of whether the directory exists, use the following command:""","""make_pie_charts.py -i Class.txt -l Class -o pie_charts/"""))
script_info['script_usage'].append(("""""","""Additionally, if you would like to display on a set number of taxa ("-n 10") and generate pie charts for all samples ("-s"), you can use the following command:""","""make_pie_charts.py -i Class.txt -l Class -o pie_charts/ -n 10 -s"""))
script_info['output_description']="""The script generates an output folder, which contains several files. For each pie chart there is a png and a pdf file. The best way to view all of the pie charts is by opening up the file taxonomy_summary_pie_chart.html."""
script_info['required_options']=[\
make_option('-i', '--input_files', dest='counts_fname',\
    action='store',type='string',
    help='list of files with sample counts by taxonomy [REQUIRED]'),
make_option('-l', '--labels', dest='labels',action='store',type='string',
            help='list of labels for pie chart(i.e. Phylum,Class)[REQUIRED]')
]

script_info['optional_options']=[\
make_option('-s', '--sample_flag', dest='do_sample',
    help='if True pie charts will be created for each sample',default=False,
                      action = 'store_true'),
make_option('-n', '--num', dest='num_categories', \
             help='Maximum number of individual categories in each pie chart. \
All additional categories are grouped into an "other" category. \
[default: %default]', default='20'),
make_option('-o', '--dir-prefix', dest='dir_path',\
             help='directory prefix for all analyses')
]

script_info['version']=__version__



def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    matplotlib_version = re.split("[^\d]", matplotlib.__version__)
    matplotlib_version_info = tuple([int(i) for i in matplotlib_version if \
                            i.isdigit()])
    if matplotlib_version_info != (0,98,5,3) and \
        matplotlib_version_info != (0,98,5,2):
        print "This code was only tested with Matplotlib-0.98.5.2 and \
              Matplotlib-0.98.5.3"

    if not opts.counts_fname:
        option_parser.error("A list of input files must be specified")
    if not opts.labels:
        option_parser.error("A list of label names cooresponding to files must\
 be specified")

    do_sample = opts.do_sample
    counts_fname = opts.counts_fname
    labels = opts.labels
    data = [(label,f.strip()) \
            for f,label in zip(counts_fname.split(","),labels.split(","))]
    filepath=data[0][1]
    filename=filepath.strip().rpartition('/')[0]
    num_categories = int(opts.num_categories)

    qiime_dir=get_qiime_project_dir()

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


    data_dir_path = get_random_directory_name(output_dir=dir_path)

    charts_path = os.path.join(data_dir_path,'pie_charts')
    try:
        os.mkdir(charts_path)
    except OSError:     #raised if dir exists
        pass


    javascript_path = \
            os.path.join(data_dir_path,'js')
    try:
        os.mkdir(javascript_path)
    except OSError:     #raised if dir exists
        pass
    shutil.copyfile(os.path.join(qiime_dir,'qiime','support_files','js/overlib.js'),\
                                     os.path.join(javascript_path,'overlib.js'))
    css_path = \
            os.path.join(data_dir_path,'css')
    try:
        os.mkdir(css_path)
    except OSError:     #raised if dir exists
        pass
    shutil.copyfile(os.path.join(qiime_dir,'qiime','support_files','css/qiime_style.css'),\
                                os.path.join(css_path,'qiime_style.css'))


    make_all_pie_charts(data,data_dir_path,filename,num_categories, do_sample,args)

if __name__ == "__main__":
    main()
