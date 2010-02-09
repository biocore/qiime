#!/usr/bin/env python
# File created on 09 Feb 2010
#file make_otu_heatmap_html.py

from __future__ import division

__author__ = "Jesse Stombaugh"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Jesse Stombaugh"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Pre-release"


from qiime.util import parse_command_line_parameters
from optparse import make_option
from qiime.make_otu_heatmap_html import generate_heatmap_plots,get_otu_counts
from qiime.make_3d_plots import create_dir
import shutil
import os
from qiime.util import get_qiime_project_dir


script_description = """This script takes an otu table and generates an OTU \
count Heatmap. This Heatmap is a dynamic web application, which allows the user \
to filter the table by counts, transpose the table and several mouseover \
features."""

script_usage = """
Example 1: Create html file and javascript array from otu counts table:

Usage: make_otu_heatmap_html.py -i otu_counts.txt

Example 2: Create html file, then javascript array in specified directory:

Usage: make_otu_heatmap_html.py -i otu_counts.txt -o ./test/

Example 3: Create html file, then javascript array where the number of hits
per otu are specified:

Usage: make_otu_heatmap_html.py -i otu_counts.txt -n 50"""

required_options = [\
 # Example required option
 #make_option('-i','--input_dir',help='the input directory'),\
 make_option('-i', '--otu_count_fname', help='This is the path to the OTU table \
(i.e., the resulting OTU table from make_otu_table.py or filter_otu_table.py)')
]

optional_options = [\
 # Example optional option
 #make_option('-o','--output_dir',help='the output directory [default: %default]'),\
 make_option('-n', '--num_otu_hits', help='This is the minimum number of \
Samples that an OTU is present in, for an OTU to be kept in the OTU table \
[default: %default]',default=5),
 make_option('-o', '--dir_path',help='This is the location where the resulting \
output should be written [default: %default]',default='')
]


def main():
    option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)
      
    data = {}

    #Open and get coord data
    data['otu_counts'] = get_otu_counts(opts, data)

    filepath=opts.otu_count_fname
    filename=filepath.strip().split('/')[-1].split('.')[0]

    dir_path = create_dir(opts.dir_path,'otu_heatmap_')

    if dir_path and not dir_path.endswith('/'):
        dir_path=dir_path+'/'

    js_dir_path = create_dir(os.path.join(dir_path,'js/'),'')

    qiime_dir=get_qiime_project_dir()

    js_path=os.path.join(qiime_dir,'support_files/js')

    shutil.copyfile(os.path.join(js_path,'overlib.js'), js_dir_path+'overlib.js')
    shutil.copyfile(os.path.join(js_path,'otu_count_display.js'), js_dir_path+\
              'otu_count_display.js')
    shutil.copyfile(os.path.join(js_path,'jquery.js'), js_dir_path+'jquery.js')
    shutil.copyfile(os.path.join(js_path,'jquery.tablednd_0_5.js'), js_dir_path+\
              'jquery.tablednd_0_5.js')

    try:
        action = generate_heatmap_plots
    except NameError:
        action = None
    #Place this outside try/except so we don't mask NameError in action
    if action:
        action(opts,data, dir_path,js_dir_path,filename)

if __name__ == "__main__":
    main()