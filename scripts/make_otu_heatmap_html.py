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


from qiime.util import parse_command_line_parameters, get_options_lookup
from optparse import make_option
from qiime.make_otu_heatmap_html import generate_heatmap_plots,get_otu_counts
from qiime.make_3d_plots import create_dir
import shutil
import os
from qiime.util import get_qiime_project_dir

options_lookup = get_options_lookup()

#make_otu_heatmap_html.py
script_info={}
script_info['brief_description']="""Make heatmap of OTU table"""
script_info['script_description']="""Once the OTU table has been generated, the user can create an interactive OTU heatmap. This script parses the OTU count table and filters the table by counts per otu (user-specified), then converts the table into a javascript array, which can be loaded into a web application. The OTU heatmap displays raw OTU counts per sample, where the counts are colored based on the contribution of each OTU to the total OTU count present in that sample (blue: contributes low percentage of OTUs to sample; red: contributes high percentage of OTUs). This web application allows the user to filter the otu table by number of counts per otu. The user also has the ability to view the table based on taxonomy assignment. Additional features include: the ability to drag rows (up and down) by clicking and dragging on the row headers; and the ability to zoom in on parts of the heatmap by clicking on the counts within the heatmap."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Examples:""","""By using the default values ("-n 5), you can then use the code as follows:""","""%prog -i otu_table.txt"""))
script_info['script_usage'].append(("""""","""If you would like to filter the OTU table by a different number of counts per OTU (i.e., 10), you can use the following code:""","""%prog -i otu_table.txt -n 10"""))
script_info['script_usage'].append(("""""","""If you would like to specify a different output directory (i.e., "otu_heatmap"), you can use the following code:""","""%prog -i otu_table.txt -o otu_heatmap"""))
script_info['output_description']="""The interactive heatmap is located in a randomly generated folder where the name of the folder starts with "otu_table". The resulting folder contains the interactive heatmap (html file) along with a javascript library folder. This web application has been tested in Mozilla Firefox and Safari. Safari is recommended for viewing the OTU Heatmap, since the HTML table generation is much faster."""
script_info['required_options']=[\
 options_lookup['otu_table_as_primary_input']
]
script_info['optional_options']=[\
options_lookup['output_dir'],
 make_option('-n', '--num_otu_hits', help='This is the minimum number of \
Samples that an OTU is present in, for an OTU to be kept in the OTU table \
[default: %default]',default=5)
]

script_info['version'] = __version__





def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
      
    data = {}

    #Open and get coord data
    data['otu_counts'] = get_otu_counts(opts.otu_table_fp, data)

    filepath=opts.otu_table_fp
    filename=filepath.strip().split('/')[-1].split('.')[0]

    dir_path = create_dir(opts.output_dir,'otu_heatmap_')

    if dir_path and not dir_path.endswith('/'):
        dir_path=dir_path+'/'

    js_dir_path = create_dir(os.path.join(dir_path,'js/'),'')

    qiime_dir=get_qiime_project_dir()

    js_path=os.path.join(qiime_dir,'qiime/support_files/js')

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