#!/usr/bin/env python
# File created on 09 Feb 2010
#file make_3d_plots.py

from __future__ import division

#!/usr/bin/env python
#file compare_3d_plots.py

__author__ = "Dan Knights"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Dan Knights"] #remember to add yourself
__license__ = "GPL"
__version__ = "0.92-dev"
__maintainer__ = "Dan Knights"
__email__ = "daniel.knights@colorado.edu"
__status__ = "Pre-release"

from qiime.util import parse_command_line_parameters
from optparse import make_option
from qiime.make_3d_plots import generate_3d_plots,\
get_coord,get_map,remove_unmapped_samples,\
process_coord_filenames,get_multiple_coords,\
process_colorby,create_dir
from qiime.parse import parse_coords,group_by_field,group_by_fields
import shutil
import os
from random import choice
from time import strftime
from qiime.util import get_qiime_project_dir
                                

#from qiime.make_3d_plots import get_coord, get_map, remove_unmapped_samples,\
#process_colorby, process_custom_axes, process_coord_filenames, \
#get_multiple_coords, get_custom_coords, remove_nans,\
#scale_custom_coords, create_dir, _do_3d_plots
#import shutil
#from optparse import OptionParser
#from random import choice
#from time import strftime
#import sys
#import os.path
#import numpy as np

script_info={}
script_info['brief_description']="""Plot two PCoA files on the same 3D plot"""
script_info['script_description']="""This script generates a 3D plot comparing two sets of principle coordinates using as input two principle coordinates files. The principle coordinates files are obtained by applying "principle_coordinates.py" to a file containing beta diversity measures. The beta diversity files are optained by applying "beta_diversity.py" to an OTU table. One may apply "transform_coordinate_matrices.py" to the principle coordinates files before using this script to compare them."""
script_info['script_usage']=[]
script_info['script_usage'].append(("Example 1","""Compare two pca/pcoa files in the same 3d plot where each sample ID is assigned its own color:""","""compare_3d_plots.py -i 'raw_pca_data1.txt,raw_pca_data2.txt'"""))
script_info['script_usage'].append(("Example 2","""Compare two pca/pcoa files in the same 3d plot with two coloring schemes (Day and Type):""","""compare_3d_plots.py -i 'raw_pca_data1.txt,raw_pca_data2.txt' -m input_map.txt -b 'Day,Type'"""))
script_info['script_usage'].append(("Example 3","""Compare two pca/pcoa files in the same 3d plot for a combination of label headers from a mapping file: ""","""compare_3d_plots.py -i 'raw_pca_data1.txt,raw_pca_data2.txt' -m input_map.txt -b 'Type&&Day' -o ./test/"""))
script_info['output_description']="""This script results in a folder containing an html file which displays the 3D Plots generated."""
script_info['required_options']= [\
    make_option('-i', '--coord_fnames', \
        help='This is comma-separated list of the paths to the principal \
coordinates files (i.e., resulting file \
from principal_coordinates.py), e.g \'pcoa1.txt,pcoa2.txt\'')
]

script_info['optional_options']= [\
 make_option('-m', '--map_fname', dest='map_fname', \
     help='This is the user-generated mapping file [default=%default]'),
 make_option('-b', '--colorby', dest='colorby',\
     help='This is a list of the categories to color by in the plots from the \
user-generated mapping file. The categories must match the name of a column \
header in the mapping file exactly and multiple categories can be list by comma \
separating them without spaces. The user can also combine columns in the \
mapping file by separating the categories by "&&" without spaces \
[default=%default]'),
 make_option('-a', '--custom_axes',help='This is a category or list of \
categories from the user-generated mapping file to use as a custom axis in the \
plot.  For instance, if there is a pH category and one would like to see \
the samples plotted on that axis instead of PC1, PC2, etc., one can use \
this option.  It is also useful for plotting time-series data \
[default: %default]'),
 make_option('-p', '--prefs_path',help='This is the user-generated preferences \
file. NOTE: This is a file with a dictionary containing preferences for the \
analysis. See make_3d_plot_prefs_file.py. [default: %default]'),
 make_option('-o', '--dir_path', dest='dir_path',help='This is the location \
where the resulting output should be written [default=%default]',default='')
]
script_info['version'] = __version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)

    if len(opts.coord_fnames.split(',')) < 2:
        parser.error('Please provide at least two coordinate files')

    data = {}

    #Open and get coord data (for multiple coords files)
    coord_files = process_coord_filenames(opts.coord_fnames)
    num_coord_files = len(coord_files)
    data['edges'], data['coord'] = get_multiple_coords(coord_files)

    #Open and get mapping data, if none supplied create a pseudo mapping \
    #file
    if opts.map_fname:
        data['map'] = get_map(opts, data)
    else:
        data['map']=(([['#SampleID','Sample']]))
        for i in range(len(data['coord'][0])):
            data['map'].append([data['coord'][0][i],'Sample'])

    # duplicate samples in mapping file for all coord files being compared
    newmap = [data['map'][0]]
    for i in xrange(len(coord_files)):
        for sample in data['map'][1:]:
            newsample = ['%s_%d' %(sample[0],i)]
            newsample.extend(sample[1:])
            newmap.append(newsample)
    data['map'] = newmap


    # remove any samples not present in mapping file
    remove_unmapped_samples(data['map'],data['coord'])


    #Determine which mapping headers to color by, if none given, color by all 
    # columns in map file
    if opts.prefs_path:
        prefs = eval(open(opts.prefs_path, 'U').read())
        prefs, data=process_colorby(None, data, prefs)
    elif opts.colorby:
        prefs,data=process_colorby(opts.colorby,data)
    else:
        default_colorby = ','.join(data['map'][0])
        prefs,data=process_colorby(default_colorby,data)
        prefs={'Sample':{'column':'#SampleID'}}

    # process custom axes, if present.
    custom_axes = None
    if opts.custom_axes:
        custom_axes = process_custom_axes(options.custom_axes)
        get_custom_coords(custom_axes, data['map'], data['coord'])
        remove_nans(data['coord'])
        scale_custom_coords(custom_axes,data['coord'])

    # Generate random output file name and create directories
    dir_path = opts.dir_path
    if dir_path and not dir_path.endswith('/'):
        dir_path = dir_path + '/'
    dir_path = create_dir(dir_path,'3d_plots_') 

    alphabet = "ABCDEFGHIJKLMNOPQRSTUZWXYZ"
    alphabet += alphabet.lower()
    alphabet += "01234567890"

    qiime_dir=get_qiime_project_dir()

    jar_path=os.path.join(qiime_dir,'qiime/support_files/jar/')

    data_file_path=''.join([choice(alphabet) for i in range(10)])
    data_file_path=strftime("%Y_%m_%d_%H_%M_%S")+data_file_path
    data_file_dir_path = dir_path+data_file_path

    data_file_dir_path=create_dir(data_file_dir_path,'')
    jar_dir_path = create_dir(os.path.join(dir_path,'jar/'),'')
    shutil.copyfile(os.path.join(jar_path,'king.jar'), jar_dir_path+'king.jar')

    filepath=opts.coord_fnames[0]
    filename=filepath.strip().split('/')[-1]
    
    try:
        action = generate_3d_plots
    except NameError:
        action = None
    #Place this outside try/except so we don't mask NameError in action
    if action:
        action(prefs, data, custom_axes, dir_path, data_file_path,filename)

if __name__ == "__main__":
    main()
