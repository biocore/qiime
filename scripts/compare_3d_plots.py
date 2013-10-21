#!/usr/bin/env python
# File created on 09 Feb 2010
#file compare_3d_plots.py

from __future__ import division

#!/usr/bin/env python
#file compare_3d_plots.py

__author__ = "Dan Knights"
__copyright__ = "Copyright 2011, The QIIME Project"
__credits__ = ["Dan Knights", "Antonio Gonzalez Pena",
                "Jose Antonio Navas Molina"] #remember to add yourself
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Dan Knights"
__email__ = "daniel.knights@colorado.edu"
__status__ = "Development"

from qiime.util import parse_command_line_parameters, get_options_lookup, create_dir
from qiime.util import make_option
from qiime.make_3d_plots import generate_3d_plots,\
get_coord,remove_unmapped_samples,\
get_multiple_coords,\
process_colorby, process_custom_axes, get_custom_coords, remove_nans, scale_custom_coords,\
validate_coord_files
from qiime.parse import parse_coords,group_by_field,group_by_fields
from qiime.colors import get_map
from qiime.colors import sample_color_prefs_and_map_data_from_options
import shutil
import os
from random import choice
from time import strftime
from qiime.util import get_qiime_project_dir
from cogent.util.misc import get_random_directory_name
options_lookup = get_options_lookup()                                

script_info={}
script_info['brief_description']="""Plot several PCoA files on the same 3D plot"""
script_info['script_description']="""This script generates a 3D plot comparing two or more sets of principal coordinates using as input two or more principal coordinates files. Edges are drawn in the plot connecting samples with the same ID across different principal coordinates files. The user can also include a file listing the edges to be drawn in the plot, in which case the user may submit any number of principal coordinates files (including one). If the user includes the edges file, the sample IDs need not match between principal coordinates files.

The principal_coordinates coordinates files are obtained by applying "principal_coordinates.py" to a file containing beta diversity measures. The beta diversity files are optained by applying "beta_diversity.py" to an OTU table. One may apply "transform_coordinate_matrices.py" to the principal_coordinates coordinates files before using this script to compare them."""
script_info['script_usage']=[]
script_info['script_usage'].append(("Example 1","""Compare two pca/pcoa files in the same 3d plot where each sample ID is assigned its own color:""","""%prog -i $PWD/raw_pca_data1.txt,$PWD/raw_pca_data2.txt -m $PWD/Fasting_Map.txt"""))
script_info['script_usage'].append(("Example 2","""Compare two pca/pcoa files in the same 3d plot with two coloring schemes (Treatment and DOB):""","""%prog -i $PWD/raw_pca_data1.txt,$PWD/raw_pca_data2.txt -m $PWD/Fasting_Map.txt -b 'Treatment,DOB'"""))
script_info['script_usage'].append(("Example 3","""Compare two pca/pcoa files in the same 3d plot for a combination of label headers from a mapping file: ""","""%prog -i $PWD/raw_pca_data1.txt,$PWD/raw_pca_data2.txt -m $PWD/Fasting_Map.txt -b 'Treatment&&DOB' -o $PWD/test/"""))
script_info['script_usage'].append(("Example 4","""Pass in a list of desired edges and only one pca/pcoa file: ""","""%prog -i $PWD/raw_pca_data1.txt -e $PWD/edges.txt -m Fasting_Map.txt -b 'Treatment&&DOB' -o $PWD/test2/"""))
script_info['script_usage'].append(("Example 5","""Pass in a list of desired edges and only one pca/pcoa file: ""","""%prog -i $PWD/raw_pca_data1.txt,$PWD/raw_pca_data2.txt -e $PWD/edges.txt -m $PWD/Fasting_Map.txt -b 'Treatment&&DOB' -o $PWD/test3/"""))
script_info['output_description']="""This script results in a folder containing an html file which displays the 3D Plots generated."""
script_info['required_options']= [\
    make_option('-i', '--coord_fnames',type='existing_filepaths',\
        help='This is comma-separated list of the paths to the principal \
coordinates files (i.e., resulting file \
from principal_coordinates.py), e.g \'pcoa1.txt,pcoa2.txt\''),
 make_option('-m', '--map_fname', dest='map_fname',type='existing_filepath', \
     help='This is the user-generated mapping file [default=%default]'),
]

script_info['optional_options']= [\
 make_option('-b', '--colorby', dest='colorby',type='string',\
     help='This is a list of the categories to color by in the plots from the \
user-generated mapping file. The categories must match the name of a column \
header in the mapping file exactly and multiple categories can be list by comma \
separating them without spaces. The user can also combine columns in the \
mapping file by separating the categories by "&&" without spaces \
[default=%default]'),
 make_option('-a', '--custom_axes',type='string',help='This is a category or list of \
categories from the user-generated mapping file to use as a custom axis in the \
plot.  For instance, if there is a pH category and one would like to see \
the samples plotted on that axis instead of PC1, PC2, etc., one can use \
this option.  It is also useful for plotting time-series data \
[default: %default]'),
 make_option('-p', '--prefs_path',help='This is the user-generated preferences \
file. NOTE: This is a file with a dictionary containing preferences for the \
analysis. See make_prefs_file.py. [default: %default]',type='existing_filepath'),
 make_option('-k', '--background_color', type='choice', choices=['black', 'white'],
    help='This is the background color to use in the plots (Options are \
\'black\' or \'white\'. [default: %default]'),
 make_option('-e', '--edges_file',help='A file where each line contains two \
sample IDs separated by a whitespace character; for each pair of sample IDs, \
an edge will be drawn from the first sample to the second sample. \
[default: %default]',default=None,type='existing_filepath'),
 make_option('--serial',action="store_true", \
 help='Connect the 1st set of points to the 2nd, the 2nd to the 3rd, etc. \
Default behavior is to connect each set of points back to the 1st set. This \
flag is ignored if the user specifies an edges file.'),
 options_lookup['output_dir']
]
script_info['version'] = __version__

def main():
    print "\nThis script is being deprecated in favor of make_emperor.py, and will no longer be available in QIIME 1.8.0-dev.\n"

    option_parser, opts, args = parse_command_line_parameters(**script_info)

    prefs, data, background_color, label_color, ball_scale, arrow_colors = \
                            sample_color_prefs_and_map_data_from_options(opts)
    
    if len(opts.coord_fnames) < 2 and opts.edges_file is None:
        option_parser.error('Please provide at least two ' +\
                     'coordinate files or a custom edges file')

    #Open and get coord data (for multiple coords files)
    coord_files = opts.coord_fnames
    coord_files_valid = validate_coord_files(coord_files)
    if not coord_files_valid:
        option_parser.error('Every line of every coord file must ' +\
                            'have the same number of columns.')
    num_coord_files = len(coord_files)
    data['edges'], data['coord'] = \
        get_multiple_coords(coord_files, opts.edges_file, opts.serial)
    
    # if the edges file wasn't supplied, we appended _i to each file's samples
    # therefore we now add duplicated samples with _0, _1,... to mapping file
    if opts.edges_file is None:
        newmap = [data['map'][0]]
        for i in xrange(len(coord_files)):
            for sample in data['map'][1:]:
                newsample = ['%s_%d' %(sample[0],i)]
                newsample.extend(sample[1:])
                newmap.append(newsample)
        data['map'] = newmap

    # remove any samples not present in mapping file
    remove_unmapped_samples(data['map'],data['coord'],data['edges'])

    if(len(data['coord'][1]) == 0):
        raise ValueError, '\n\nError: None of the sample IDs in the coordinates files were present in the mapping file.\n'
    
    # process custom axes, if present.
    custom_axes = None
    if opts.custom_axes:
        custom_axes = process_custom_axes(opts.custom_axes)
        get_custom_coords(custom_axes, data['map'], data['coord'])
        remove_nans(data['coord'])
        scale_custom_coords(custom_axes,data['coord'])

    

    # Generate random output file name and create directories
    if opts.output_dir:
        create_dir(opts.output_dir)
        dir_path = opts.output_dir
    else:
        dir_path='./'
    
    qiime_dir=get_qiime_project_dir()

    jar_path=os.path.join(qiime_dir,'qiime/support_files/jar/')

    data_dir_path = get_random_directory_name(output_dir=dir_path,
                                              return_absolute_path=False)    

    try:
        os.mkdir(data_dir_path)
    except OSError:
        pass

    jar_dir_path = os.path.join(dir_path,'jar')
    
    try:
        os.mkdir(jar_dir_path)
    except OSError:
        pass
    
    shutil.copyfile(os.path.join(jar_path,'king.jar'), os.path.join(jar_dir_path,'king.jar'))

    filepath=coord_files[0]
    filename=filepath.strip().split('/')[-1]
    
    try:
        action = generate_3d_plots
    except NameError:
        action = None

    #Place this outside try/except so we don't mask NameError in action
    if action:
        generate_3d_plots(prefs, data, custom_axes,
               background_color, label_color,
               dir_path, data_dir_path, filename,
               ball_scale=ball_scale, arrow_colors=arrow_colors,
               user_supplied_edges=not(opts.edges_file is None))

if __name__ == "__main__":
    main()
