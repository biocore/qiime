#!/usr/bin/env python
from __future__ import division

__author__ = "Greg Caporaso"
__copyright__ = "Copyright 2010, The QIIME Project"
__credits__ = ["Greg Caporaso", "Jesse Stombaugh","Jeremy Widmann"]
__license__ = "GPL"
__version__ = "1.2.0-dev"
__maintainer__ = "Jesse Stombaugh"
__email__ = "jesse.stombaugh@colorado.edu"
__status__ = "Development"

from optparse import make_option
from qiime.util import parse_command_line_parameters, get_options_lookup
from qiime.format import build_prefs_string
from qiime.colors import get_map
from qiime.parse import parse_otu_table

options_lookup = get_options_lookup()

script_info={}
script_info['brief_description']="""Generate preferences file"""
script_info['script_description']="""This script generates a preferences (prefs) file, which can be passed to make_distance_histograms.py, make_2d_plots.py and make_3d_plots.py. The prefs file allows for defining the monte_carlo distance and gradient coloring of continuous values in the 2D and 3D plots. Currently there is only one color gradient: red to blue."""
script_info['script_usage']=[]
script_info['script_usage'].append(("""Examples:""","""To make a prefs file, the user is required to pass in a user-generated mapping file using "-m" and an output filepath, using "-o". When using the defaults, the script will use ALL categories from the mapping file, set the background to black and the monte_carlo distances to 10.""",""""make_prefs_file.py -m mapping.txt -o prefs_out.txt"""))
script_info['script_usage'].append(("""""","""If the user would like to use specified categories ('SampleID,Individual') or combinations of categories ('SampleID&&Individual'), they will need to use the -b option, where each category is comma delimited, as follows:""","""make_prefs_file.py -b "SampleID,Individual,SampleID&&Individual" -o prefs_out.txt"""))
script_info['script_usage'].append(("""""","""If the user would like to change the background color for their plots, they can pass the '-k' option, where the colors: black and white can be used for 3D plots and many additional colors can be used for the 2D plots, such as cyan, pink, yellow, etc.: ""","""make_prefs_file.py -k white -o prefs_out.txt"""))
script_info['script_usage'].append(("""""","""If the user would like to change the monte_carlo distances, they can pass the '-d' option as follows: ""","""make_prefs_file.py -d 15 -o prefs_out.txt"""))
script_info['script_usage'].append(("""""","""If the user would like to add a list of taxons they can pass the '-i' option, which is the resulting taxa file from summarize_taxa.py, as follows: ""","""make_prefs_file.py -i taxa_level_3.txt -o prefs_out.txt"""))
script_info['output_description']="""The result of this script is a text file, containing coloring preferences to be used by make_distance_histograms.py, make_2d_plots.py and make_3d_plots.py."""
script_info['optional_options']=[]


script_info['required_options']=[\
    make_option('-m', '--map_fname', dest='map_fname', \
     help='This is the metadata mapping file [default=%default]'),
    options_lookup['output_fp']
]

script_info['optional_options']=[\
    make_option('-b','--mapping_headers_to_use',action='store',\
          type='string',dest='mapping_headers_to_use',help='mapping fields to'+\
          'use in prefs file [default: %default]', default='ALL'),\
    make_option('-k', '--background_color',help='This is the background'+ \
          'color to  use in the plots. [default: %default]',default='black'),
    make_option('-d','--monte_carlo_dists',action='store',\
          type='string',dest='monte_carlo_dist',help='monte carlo distance'+\
          'to use for each sample header [default: %default]',default=10),\
    make_option('-i', '--input_taxa_file', dest='input_taxa_file',\
      action='store',type='string', help='summarized taxa file with sample' + \
            'counts by taxonomy (resulting file from summarize_taxa.py)'),
]

script_info['version']=__version__

def main():
    option_parser, opts, args = parse_command_line_parameters(**script_info)
    
    data={}
    mapping,headers,comments = get_map(opts, data)
    
    mapping_headers_to_use=opts.mapping_headers_to_use
    background_color=opts.background_color
    monte_carlo_dist=opts.monte_carlo_dist
    
    taxonomy_count_file = opts.input_taxa_file
    
    if taxonomy_count_file:
        try:
            counts_f = open(taxonomy_count_file, 'U').readlines()
            sample_ids, otu_ids, otu_table, lineages = \
                       parse_otu_table(counts_f,count_map_f=float)
        except (TypeError, IOError):
            raise ValueError, 'Summarized taxa file could not be parsed.'
    else:
        otu_ids=None
        
    out = build_prefs_string(mapping_headers_to_use, background_color, \
                                monte_carlo_dist, headers, otu_ids)
                                
    f = open(opts.output_fp,'w')
    f.write(out)
    f.close()

if __name__ == "__main__":
    main()
