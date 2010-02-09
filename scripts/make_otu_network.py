#!/usr/bin/env python
# File created on 09 Feb 2010
from __future__ import division

__author__ = "Julia Goodrich"
__copyright__ = "Copyright 2010, The QIIME project"
__credits__ = ["Julia Goodrich"]
__license__ = "GPL"
__version__ = "1.0-dev"
__maintainer__ = "Julia Goodrich"
__email__ = "julia.goodrich@colorado.edu"
__status__ = "Pre-release"

"""
This script generates the otu networks and statistics

"""

from qiime.util import parse_command_line_parameters
from optparse import make_option
from qiime.make_otu_network import create_dir, create_network_and_stats


script_description = """This script generates the otu network files to be \
passed into cytoscape and statistics for those networks"""

script_usage = """Create network cytoscape and statistic files in folder \
otu_network within user specified output folder(-o). It uses the OTU file \
otus.txt (-i) and the user mapping file input_map.txt (-m). Output will be \
edge and node files to be loaded into cytoscape and props files labeled by \
category used for coloring. All of the files will be in folder otu_network. \
By default this folder will be put into the current working dir (overwrite \
with -o).

python ~/code/Qiime/trunk/qiime/make_otu_network.py -i otus.txt -m \
input_map.txt -o /Users/bob/qiime_run/
"""

required_options = [\
make_option('-m', '--mapping_file',action='store',type='string',\
				dest='map_file',\
                help='name of input map file [REQUIRED]'),
make_option('-i', '--input_file',action='store',type='string',\
			dest='counts_file',\
            help='name of otu table file [REQUIRED]')
]

optional_options = [\
make_option('-o', '--output_dir', action='store',type='string',\
			   dest='dir_path',\
               help='output directory for all analyses [default: cwd]')
]


def main():
	option_parser, opts, args = parse_command_line_parameters(
      script_description=script_description,
      script_usage=script_usage,
      version=__version__,
      required_options=required_options,
      optional_options=optional_options)
	if not opts.counts_file:
		parser.error("An otu table file must be specified")

	if not opts.map_file:
		parser.error("A Map file must be specified")
	dir_path = create_dir(opts.dir_path)
	map_lines = open(opts.map_file,'U').readlines()
	otu_sample_lines = open(opts.counts_file, 'U').readlines()
	create_network_and_stats(dir_path,map_lines,otu_sample_lines)

if __name__ == "__main__":
    main()
