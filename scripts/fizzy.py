#!/usr/bin/env python
# File created on 02 May 2013
from __future__ import division
from qiime.util import parse_command_line_parameters, make_option
import qiime.fizzy as fizzy 

__author__ = "Gregory Ditzler"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Gregory Ditzler","Calvin Morrison","Gail Rosen"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Gregory Ditzler"
__email__ = "gregory.ditzler@gmail.com"
__status__ = "Development"
 

script_info = {}
script_info['brief_description'] = """Run feature selection on 
  aubundance data contained in a Biom file."""
script_info['script_description'] ="""This script will run a 
	feature selection algorithm on abundance data contained in a 
	Biom file given a label file, which is stored in tab-delimited 
	format. The current feature selection methods uses a forward 
	search algorithm to select the features. The objective functions 
	are based on information theory. At the moment, users are limited 
	to the objective functions implemented in the PyFeast feature 
	selection module. """
script_info['script_usage'] = [(
	"""Run JMI feature selection on a Biom file:""",
	"""To perform feature selection the biom file, map file must 
	be specified in advance, and the label column in the map file. 
	Here we use JMI and select 15 features. """,
	"""%prog -i data.biom -m map.txt -c Class_Column -f jmi -k 15""")]
script_info['output_description']= """Text file containing the top features 
	selected by the algorithm. """

# set the required options. for fizzy we need to have an exisiting input 
# file, map file, and label column in the map file specified. other than
# that we can set some default parameters to run a feature selection 
# algorithm.
script_info['required_options'] = [
	make_option('-i', \
		'--input_path',\
		type="existing_filepath",
		help='input sparse biome file'),\
	make_option('-m',\
		'--map_path',\
		type="existing_filepath",\
		help='map file with labeling scheme'),
	make_option('-c',\
		'--column_label',\
		type="string",\
		help='column indicating the labels in the map file.')
]
script_info['optional_options'] = [
	make_option('-o',\
		'--output_path',\
		type="new_dirpath",\
		help='the output directory [default: %default]',\
		default='output.txt'),\
	make_option('-k',\
		'--n_select',\
		type="int",\
		help='number of feature to select [default: %default]',\
		default=15),\
	make_option('-f',\
		'--fs_method', \
		type='choice', \
		help='feature selection method. valid options are ' \
		+ ', '.join(fizzy.get_fs_methods()) + '. [default: %default]',
		choices=fizzy.get_fs_methods(),
		default='MIM')
]
script_info['version'] = __version__



def main():
	"""
		main()
	"""
	option_parser, opts, args = parse_command_line_parameters(**script_info)

	# run the fizzy feature selection routine
	fizzy.run_feature_selection( \
		open(opts.input_path,'U'), \
		open(opts.map_path,'U'), \
		opts.column_label, \
		opts.output_path, \
		opts.fs_method, \
		opts.n_select)


if __name__ == "__main__":
	main()

