#!/usr/bin/env python
# File created on 02 May 2013
from __future__ import division

__author__ = "Gregory Ditzler"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Gregory Ditzler","Calvin Morrison","Gail Rosen"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Gregory Ditzler"
__email__ = "gregory.ditzler@gmail.com"
__status__ = "Development"
 


from qiime.util import parse_command_line_parameters, make_option
import qiime.fizzy as fizzy 


feature_selection_choices = ['cife', 'condmi','cmim','condred','icap','jmi','mim','mifs','mrmr']

script_info = {}
script_info['brief_description'] = """Run feature selection on aubundance data \
	contained in a Biom file."""
script_info['script_description'] ="""This script will run a feature selection 
	algorithm on abundance data contained in a Biom file given a label file, which \
	is stored in tab-delimited format. The current feature selection methods uses \
	a forward search algorithm to select the features. The objective functions are \
	based on information theory. At the moment, users are limited to the objective \
	functions implemented in the PyFeast feature selection module. """
script_info['script_usage'] = [(\
	"""Run JMI feature selection on a Biom file:""",\
	"""To perform feature selection the biom file and TSV file must be specified \
	in advance. Here we use JMI and select 15 features. """,\
	"""%prog -i data.biom -m label.map -c Class_Column -f jmi -k 15""")]
script_info['output_description']= """Text file containing the top features \
	selected by the algorithm. """

script_info['required_options'] = [
	# Example required option
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
	# Example optional option
	make_option('-o',\
		'--output_path',\
		type="new_dirpath",\
		help='the output directory [default: %default]'),\
	make_option('-k',\
		'--n_select',\
		type="int",\
		help='number of feature to select [default: 15]'),\
	#make_option('-f',\
	#	'--fs_method',\
	#	type="string",\
	#	help='feature selection method [options: jmi, mrmr, mim, mifs]')
	make_option('-f',\
		'--fs_method', \
		type='choice', \
		help='feature selection method. valid options are ' \
		+ ', '.join(feature_selection_choices) + '. [default: %default]',
		choices=feature_selection_choices,
		default='mim')
]
script_info['version'] = __version__



def main():
	"""
		main()
	"""
	option_parser, opts, args = parse_command_line_parameters(**script_info)

	# check to see if the user has set the feature selection. if not set
	# then use the mutual information maximization algorithm. 
	if not opts.fs_method:
		fs_method = 'mim'
	else:
		fs_method = opts.fs_method

	# check to see if the user has see the output file. if the output file
	# has not been set then we will make a generic file with a "meaningless"
	# file name
	if not opts.output_path:
		out_fmt = 'output.txt'
	else:
		out_fmt = opts.output_path

	if not opts.n_select:
		n_select = 15
	else:
		n_select = int(opts.n_select)

	# run the fizzy feature selection routine
	fizzy.run_feature_selection( \
		open(opts.input_path,'U'), \
		open(opts.map_path,'U'), \
		opts.column_label, \
		out_fmt, \
		fs_method, \
		n_select)


if __name__ == "__main__":
	main()

