#!/usr/bin/env python
# File created on 02 May 2013
from __future__ import division

__author__ = "Gregory Ditzler"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Gregory Ditzler", "Calvin Morrison"]
__license__ = "GPL"
__version__ = "1.3.0"
__maintainer__ = "Gregory Ditzler"
__email__ = "gregory.ditzler@gmail.com"
__status__ = "Release"
 


from qiime.util import parse_command_line_parameters, make_option
import qiime.fizzy as fizzy 

script_info = {}
script_info['brief_description'] = ""
script_info['script_description'] = ""
script_info['script_usage'] = [("","","")]
script_info['output_description']= ""
script_info['required_options'] = [
	# Example required option
	make_option('-i','--input_file',type="existing_file",help='input sparse biome file'),\
	make_option('-c','--label_file',type="existing_file",help='csv file with labeling scheme')
]
script_info['optional_options'] = [
	# Example optional option
	make_option('-o','--output_file',type="new_dirpath",help='the output directory [default: %default]'),\
	make_option('-f','--fs_method',type="alg_parameter",help='feature selection method [options: jmi, mrmr, mim, mifs]')
]
script_info['version'] = __version__



def main():
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
	if not opts.output_file:
		out_fmt = 'output.txt'
	else:
		out_fmt = opts.output_file

	# run the fizzy feature selection routine
	fizzy.run_feature_selection( \
		open(opts.input_file,'U'), \
		open(opts.label_file,'U'), \
		fs_method, out_fmt)


if __name__ == "__main__":
	main()