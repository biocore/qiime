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


def parse_biome(fname): 
	"""
		parse_biome(fname)
		@fname
	"""
	import json 
	import numpy 

	json_string = fname.read()
	obj = json.loads(json_string)

	# the biom file will have the names of the features and the 
	# otu counts that we need to use with our algorithm. 
	variable_names = []
	variables = obj['rows']
	for var in variables:
		variable_names.append( var['id'] )

	# extract the data from the json structure. we are assuming that
	# the file in a sparse format. 
	data_matrix = numpy.zeros(obj['shape'])
	for index_set in obj['data']:
		data_matrix[index_set[0], index_set[1]] = index_set[2]

	return data_matrix, variable_names

def parse_csv(fname):
	"""
		parse_csv(fname)
		@fname
	"""
	return None

def run_feature_selection(file_biome, file_csv, method='mim', fmt='biome'):
	"""
		run_feature_selection(fname_biome, fname_csv, method)
		@fname_biome
		@fname_csv 
		@method - feature selection method [see PyFeast docs]
		@fmt - input file format (['biome', 'table'])
	"""
	data_matrix, variable_names = parse_biome(file_biome)
	parse_csv(file_csv)


