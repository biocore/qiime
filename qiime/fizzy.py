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


import numpy

def parse_biome(fname): 
	"""
		parse_biome(fname)
		@fname
	"""
	import json 

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

	return data_matrix.transpose(), variable_names

def parse_csv(fname):
	"""
		parse_csv(fname)
		@fname
	"""
	import csv 
	tab_data = csv.reader(fname, delimiter='\t')
	for line in tab_data:
		labels = line

	for n in range(len(labels)):
		labels[n] = float(labels[n])

	return numpy.array(labels)


def run_pyfeast(data, labels, features, method='mim', n_select=15):
	"""
		run_pyfeast(data, labels, method)
		@data
		@labels
		@method
	"""
	
	if method == "cife":
		from feast import CIFE as fs_method
	elif method == "cmim":
		from feast import CMIM as fs_method
	elif method == "condmi":
		from feast import CondMI as fs_method
	elif method == "condred":
		from feast import Condred as fs_method
	elif method == "disr":
		from feast import DISR as fs_method
	elif method == "icap":
		from feast import ICAP as fs_method
	elif method == "jmi":
		from feast import JMI as fs_method
	elif method == "mim":
		from feast import MIM as fs_method
	elif method == "mifs":
		from feast import MIFS as fs_method
	elif method == "mrmr":
		from feast import mRMR  as fs_method

	sf = fs_method(data, labels, n_select)
	reduced_set = []
	for k in range(len(sf)):
		reduced_set.append(features[int(sf[k])])
	return reduced_set


def write_output_file(selected_features, f_out):
	"""
		write_output_file(selected_features, f_out)
		@selected_features
		@f_out
	"""
	f = open(f_out, 'w')
	for sf in selected_features:
		f.write(sf + '\n')
	f.close()

def run_feature_selection(file_biome, file_csv, out_file, method='mim', n_select=15):
	"""
		run_feature_selection(fname_biome, fname_csv, method)
		@fname_biome
		@fname_csv 
		@method - feature selection method [see PyFeast docs]
		@fmt - input file format (['biome', 'table'])
	"""
	data_matrix, variable_names = parse_biome(file_biome)
	label_vector = parse_csv(file_csv)
	reduced_set = run_pyfeast(data_matrix, label_vector, variable_names, method, n_select)
	write_output_file(reduced_set, out_file)
	# thats all folks

