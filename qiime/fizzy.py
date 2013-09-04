#!/usr/bin/env python
# File created on 02 May 2013
from __future__ import division
from cogent.app.util import ApplicationNotFoundError
import numpy


__author__ = "Gregory Ditzler"
__copyright__ = "Copyright 2011, The QIIME project"
__credits__ = ["Gregory Ditzler", "Calvin Morrison", "Gail Rosen"]
__license__ = "GPL"
__version__ = "1.7.0-dev"
__maintainer__ = "Gregory Ditzler"
__email__ = "gregory.ditzler@gmail.com"
__status__ = "Development"




def get_fs_methods():
	"""
		get_fs_methods()
		return the feature selection methods that are 
		available for use in a list. note that the options
		are case sensitive. 
	"""
	return ['CIFE','CMIM','CondMI','Condred','ICAP','JMI','MIM','MIFS','mRMR']

def parse_biom(fname): 
	"""
		parse_biom(fname)
		@fname - this is the file handle. Ex. Use something like:
			data, features = parse_biom(open('file.biom','U'))
		@data_matrix (return) - dense matrix for feature selection
		@variable_names (return) - feature names in a list
		@observation_names (return) - names of the samples in the 
			database found in fname. 
	"""
	from biom.parse import parse_biom_table
	biom_table = parse_biom_table(fname)
	observation_names = list(biom_table.SampleIds)
	variable_names = list(biom_table.ObservationIds)
	data_matrix = []
	for data in biom_table.iterObservationData():
		data_matrix.append(data)
	data_matrix = numpy.array(data_matrix)

	return data_matrix.transpose(), variable_names, observation_names


def parse_map_file(fname, column_name, observation_names):
	"""
		parse_map_file(fname)
		@fname - file handle
		@column_name - name of the column that contains the class 
			labels
		@observation_names - names of the obervatiosn in the order
			of which the samples appear in the data set. 
		@labels - numpy array of class labels
	"""
	from qiime.parse import parse_mapping_file_to_dict
	obj, comm = parse_mapping_file_to_dict(fname)
	label_full = []
	labels = []

	# grab the class labels which are likely to be in string 
	# format. 

	for id_set in observation_names:
		if(id_set not in obj):
			raise ValueError("""Unknown observation name supplied. Make sure that the
			observation name is in map file you specified""")

	for id_set in observation_names:
		if(column_name not in obj[id_set]):
			raise ValueError("""Unknown observation name supplied. Make sure that the
			observation name is in map file you specified""")

	for id_set in observation_names:
		label_full.append(obj[id_set][column_name])

	# now that we have the full class labels we need to determine
	# the number of unique classes in the the data. if the number of
	# classes is equal to the number of obervations, throw an error. 
	# its likely the user does not know what they are doing. 
	unique_classes = numpy.unique(label_full)
	if len(unique_classes) == len(observation_names):
		raise ValueError("""Number of 
			classes is the number of observations.  The number of 
			classes must be less than the number of observations in 
			map file that was specified.""")

	for str_lab in label_full:
		for uclass,n in map(None, unique_classes, range(len(unique_classes))):
			if str_lab == uclass:
				labels.append(float(n))
				break
	return numpy.array(labels)


def run_pyfeast(data, labels, features, method='MIM', n_select=15):
	"""
		run_pyfeast(data, labels, method)
		@data - numpy data (dense)
		@labels - vector of class labels (discrete)
		@method - feature selection method

		The feature selection method is based off of the FEAST 
		C variable selection toolbox. 

		Reference:
		Gavin Brown, Adam Pocock, Ming-Jie Zhao, and Mikel Lujan, 
			"Conditional Likelihood Maximisation: A Unifying Framework 
			for Information Theoretic Feature Selection," Journal of 
			Machine Learning Research, vol. 13, pp. 27--66, 2012.
			(http://jmlr.csail.mit.edu/papers/v13/brown12a.html)
	"""
	
	try:
		import feast
	except ImportError:
		raise ApplicationNotFoundError("""Error loading 
			the PyFeast module. It is likely that  
			you do not have PyFeast installed locally.""")

	try:
		fs_method = getattr(feast, method)
	except AttributeError:
		raise AttributeError("""Unknown feature selection method
			is being specified for PyFeast. Make sure the feature 
			selection method being selected is a valid one. """)

	sf = fs_method(data, labels, n_select)
	reduced_set = []
	for k in range(len(sf)):
		reduced_set.append(features[int(sf[k])])
	return reduced_set


def write_output_file(selected_features, f_out):
	"""
		write_output_file(selected_features, f_out)
		@selected_features - list of strings contaning feature names
		@f_out - output file name. this is a string not a handle!
	"""
	f = open(f_out, 'w')
	for sf in selected_features:
		f.write(sf + '\n')
	f.close()

def run_feature_selection(file_biom, file_map, column_name, out_file, method='MIM', n_select=15):
	"""
		run_feature_selection(fname_biom, fname_csv, method)
		@file_biom - handle of the biom file
		@file_map - handle of the csv file
		@out_file - result destination (string)
		@column_name - column name containing the class labels found 
			in the map file. 
		@method - feature selection method [see PyFeast docs]
		@n_select - number of features to selectio (integer)
	"""
	data_matrix, variable_names, observation_names = parse_biom(file_biom)
	label_vector = parse_map_file(file_map, column_name, observation_names)
	reduced_set = run_pyfeast(data_matrix, label_vector, variable_names, method, n_select)
	write_output_file(reduced_set, out_file)
